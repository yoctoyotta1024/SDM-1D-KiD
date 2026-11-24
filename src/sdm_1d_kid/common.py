"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: common.py
Project: sdm_1d_kid
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
common functions between KiD setup of PySDM (and its PyMPDATA through PySDMKiDDynamics)
and other microphysics schemes (and through their PyMPDATA KiDDynamics)
"""


def get_options(is_accurate=False):
    from PyMPDATA import Options

    if is_accurate:
        return Options(
            n_iters=3,
            infinite_gauge=True,
            nonoscillatory=True,
            third_order_terms=True,
        )
    else:
        return Options(n_iters=3, nonoscillatory=True)


def qvap_initial_condition(particle_reservoir_depth=None):
    from scipy.interpolate import interp1d

    if particle_reservoir_depth is not None:
        return interp1d(
            (-max(particle_reservoir_depth, 1), 0, 740, 3260),
            (0.015, 0.015, 0.0138, 0.0024),
            fill_value="extrapolate",
        )
    else:
        return interp1d(
            (0, 740, 3260), (0.015, 0.0138, 0.0024), fill_value="extrapolate"
        )


def theta_initial_condition():
    from scipy.interpolate import interp1d

    return interp1d(
        (0.0, 740.0, 3260.00),
        (297.9, 297.9, 312.66),
        fill_value="extrapolate",
    )


def hydrostatic_equilibrium(temperature, pressure, thd, qv0, is_approx_drhod_dz):
    """
    returns function to integrate (with solve_ivp) to get hydrostatic equilibrium
    for initial thermodynamic profiles

    if ``is_approx_drhod_dz==True``, then effect of water vapour on density/pressure
    profile is neglected.
    """
    from numdifftools import Derivative
    from PyMPDATA_examples.Shipway_and_Hill_2012 import formulae
    from PyMPDATA_examples.Shipway_and_Hill_2012.formulae import const

    def drhod_dz(z, rhod):
        if z < 0:
            return 0
        T = temperature(rhod[0], thd(z, qv0(z)))
        p = pressure(rhod[0], T, qv0(z))
        drhod_dz = formulae.drho_dz(const.g, p, T, qv0(z), const.lv)
        if not is_approx_drhod_dz:  # to resolve issue #335
            qv = qv0(z)
            dqv_dz = Derivative(qv0)(z)
            drhod_dz = drhod_dz / (1 + qv) - rhod * dqv_dz / (1 + qv) ** 2
        return drhod_dz

    return drhod_dz


def exner_function_thermodynamics(
    p_surf, _th, is_exner_novapour_uniformrho, z_points=None
):
    """
    Returns function for initial thermodynamics profiles as in
    Shipway and Hill (2012) and Hill et al. (2023), except for pressure
    where we account for vapour to match PySDM (see below).

    If ``is_exner_novapour_uniformrho==False``, profiles calculated as in
    original KiD setup of Shipway and Hill (2012).

    If ``is_exner_novapour_uniformrho==True``, over-writes density profile
    to be a constant = 1.0 throughout the column, meanwhile the other
    thermodynamics is still calculated using Exner function. This is what was
    done by Hill et al. (2023).

    *Note:* these formulae/choices are not documented well in their papers, see:
    https://github.com/BShipway/KiD/blob/bad81aa6efa4b7e4743b6a1867382fc74c10a884/src/test_cases.f90#L784
    """
    import numpy as np
    from PyMPDATA_examples.Shipway_and_Hill_2012.formulae import const
    from scipy.interpolate import interp1d

    def z2exner(z, return_rho=False, return_press=False, return_temp=False):
        # z2exner function requires z >= 0.0
        # therefore where z < 0.0, treat as if z == 0.0
        z = np.where(z < 0.0, 0.0, z)

        nz = len(z)
        p_ref = 100000.0
        r_on_cp = const.Rd / const.c_pd
        theta = _th(z)

        dexner = np.zeros(nz)
        dexner[0] = const.g * z[0] / (const.c_pd * theta[0])
        for k in range(1, nz):
            dexner[k] = (
                const.g
                * (z[k] - z[k - 1])
                / (const.c_pd * 0.5 * (theta[k] + theta[k - 1]))
            )

        exner = np.zeros(nz)
        exner[-1] = (p_surf / p_ref) ** r_on_cp - np.sum(dexner)
        for k in range(nz - 2, -1, -1):
            exner[k] = exner[k + 1] + dexner[k]

        if return_rho:
            rho = (p_ref * exner ** (1.0 / r_on_cp - 1.0)) / (const.Rd * theta)
            return rho
        if return_press:
            return p_ref * exner ** (1.0 / r_on_cp)
        if return_temp:
            return theta * exner
        else:
            return exner

    def temp(z):
        return z2exner(z, return_temp=True)

    def press(z, qv):
        # Exner gives dry-air pressure; to be consistent with PySDM
        # (see state_variable_triplet.p), we add the vapour partial pressure here
        # (q_v = eps * p_v / p_d) to calculate the total pressure.
        # *Note:* this is different to the original KiD setup which uses the
        # dry air pressure (z2exner(z, return_press=True)) as the total pressure.
        return z2exner(z, return_press=True) * (1 + qv / const.eps)

    if is_exner_novapour_uniformrho:
        rhod0 = 1.0 * np.ones_like(z_points)  # [kg / m^3]

        def rhod(z):
            z = np.where(z < 0.0, 0.0, z)
            return interp1d(z_points, rhod0)(z)
    else:

        def rhod(z):
            return z2exner(z, return_rho=True)

    return temp, press, rhod


def get_initial_spectrum_and_size_range(norm_factor, m_mode, s_geom):
    from PySDM.initialisation import spectra

    spectrum = spectra.Lognormal(norm_factor=norm_factor, m_mode=m_mode, s_geom=s_geom)

    ### *(!)Note:* if changing default_cdf_range / size_range, be mindful that
    ### MINRADIUS/MAXRADIUS in Cleo's config.yaml may also need to change
    default_cdf_range = (0.00001, 0.99999)
    size_range = spectrum.percentiles(
        default_cdf_range
    )  # min and max range of radii to sample [m]

    return spectrum, size_range


def get_alpha_sampling(probdistrib, alpha, size_range=None):
    import numpy as np
    from PySDM.initialisation.sampling.spectral_sampling import AlphaSampling

    return AlphaSampling(
        probdistrib,
        alpha=alpha,
        size_range=size_range,
        interp_points=100000,
        dist_1_inv=lambda y, size_range: np.exp(
            (np.log(size_range[1]) - np.log(size_range[0])) * y + np.log(size_range[0])
        ),
    )


def get_spectral_sampling_method(alpha_sampling):
    return alpha_sampling.sample_quasirandom
