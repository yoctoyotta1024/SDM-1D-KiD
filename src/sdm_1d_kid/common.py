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
