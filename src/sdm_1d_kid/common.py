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


# %% -------------------------- KiD Driver Setup ----------------------------- #
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


# %% ------------------- Initial Superdroplet Conditions --------------------- #
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


# %% -------------------------- Post-Processing ------------------------------ #
### --- MetPy calculations on numpy blocks
# (applied lazily with xarray.apply_ufunc) -- ###
def _lazy(func, *arrays):
    """Apply numpy function to (dask-backed) DataArrays block-wise without loading them."""
    import xarray as xr

    return xr.apply_ufunc(func, *arrays, dask="parallelized", output_dtypes=[float])


def _vapour_pressure_hpa(press_hpa, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.vapor_pressure(press_hpa * units.hPa, qvap_kgkg * units("kg/kg"))
        .to("hPa")
        .magnitude
    )


def _theta_k(press_hpa, temp_k):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.potential_temperature(press_hpa * units.hPa, temp_k * units.kelvin)
        .to("K")
        .magnitude
    )


def _theta_virtual_k(press_hpa, temp_k, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.virtual_potential_temperature(
            press_hpa * units.hPa, temp_k * units.kelvin, qvap_kgkg * units("kg/kg")
        )
        .to("K")
        .magnitude
    )


def _density_kgm3(press_hpa, temp_k, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.density(
            press_hpa * units.hPa, temp_k * units.kelvin, qvap_kgkg * units("kg/kg")
        )
        .to("kg/m^3")
        .magnitude
    )


def _relative_humidity_percent(press_hpa, temp_k, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.relative_humidity_from_mixing_ratio(
            press_hpa * units.hPa, temp_k * units.kelvin, qvap_kgkg * units("kg/kg")
        )
        .to("percent")
        .magnitude
    )


def add_thermodynamics(ds):
    """Add vapour and dry-air pressures, densities, potential temperatures and relative humidity.

    Densities and virtual potential temperature include water vapour but not condensate.
    The dry-air potential temperature uses the partial pressure of dry air.

    Args:
        ds (xarray.Dataset): Dataset with CLEO thermodynamics ("press" [hPa], "temp" [K]
          and "qvap" [g/kg]).
    """
    qvap = ds.qvap / 1000  # [kg/kg]
    press_vapour = _lazy(_vapour_pressure_hpa, ds.press, qvap)
    press_dry = ds.press - press_vapour  # [hPa]
    rho = _lazy(_density_kgm3, ds.press, ds.temp, qvap)

    return ds.assign(
        relh=_lazy(_relative_humidity_percent, ds.press, ds.temp, qvap).assign_attrs(
            units="%", long_name="relative humidity"
        ),
        theta=_lazy(_theta_k, ds.press, ds.temp).assign_attrs(
            units="K", long_name="potential temperature"
        ),
        theta_dry=_lazy(_theta_k, press_dry, ds.temp).assign_attrs(
            units="K",
            long_name="dry air potential temperature",
            description="(potential temperature from partial pressure of dry air)",
        ),
        theta_virtual=_lazy(_theta_virtual_k, ds.press, ds.temp, qvap).assign_attrs(
            units="K", long_name="virtual potential temperature"
        ),
        rho=rho.assign_attrs(
            units="kg m^-3", long_name="density of moist air (dry air + vapour)"
        ),
        rho_dry=(rho / (1 + qvap)).assign_attrs(
            units="kg m^-3", long_name="density of dry air"
        ),
        press_vapour=press_vapour.assign_attrs(
            units="hPa", long_name="partial pressure of water vapour"
        ),
        press_dry=press_dry.assign_attrs(
            units="hPa", long_name="partial pressure of dry air"
        ),
    )


### --- Generic post-processing functions --- ###
def order_dimensions(ds):
    """Order dimensions of every variable as (ensemble, height, time, time_sds, radius_bin)."""
    return ds.transpose(
        "ensemble",
        "height",
        "time",
        "time_sds",
        "radius_bin",
        ...,
        missing_dims="ignore",
    )


def metadata_corrections(ds):
    """Correct non-standard "units" attributes of all variables and coordinates.

    - "" -> "1": CF convention for dimensionless quantities
    - "micro-m" -> "um": parseable by UDUNITS and pint
    - "g/Kg" -> "g/kg": typo in CLEO's output (v0.71.0, libs/observers/thermo_observer.hpp)

    Only attributes are changed (the data are not loaded) and the input dataset is not modified.
    """
    corrections = {"": "1", "micro-m": "um", "g/Kg": "g/kg"}

    ds = ds.copy()  # shallow copy: new attributes dictionaries, same (lazy) data
    for var in ds.variables.values():
        units = var.attrs.get("units")
        if units in corrections:
            var.attrs["units"] = corrections[units]
    return ds


### --- Checks on final post-processed datasets --- ###
def generic_check_variable(ds, logger, var, var2, source, source2, rtol):
    """Log whether 'var2' matches the reference 'var' to within relative tolerance 'rtol'."""
    import numpy as np

    reldiff = float((np.abs(ds[var] - ds[var2]) / ds[var]).max())
    if reldiff > rtol:
        logger.warning(
            f"{var} ({source}) differs from {var2} ({source2}) by up to {reldiff:.2e} "
            f"(relative to {var})"
        )
        return False
    logger.info(f"{var} ({source}) matches {var2} ({source2}) ({reldiff:.2e})")
    return True


def check_dry_density(ds, logger, source="MetPy", source2="KiD Driver", rtol=1e-3):
    """Log whether the dry air density recomputed by MetPy ("rho_dry") matches
    the density prescribed by a KiD driver ("rhod").

    The two are calculated differently (the driver's formulae and constants vs. MetPy's
    from p, T and qvap), so they need not agree exactly; a large difference indicates an
    inconsistency between the driver's thermodynamics and the post-processing (e.g. units).
    """
    return generic_check_variable(ds, logger, "rho_dry", "rhod", source, source2, rtol)


def check_relative_humidity(
    ds, logger, source="MetPy", source2="KiD Driver", rtol=5e-3
):
    """Same as check_dry_density but for relative humidity"""
    return generic_check_variable(ds, logger, "relh", "rh", source, source2, rtol)


def check_dry_potential_temperature(
    ds, logger, source="MetPy", source2="KiD Driver", rtol=1e-3
):
    """Same as check_dry_density but for dry-air potential temperature"""
    return generic_check_variable(ds, logger, "theta_dry", "thd", source, source2, rtol)
