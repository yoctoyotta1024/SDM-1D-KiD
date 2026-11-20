"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: formulae.py
Project: thermo
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
functions for calculations of some quantities e.g. potential temperature(s)
"""

import numpy as np


def dry_potential_temperature(temp: np.ndarray, press: np.ndarray):
    r"""Calculate the potential temperature for dry air.

    This function calculates the potential temperature for dry air given the temperature and pressure.

    .. math::
          \theta_{\rm{dry}} = T \left( \frac{P_{\rm ref}}{P} \right)
              ^{ \frac{R_{\rm{dry}}}{c_{\rm{p, dry}}} }

    where :math:`P_{\rm ref}` is the first pressure in press array.

    Args:
      temp (array-like):
          Temperature values (K).
      press (array-like):
          Pressure values (Pa).

    Returns:
        array-like: The dry potential temperature (K).
    """
    from metpy import calc
    from metpy.units import units

    theta_dry = calc.potential_temperature(press * units.Pa, temp * units.kelvin)

    return theta_dry.magnitude  # Kelvin


def moist_equiv_potential_temperature(
    temp: np.ndarray, press: np.ndarray, qvap: np.ndarray
):
    """
    Calculate the moist potential temperature.

    .. math::
          \theta_e = \theta \\cdot \\exp\\left(\frac{L_v \\cdot q}{c_p \\cdot T}\right)

    Args:
        temp (array-like):
            Temperature values (K).
        press (array-like):
            Pressure values (Pa).
        qvap (array-like):
            Mixing ratio of water vapour (kg/kg)
    Returns:
        array-like: The moist potential temperature (K).
    """
    from metpy import calc
    from metpy.units import units

    relh = calc.relative_humidity_from_mixing_ratio(
        press * units.Pa, temp * units.kelvin, qvap
    )
    dewpoint = calc.dewpoint_from_relative_humidity(temp * units.kelvin, relh)

    theta_equiv = calc.equivalent_potential_temperature(
        press * units.Pa, temp * units.kelvin, dewpoint
    )

    return theta_equiv.magnitude  # Kelvin


def moist_static_energy(
    temp: np.ndarray, qvap: np.ndarray, height: np.ndarray | None = None
):
    """
    Calculate the moist static energy [kilojoule / kilogram]

    .. math::
          \theta_e = \theta \\cdot \\exp\\left(\frac{L_v \\cdot q}{c_p \\cdot T}\right)

    Args:
        temp (array-like):
            Temperature values (K).
        press (array-like):
            Pressure values (Pa).
        qvap (array-like):
            Mixing ratio of water vapour (kg/kg)
    Returns:
        array-like: The moist potential temperature (K).
    """
    from metpy import calc
    from metpy.units import units

    if height is None:
        height = np.zeros(temp.shape)

    specific_humidity = calc.specific_humidity_from_mixing_ratio(qvap)
    mse = calc.moist_static_energy(
        height * units.meters, temp * units.kelvin, specific_humidity
    )

    return mse.magnitude  # [kilojoule / kilogram]


def supersaturation(temp: np.ndarray, press: np.ndarray, qvap: np.ndarray):
    """
    Calculate supersaturation based on the method described in PyMPDATA-examples

    This function uses the calculations in the Shipway and Hill (2012) example from
    PyMPDATA-examples library to compute the supersaturation.

    Parameters:
    temp (float): Temperature in Kelvin.
    press (float): Pressure in Pascals.
    qvap (float): Specific humidity (kg/kg).

    Returns:
    float: Supersaturation value.
    """
    from PyMPDATA_examples import Shipway_and_Hill_2012 as kid

    pvs = kid.formulae.pvs_Celsius(temp - kid.const.T0)
    relh = kid.formulae.pv(press, qvap) / pvs

    return relh - 1
