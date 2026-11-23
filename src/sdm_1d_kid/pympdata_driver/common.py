"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: common.py
Project: pympdata_driver
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
common functions between PySDM KiDDynamics and PyMPDATA KiDDynamics
both should be the same to run the same KiD
"""


def get_options(is_pysdm=False):
    from PyMPDATA import Options

    if is_pysdm:
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
