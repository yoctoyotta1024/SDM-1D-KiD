"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: settings.py
Project: pympdata_driver
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
File is adaptation from PyMPDATA Shipway and Hill 2012 example, copied from:
https://github.com/open-atmos/PyMPDATA/blob/main/examples/PyMPDATA_examples/Shipway_and_Hill_2012/settings.py
between v1.6.1 and v1.6.2
"""

import numpy as np
from PyMPDATA_examples.Shipway_and_Hill_2012 import formulae
from PyMPDATA_examples.Shipway_and_Hill_2012.formulae import si
from pystrict import strict
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

import sdm_1d_kid.common as cmn
from sdm_1d_kid.pympdata_driver.arakawa_c import arakawa_c


@strict
class Settings:
    def __init__(
        self,
        dt: float,
        dz: float,
        wmax_const: float = 3 * si.m / si.s * si.kg / si.m**3,
        tscale_const: float = 600 * si.s,
        t_max: float = 15 * si.minutes,
        z_min: float = 0 * si.metres,
        z_max: float = 3000 * si.metres,
        p_surf: float = 1000 * si.hPa,
        is_exner_novapour: bool | None = False,
        is_exner_novapour_uniformrho: bool | None = False,
        is_approx_drhod_dz: bool | None = False,
    ):
        ### General Settings
        self.dt = dt
        self.dz = dz

        self.z_min = z_min
        self.z_max = z_max
        self.t_max = t_max

        ### Set Time-Varying Updraught
        rhod_w_const = wmax_const
        t_1 = tscale_const
        self.rhod_w = lambda t: rhod_w_const * np.sin(np.pi * t / t_1) if t < t_1 else 0

        ### Set (Initial) Thermodynamic Profile
        self.qv0 = cmn.qvap_initial_condition(particle_reservoir_depth=-self.z_min)
        self._th = cmn.theta_initial_condition()

        self.thd = lambda z, qv: formulae.th_dry(self._th(z), qv)

        ### z position and function to deal with reservoir below surface
        z_points = np.arange(0, self.z_max + self.dz / 2, self.dz / 2)

        ### function to deal with z < 0 like it is z == 0
        ### so returns z height either above or at surface [m]
        zpositive = lambda z: np.where(z < 0.0, 0.0, z)

        if is_exner_novapour or is_exner_novapour_uniformrho:
            temp, press, rhod = cmn.exner_function_thermodynamics(
                p_surf, self._th, is_exner_novapour_uniformrho, z_points=z_points
            )
            self.temp = temp
            self.press = press
            self.rhod = rhod
        else:
            drhod_dz = cmn.hydrostatic_equilibrium(
                formulae.temperature,
                formulae.pressure,
                self.thd,
                self.qv0,
                is_approx_drhod_dz,
            )

            rhod0 = formulae.rho_d(p_surf, self.qv0(0), self._th(0))
            rhod_solution = solve_ivp(
                fun=drhod_dz,
                t_span=(0, self.z_max),
                y0=np.asarray((rhod0,)),
                t_eval=z_points,
                max_step=self.dz / 2,
            )
            assert rhod_solution.success

            self.rhod = lambda z: interp1d(z_points, rhod_solution.y[0])(zpositive(z))
            self.temp = lambda z: formulae.temperature(
                self.rhod(zpositive(z)), self.thd(zpositive(z), self.qv0(zpositive(z)))
            )
            self.press = lambda z, qv: formulae.pressure(
                self.rhod(zpositive(z)), self.temp(zpositive(z)), qv
            )

        self.nz_vec = arakawa_c.z_vector_coord((self.nz,))
        self.press0 = lambda z: self.press(z, self.qv0(zpositive(z)))

    @property
    def nz(self):
        nz = (self.z_max - self.z_min) / self.dz
        assert nz == int(nz)
        return int(nz)

    @property
    def nt(self):
        nt = self.t_max / self.dt
        assert nt == int(nt)
        return int(nt)
