"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: settings.py
Project: helpers
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
File is adaptation from Shipway and Hill 2012 example from PySDM repo, copied from:
https://github.com/open-atmos/PySDM/blob/main/examples/PySDM_examples/featured/Shipway_and_Hill_2012/settings.py
between v3.1.0 and v3.1.1
"""

import numpy as np
from PySDM import Formulae
from PySDM.dynamics import condensation
from PySDM.dynamics.collisions.collision_kernels import Geometric
from PySDM.physics import si
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

import sdm_1d_kid.pympdata_driver.common as cmn


class Settings:
    def __init__(
        self,
        dt: float,
        dz: float,
        wmax_const: float = 3 * si.m / si.s * si.kg / si.m**3,
        tscale_const: float = 600 * si.s,
        t_max: float = 15 * si.minutes,
        z_max: float = 3000 * si.m,
        p_surf: float = 1000 * si.hPa,
        is_approx_drhod_dz: bool = False,
        n_sd_per_gridbox: int = 128,
        particles_per_volume_STP: int = 50 / si.cm**3,
        m_mode=0.08 / 2 * si.um,
        s_geom=1.4,
        kappa: float = 0.9,
        z_part: tuple | None = None,
        enable_condensation: bool = True,
        precip: bool = True,
        formulae: Formulae = None,
        collision_kernel=None,
        save_spec_and_attr_times=(),
        cloud_water_radius_range=(1 * si.um, 50 * si.um),
        rain_water_radius_range=(50 * si.um, np.inf * si.um),
    ):
        ### General Settings
        self.dt = dt
        self.dz = dz

        self.z_max = z_max
        self.t_max = t_max

        ### General PySDM SDM Settings
        self.formulae = formulae or Formulae()
        self.n_sd_per_gridbox = n_sd_per_gridbox
        self.kappa = kappa
        self.particles_per_volume_STP = particles_per_volume_STP
        self.precip = precip
        self.enable_condensation = enable_condensation
        self.z_part = z_part
        self.collision_kernel = collision_kernel or Geometric(collection_efficiency=1)
        self.cloud_water_radius_range = cloud_water_radius_range
        self.rain_water_radius_range = rain_water_radius_range
        self.save_spec_and_attr_times = save_spec_and_attr_times

        ### Set Time-Varying Updraught
        rhod_w_const = wmax_const
        t_1 = tscale_const
        self.rhod_w = lambda t: rhod_w_const * np.sin(np.pi * t / t_1) if t < t_1 else 0

        ### Initial PySDM SDM Superdroplets Conditions
        apprx_w1 = rhod_w_const / self.formulae.constants.rho_STP
        self.particle_reservoir_depth = (
            (2 * apprx_w1 * t_1 / np.pi) // self.dz + 1
        ) * self.dz
        a, b = cmn.get_initial_spectrum_and_size_range(
            norm_factor=particles_per_volume_STP / self.formulae.constants.rho_STP,
            m_mode=m_mode,
            s_geom=s_geom,
        )
        self.spectrum = a  # wet_radius_spectrum_per_mass_of_dry_air
        self.size_range = b  # [min, max] radii to sample

        ### Set (Initial) Thermodynamic Profile
        self.qv0 = cmn.qvap_initial_condition(
            particle_reservoir_depth=self.particle_reservoir_depth
        )
        self._th = cmn.theta_initial_condition()

        def thd(z, qv=None):
            if qv is not None:
                return self.formulae.state_variable_triplet.th_dry(self._th(z), qv)
            else:
                return self.formulae.state_variable_triplet.th_dry(
                    self._th(z), self.qv0(z)
                )

        self.thd = thd

        ### z position and function to deal with reservoir below surface
        self.z_min = -self.particle_reservoir_depth
        z_points = np.linspace(self.z_min, self.z_max, 2 * self.nz + 1)

        ### function to convert dimensionless vertical coordinate into
        ### height above reservoir, i.e. relative to surface [m]
        self.zposition = lambda zZ: zZ * (self.nz * self.dz) + self.z_min

        drhod_dz = cmn.hydrostatic_equilibrium(
            self.formulae.state_variable_triplet.T,
            self.formulae.state_variable_triplet.p,
            self.thd,
            self.qv0,
            is_approx_drhod_dz,
        )

        self.rhod0 = self.formulae.state_variable_triplet.rho_d(
            p=p_surf,
            water_vapour_mixing_ratio=self.qv0(0 * si.m),
            theta_std=self._th(0 * si.m),
        )
        rhod_solution = solve_ivp(
            fun=drhod_dz,
            t_span=(self.z_min, self.z_max),
            y0=np.asarray((self.rhod0,)),
            t_eval=z_points,
            max_step=dz / 2,
        )
        assert rhod_solution.success

        self.rhod = interp1d(z_points, rhod_solution.y[0])

        ### PySDM SDM Settings
        self.condensation_rtol_x = condensation.DEFAULTS.rtol_x
        self.condensation_rtol_thd = condensation.DEFAULTS.rtol_thd
        self.condensation_adaptive = True
        self.condensation_update_thd = False
        self.coalescence_adaptive = True

        ### PySDM Output Settings
        self.number_of_bins = 100
        self.r_bins_edges_dry = np.logspace(
            np.log10(0.001 * si.um),
            np.log10(1 * si.um),
            self.number_of_bins + 1,
            endpoint=True,
        )
        self.r_bins_edges = np.logspace(
            np.log10(0.001 * si.um),
            np.log10(10 * si.mm),
            self.number_of_bins + 1,
            endpoint=True,
        )

    @property
    def n_sd(self):
        return self.nz * self.n_sd_per_gridbox

    @property
    def nz(self):
        assert (
            self.particle_reservoir_depth / self.dz
            == self.particle_reservoir_depth // self.dz
        )
        nz = (self.z_max + self.particle_reservoir_depth) / self.dz
        assert nz == int(nz)
        return int(nz)

    @property
    def nt(self):
        nt = self.t_max / self.dt
        assert nt == int(nt)
        return int(nt)
