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

from collections.abc import Iterable

import numpy as np
from PySDM import Formulae
from PySDM.dynamics import condensation
from PySDM.dynamics.collisions.collision_kernels import Geometric
from PySDM.initialisation import spectra
from PySDM.physics import si
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

import sdm_1d_kid.pympdata_driver.common as cmn


class Settings:
    def __dir__(self) -> Iterable[str]:
        return (
            "n_sd_per_gridbox",
            "p_surf",
            "kappa",
            "rho_times_w_1",
            "particles_per_volume_STP",
            "dt",
            "dz",
            "precip",
            "z_max",
            "t_max",
            "cloud_water_radius_range",
            "rain_water_radius_range",
            "r_bins_edges_dry",
            "r_bins_edges",
        )

    def __init__(
        self,
        *,
        n_sd_per_gridbox: int,
        p_surf: float = 1000 * si.hPa,
        kappa: float = 1,
        rho_times_w_1: float = 2 * si.m / si.s * si.kg / si.m**3,
        is_approx_drhod_dz: bool = False,
        particles_per_volume_STP: int = 50 / si.cm**3,
        dt: float = 1 * si.s,
        dz: float = 25 * si.m,
        z_max: float = 3000 * si.m,
        z_part: tuple | None = None,
        t_max: float = 60 * si.minutes,
        precip: bool = True,
        enable_condensation: bool = True,
        formulae: Formulae = None,
        save_spec_and_attr_times=(),
        collision_kernel=None,
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

        ### Set Time-Varying Updraught
        self.rho_times_w_1 = rho_times_w_1
        t_1 = 600 * si.s
        self.rho_times_w = lambda t: (
            rho_times_w_1 * np.sin(np.pi * t / t_1) if t < t_1 else 0
        )
        apprx_w1 = rho_times_w_1 / self.formulae.constants.rho_STP

        ### Initial PySDM SDM Superdroplets Conditions
        self.particle_reservoir_depth = (
            (2 * apprx_w1 * t_1 / np.pi) // self.dz + 1
        ) * self.dz
        self.wet_radius_spectrum_per_mass_of_dry_air = spectra.Lognormal(
            norm_factor=particles_per_volume_STP / self.formulae.constants.rho_STP,
            m_mode=0.08 / 2 * si.um,
            s_geom=1.4,
        )

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

        z_span = (-self.particle_reservoir_depth, self.z_max)
        z_points = np.linspace(*z_span, 2 * self.nz + 1)

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
            t_span=z_span,
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
        self.cloud_water_radius_range = [1 * si.um, 50 * si.um]
        self.rain_water_radius_range = [50 * si.um, np.inf * si.um]
        self.save_spec_and_attr_times = save_spec_and_attr_times

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
