"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: pysdm_kid_dynamics.py
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
https://github.com/open-atmos/PySDM/blob/main/examples/PySDM_examples/featured/Shipway_and_Hill_2012/simulation.py
between v3.1.0 and v3.1.1
"""

from collections import namedtuple
from typing import Any

import numpy as np
import PySDM.products as PySDM_products
from PyMPDATA_examples.Shipway_and_Hill_2012 import si
from PySDM import Particulator
from PySDM.dynamics import (
    AmbientThermodynamics,
    Coalescence,
    Condensation,
    Displacement,
    EulerianAdvection,
)
from PySDM.impl.mesh import Mesh
from PySDM.initialisation.sampling import spatial_sampling

import sdm_1d_kid.common as cmn
from sdm_1d_kid.pysdm.helpers.kinematic_1d import Kinematic1D
from sdm_1d_kid.pysdm.helpers.mpdata_1d import MPDATA_1D
from sdm_1d_kid.pysdm.helpers.settings import Settings


class PySDMKiDDynamics:
    def __init__(self, kid_config, sds_init, sdm_settings, outputdata):
        ### read in config directly to settings
        self.do_save_spectrum = outputdata["do_save_spectrum"]
        self.do_save_attrs = outputdata["do_save_attrs"]
        alpha = float(sds_init["alpha"])
        alpha_seed = int(sds_init["alpha_seed"])
        xi_by_dry_air_density = sds_init["xi_by_dry_air_density"]
        dt_out = float(outputdata["output_timestep"])
        save_spec_and_attr_times = np.arange(
            0, float(kid_config["t_end"]) + dt_out, dt_out
        )
        cloud_water_radius_range = tuple(
            float(r) * si.m for r in outputdata["cloud_water_radius_range"]
        )
        rain_water_radius_range = tuple(
            float(r) * si.m for r in outputdata["rain_water_radius_range"]
        )
        settings = Settings(
            dt=float(kid_config["timestep"]) * si.s,
            dz=float(kid_config["zdelta"]) * si.m,
            wmax_const=float(kid_config["wmax_const"]) * si.m / si.s * si.kg / si.m**3,
            tscale_const=float(kid_config["tscale_const"]) * si.s,
            t_max=float(kid_config["t_end"]) * si.s,
            z_max=float(kid_config["zmax"]) * si.m,
            p_surf=float(kid_config["psurf_hpa"]) * si.hPa,
            is_exner_novapour=kid_config["is_exner_novapour"],
            is_exner_novapour_uniformrho=kid_config["is_exner_novapour_uniformrho"],
            is_approx_drhod_dz=kid_config["is_approx_drhod_dz"],
            n_sd_per_gridbox=int(sds_init["nsupers_pergbx"]),
            particles_per_volume_STP=float(sds_init["numconc"]) / si.m**3,
            m_mode=float(sds_init["geomean"]) * si.m,
            s_geom=float(sds_init["geosig"]),
            kappa=float(sdm_settings["kappa"]),
            precip=sdm_settings["is_precip"],
            collision_kernel=sdm_settings["collision_kernel"],
            save_spec_and_attr_times=save_spec_and_attr_times,
            cloud_water_radius_range=cloud_water_radius_range,
            rain_water_radius_range=rain_water_radius_range,
        )

        # *NOTE:*(!)must match above SI units^
        self.written_kid_config_units = {
            "zmax": "m",
            "zdelta": "m",
            "timestep": "s",
            "t_end": "s",
            "wmax_const": "kg m^-2 s^-1 (peak rho_d*w)",
            "tscale_const": "s",
            "psurf_hpa": "hPa",
        }

        ### Hard-coded settings
        self.nt = settings.nt
        self.z_min = settings.z_min
        self.save_spec_and_attr_times = settings.save_spec_and_attr_times
        self.number_of_bins = settings.number_of_bins
        self.r_bins_edges = settings.r_bins_edges
        self.r_bins_edges_dry = settings.r_bins_edges_dry
        self.particulator = None
        self.output_attributes = None
        self.output_products = None
        self.mesh = Mesh(
            grid=(settings.nz,),
            size=(settings.z_max + settings.particle_reservoir_depth,),
        )

        ### initialise options and MPDATA
        options = cmn.get_options(is_accurate=True)
        mpdata = MPDATA_1D(
            nz=settings.nz,
            dt=settings.dt,
            options=options,
            advector_of_t=lambda t: settings.rhod_w(t) * settings.dt / settings.dz,
            advectee_of_zZ_at_t0=lambda zZ: settings.qv0(settings.zposition(zZ)),
            g_factor_of_zZ=lambda zZ: settings.rhod(settings.zposition(zZ)),
        )

        ### initialise PySDM
        env = Kinematic1D(
            dt=settings.dt,
            mesh=self.mesh,
            thd_of_z=settings.thd,
            rhod_of_z=settings.rhod,
            alpha_seed=alpha_seed,
            z0=settings.z_min,
            solvers=mpdata,
        )

        _extra_nz = settings.particle_reservoir_depth // settings.dz
        _z_vec = settings.dz * np.linspace(
            -_extra_nz, settings.nz - _extra_nz, settings.nz + 1
        )
        self.g_factor_vec = settings.rhod(_z_vec)

        dynamics: list[Any] = [AmbientThermodynamics()]

        if settings.enable_condensation:
            dynamics.append(
                Condensation(
                    adaptive=settings.condensation_adaptive,
                    rtol_thd=settings.condensation_rtol_thd,
                    rtol_x=settings.condensation_rtol_x,
                    update_thd=settings.condensation_update_thd,
                )
            )
        dynamics.append(EulerianAdvection())

        self.products = []
        if settings.precip:
            self.add_collision_dynamic(dynamics, settings, self.products)

        dynamics.append(
            Displacement(
                enable_sedimentation=settings.precip,
                precipitation_counting_level_index=int(
                    settings.particle_reservoir_depth / settings.dz
                ),
            )
        )

        alpha_sampling = cmn.get_alpha_sampling(
            settings.spectrum, alpha, size_range=settings.size_range
        )
        self.attributes = env.init_attributes(
            spatial_discretisation=spatial_sampling.Pseudorandom(),
            spectral_sampling_method=cmn.get_spectral_sampling_method(alpha_sampling),
            kappa=settings.kappa,
            collisions_only=not settings.enable_condensation,
            z_part=settings.z_part,
            n_sd=settings.n_sd,
        )

        if not xi_by_dry_air_density:
            rhod = env["rhod"].to_ndarray()
            renormalisation_factor = (
                settings.formulae.constants.rho_STP / rhod[self.attributes["cell id"]]
            )
            self.attributes["multiplicity"] *= renormalisation_factor

        self.products += [
            PySDM_products.WaterMixingRatio(
                name="liquid water mixing ratio",
                unit="g/kg",
            ),
            PySDM_products.WaterMixingRatio(
                name="cloud water mixing ratio",
                unit="g/kg",
                radius_range=settings.cloud_water_radius_range,
            ),
            PySDM_products.WaterMixingRatio(
                name="rain water mixing ratio",
                unit="g/kg",
                radius_range=settings.rain_water_radius_range,
            ),
            PySDM_products.AmbientDryAirDensity(name="rhod"),
            PySDM_products.AmbientDryAirPotentialTemperature(name="thd"),
            PySDM_products.ParticleSizeSpectrumPerVolume(
                name="wet spectrum", radius_bins_edges=settings.r_bins_edges
            ),
            PySDM_products.ParticleConcentration(
                name="nc", radius_range=settings.cloud_water_radius_range
            ),
            PySDM_products.ParticleConcentration(
                name="nr", radius_range=settings.rain_water_radius_range
            ),
            PySDM_products.ParticleConcentration(
                name="na", radius_range=(0, settings.cloud_water_radius_range[0])
            ),
            PySDM_products.SuperDropletCountPerGridbox(),
            PySDM_products.AmbientRelativeHumidity(name="RH", unit="%"),
            PySDM_products.AmbientPressure(name="p"),
            PySDM_products.AmbientTemperature(name="T"),
            PySDM_products.AmbientWaterVapourMixingRatio(
                name="water_vapour_mixing_ratio"
            ),
        ]
        if settings.enable_condensation:
            self.products.extend(
                [
                    PySDM_products.PeakSaturation(),
                    PySDM_products.ParticleSizeSpectrumPerVolume(
                        name="dry spectrum",
                        radius_bins_edges=settings.r_bins_edges_dry,
                        dry=True,
                    ),
                    PySDM_products.ParticleVolumeVersusRadiusLogarithmSpectrum(
                        name="dry volume spectrum",
                        radius_bins_edges=settings.r_bins_edges_dry,
                        dry=True,
                    ),
                ]
            )
        if settings.precip:
            self.products.extend(
                [
                    PySDM_products.SurfacePrecipitation(),
                ]
            )
        self.particulator = Particulator(
            n_sd=settings.n_sd,
            environment=env,
            dynamics=dynamics,
            attributes=self.attributes,
            products=tuple(self.products),
        )

        self.output_attributes = {
            "cell origin": [],
            "position in cell": [],
            "radius": [],
            "multiplicity": [],
        }
        self.output_products = {}
        for k, v in self.particulator.products.items():
            if len(v.shape) == 0:
                self.output_products[k] = np.zeros(self.nt + 1)
            elif len(v.shape) == 1:
                self.output_products[k] = np.zeros((self.mesh.grid[-1], self.nt + 1))
            elif len(v.shape) == 2:
                number_of_time_sections = len(self.save_spec_and_attr_times)
                self.output_products[k] = np.zeros(
                    (self.mesh.grid[-1], self.number_of_bins, number_of_time_sections)
                )
        self.output_products["t_spectra"] = np.full(
            len(self.save_spec_and_attr_times), np.nan
        )

        ### print configuration
        key = f"dz={settings.dz}, dt={settings.dt}, options={options}"
        print(f"Simulating {settings.nt} timesteps using {key}")
        print(self._format_pysdm_kid_config(kid_config, nz=settings.nz, nt=settings.nt))

    def _format_pysdm_kid_config(self, kid_config, nz=None, nt=None):
        """Return KiD configuration as aligned 'name = value [units]' lines for printing."""
        width = max(len(key) for key in kid_config)
        lines = ["KiDDynamics configuration:"]
        for key, value in kid_config.items():
            units = self.written_kid_config_units.get(key)
            lines.append(
                f"  {key:<{width}} = {value}" + (f" [{units}]" if units else "")
            )
        if nz is not None and nt is not None:
            lines.append(f"  {'(derived)':<{width}} = {nz} gridboxes, {nt} timesteps")
        return "\n".join(lines)

    @staticmethod
    def add_collision_dynamic(dynamics, settings, _):
        dynamics.append(
            Coalescence(
                collision_kernel=settings.collision_kernel,
                adaptive=settings.coalescence_adaptive,
            )
        )

    def save_scalar(self, step):
        for k, v in self.particulator.products.items():
            if len(v.shape) > 1:
                continue
            if len(v.shape) == 1:
                self.output_products[k][:, step] = v.get()
            else:
                self.output_products[k][step] = v.get()

    def save_spectrum(self, index):
        for k, v in self.particulator.products.items():
            if len(v.shape) == 2:
                self.output_products[k][:, :, index] = v.get()

    def save_attributes(self):
        for k, v in self.output_attributes.items():
            v.append(self.particulator.attributes[k].to_ndarray())

    def save(self, step):
        self.save_scalar(step)
        time = step * self.particulator.dt
        if len(self.save_spec_and_attr_times) > 0 and (
            np.min(
                np.abs(
                    np.ones_like(self.save_spec_and_attr_times) * time
                    - np.array(self.save_spec_and_attr_times)
                )
            )
            < 0.1
        ):
            save_index = np.argmin(
                np.abs(
                    np.ones_like(self.save_spec_and_attr_times) * time
                    - np.array(self.save_spec_and_attr_times)
                )
            )
            if self.do_save_spectrum:
                self.save_spectrum(save_index)  # droplet size/mass distribs
                self.output_products["t_spectra"][save_index] = time
            if self.do_save_attrs:
                self.save_attributes()  # superdroplet attributes

    def run(self):
        mesh = self.particulator.mesh

        assert "t" not in self.output_products and "z" not in self.output_products
        self.output_products["t"] = np.linspace(
            0, self.nt * self.particulator.dt, self.nt + 1, endpoint=True
        )
        self.output_products["z"] = np.linspace(
            self.z_min + mesh.dz / 2,
            self.z_min + (mesh.grid[-1] - 1 / 2) * mesh.dz,
            mesh.grid[-1],
            endpoint=True,
        )
        spectra_edges = {
            "wet_spectrum_radius_bins_edges": self.r_bins_edges,
            "dry_spectrum_radius_bins_edges": self.r_bins_edges_dry,
        }
        for name, edges in spectra_edges.items():
            assert name not in self.output_products
            self.output_products[name] = np.asarray(edges)

        self.save(0)
        print_interval = max(1, self.nt // 10)
        for step in range(self.nt):
            if (step + 1) % print_interval == 0:
                print(f"PySDM on step {step}: {step / self.nt * 100:.2f}% complete")
            mpdata = self.particulator.environment.solvers
            mpdata.update_advector_field()
            if "Displacement" in self.particulator.dynamics:
                self.particulator.dynamics["Displacement"].upload_courant_field(
                    (mpdata.advector / self.g_factor_vec,)
                )
            self.particulator.advance(n_steps=1)
            self.save(step + 1)

        Outputs = namedtuple("Outputs", "products attributes")
        output_results = Outputs(self.output_products, self.output_attributes)
        return output_results
