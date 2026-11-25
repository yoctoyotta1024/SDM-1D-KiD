"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: kinematic_1d.py
Project: helpers
Created Date: Tuesday 24th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
File is adaptation of Kinematic1D class used by the Shipway and Hill 2012 example
from PySDM repo, copied from:
https://github.com/open-atmos/PySDM/blob/main/PySDM/environments/kinematic_1d.py
between v3.1.0 and v3.1.1

Single-column time-varying-updraft framework with moisture advection handled by
[PyMPDATA](http://github.com/open-atmos/PyMPDATA/)
"""

import numpy as np
from PySDM import Formulae
from PySDM.backends import CPU
from PySDM.environments.impl import register_environment
from PySDM.environments.impl.moist import Moist
from PySDM.impl import arakawa_c
from PySDM.initialisation.hygroscopic_equilibrium import equilibrate_wet_radii


def get_seeded_backend(rng):
    formulae = Formulae(seed=int(rng.integers(0, 2**32)))
    return CPU(formulae=formulae, n_dims=1)


@register_environment()
class Kinematic1D(Moist):
    def __init__(
        self, *, dt, mesh, thd_of_z, rhod_of_z, alpha_seed, z0=0, solvers=None
    ):
        rng = np.random.default_rng(alpha_seed)
        self.default_backend = get_seeded_backend(
            rng
        )  # for the environment and Particulator
        self.spatial_backend = get_seeded_backend(
            rng
        )  # for initial superdroplet positions
        self.spectral_backend = get_seeded_backend(
            rng
        )  # for initial superdroplet radii

        super().__init__(dt, mesh, [], backend=self.default_backend)
        self.thd0 = thd_of_z(z0 + mesh.dz * arakawa_c.z_scalar_coord(mesh.grid))

        rhod = rhod_of_z(z0 + mesh.dz * arakawa_c.z_scalar_coord(mesh.grid))
        self._values["current"]["rhod"] = self.default_backend.Storage.from_ndarray(
            rhod
        )
        self._tmp["rhod"] = self.default_backend.Storage.from_ndarray(rhod)
        self.solvers = solvers

    def get_water_vapour_mixing_ratio(self) -> np.ndarray:
        return self.solvers.advectee

    def get_thd(self) -> np.ndarray:
        return self.thd0

    def init_attributes(
        self,
        *,
        spatial_discretisation,
        spectral_sampling_method,
        kappa,
        z_part=None,
        collisions_only=False,
        n_sd=None,
    ):
        super().sync()
        self.notify()

        attributes = {}
        with np.errstate(all="raise"):
            (
                attributes["cell id"],
                attributes["cell origin"],
                attributes["position in cell"],
            ) = self.mesh.cellular_attributes(
                spatial_discretisation.sample(
                    backend=self.spatial_backend,
                    grid=self.mesh.grid,
                    n_sd=n_sd,
                    z_part=z_part,
                )
            )

            if collisions_only:
                v_wet, n_per_kg = spectral_sampling_method(
                    n_sd=n_sd, backend=self.spectral_backend
                )
                attributes["volume"] = v_wet
            else:
                r_dry, n_per_kg = spectral_sampling_method(
                    n_sd=n_sd, backend=self.spectral_backend
                )
                attributes["dry volume"] = self.default_backend.formulae.trivia.volume(
                    radius=r_dry
                )
                attributes["kappa times dry volume"] = attributes["dry volume"] * kappa
                r_wet = equilibrate_wet_radii(
                    r_dry=r_dry,
                    environment=self,
                    cell_id=attributes["cell id"],
                    kappa_times_dry_volume=attributes["kappa times dry volume"],
                )
                attributes["volume"] = self.default_backend.formulae.trivia.volume(
                    radius=r_wet
                )

            rhod = self["rhod"].to_ndarray()
            cell_id = attributes["cell id"]
            domain_volume = np.prod(np.array(self.mesh.size))

        attributes["multiplicity"] = n_per_kg * rhod[cell_id] * domain_volume

        return attributes

    @property
    def dv(self):
        return self.mesh.dv
