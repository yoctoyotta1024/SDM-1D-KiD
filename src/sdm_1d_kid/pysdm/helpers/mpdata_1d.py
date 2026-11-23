"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: mpdata_1d.py
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
https://github.com/open-atmos/PySDM/blob/main/examples/PySDM_examples/featured/Shipway_and_Hill_2012/mpdata_1d.py
between v3.1.0 and v3.1.1
"""

from functools import cached_property

import numpy as np
from PyMPDATA import ScalarField, Solver, Stepper, VectorField
from PyMPDATA.boundary_conditions import Extrapolated
from PySDM.impl import arakawa_c


class MPDATA_1D:
    def __init__(
        self,
        nz,
        dt,
        options,
        advector_of_t,
        advectee_of_zZ_at_t0,
        g_factor_of_zZ,
    ):
        self.__t = 0
        self.dt = dt

        self._options = options

        self.advector_of_t = advector_of_t

        self._grid = (nz,)
        bcs = (Extrapolated(),)
        zZ_scalar = arakawa_c.z_scalar_coord(self._grid) / nz
        self._g_factor = ScalarField(
            data=g_factor_of_zZ(zZ_scalar),
            halo=self._options.n_halo,
            boundary_conditions=bcs,
        )
        self._advector = VectorField(
            data=(np.full(nz + 1, advector_of_t(0)),),
            halo=self._options.n_halo,
            boundary_conditions=bcs,
        )
        self._advectee = ScalarField(
            data=advectee_of_zZ_at_t0(zZ_scalar),
            halo=self._options.n_halo,
            boundary_conditions=bcs,
        )

    @cached_property
    def solver(self):
        return Solver(
            stepper=Stepper(
                options=self._options, grid=self._grid, non_unit_g_factor=True
            ),
            advectee=self._advectee,
            advector=self._advector,
            g_factor=self._g_factor,
        )

    @property
    def solver_cached(self):
        return "solver" in self.__dict__

    @property
    def advectee(self):
        return (self.solver.advectee if self.solver_cached else self._advectee).get()

    @property
    def advector(self):
        return (
            self.solver.advector if self.solver_cached else self._advector
        ).get_component(0)

    def update_advector_field(self):
        self.__t += 0.5 * self.dt
        self.advector[:] = self.advector_of_t(self.__t)
        np.testing.assert_array_less(np.abs(self.advector), 1)
        self.__t += 0.5 * self.dt

    def __call__(self, _):
        self.solver.advance(1)
