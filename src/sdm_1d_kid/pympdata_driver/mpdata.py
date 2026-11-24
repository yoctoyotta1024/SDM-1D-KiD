"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: mpdata.py
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
https://github.com/open-atmos/PyMPDATA/blob/main/examples/PyMPDATA_examples/Shipway_and_Hill_2012/mpdata.py
between v1.6.1 and v1.6.2
"""

import numpy as np
from PyMPDATA import ScalarField, Solver, Stepper, VectorField
from PyMPDATA.boundary_conditions import Constant, Extrapolated
from PyMPDATA.impl.enumerations import INNER, OUTER

from sdm_1d_kid.pympdata_driver.arakawa_c import arakawa_c


class MPDATA:
    # pylint: disable=too-few-public-methods
    def __init__(
        self,
        nz,
        dt,
        qv_of_zZ_at_t0,
        g_factor_of_zZ,
        options,
    ):
        nr = 1  # nr==1 is like a bulk scheme
        self.nr = nr
        self.t = 0
        self.dt = dt
        self.fields = ("qvap", "qcond", "qice", "qrain", "qsnow", "qgrau")

        self.options = options

        self._solvers = {}
        for k in self.fields:
            grid = (nz, nr) if nr > 1 and k == "qcond" else (nz,)

            bcs_extrapol = tuple(
                Extrapolated(dim=d)
                for d in ((OUTER, INNER) if k == "qcond" and nr > 1 else (INNER,))
            )

            bcs_zero = tuple(
                Extrapolated(dim=d)
                for d in ((OUTER, INNER) if k == "qcond" and nr > 1 else (INNER,))
            )

            stepper = Stepper(
                options=self.options, n_dims=len(grid), non_unit_g_factor=True
            )

            data = g_factor_of_zZ(arakawa_c.z_scalar_coord(grid))
            g_factor = ScalarField(
                data=data, halo=self.options.n_halo, boundary_conditions=bcs_extrapol
            )

            data = (np.zeros(nz + 1),)
            advector = VectorField(
                data=data, halo=self.options.n_halo, boundary_conditions=bcs_zero
            )
            if k == "qvap":
                data = qv_of_zZ_at_t0(arakawa_c.z_scalar_coord(grid))
                bcs = (Extrapolated(dim=INNER),)
            else:
                data = np.zeros(grid)
                bcs = (Constant(value=0),)
            advectee = ScalarField(
                data=data, halo=self.options.n_halo, boundary_conditions=bcs
            )
            self._solvers[k] = Solver(
                stepper=stepper, advectee=advectee, advector=advector, g_factor=g_factor
            )

    def __getitem__(self, k):
        return self._solvers[k]
