"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: cleo_sdm_wrapper.py
Project: helpers
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
wrapper function for an instance of CleoSDM microphysics scheme so it can be
used by generic test cases and run scripts, e.g. with PyMPDATA.
NOTE: To use the wrapper, you must first export "CLEO_PYTHON_BINDINGS".
E.g. if python bindings are built in $HOME/SDM-1D-KiD/build/, do:
export CLEO_PYTHON_BINDINGS=$HOME/SDM-1D-KiD/build/_deps/cleo-build/cleo_python_bindings/
"""

import os
import sys

from sdm_1d_kid.cleo.helpers.cleo_sdm import CleoSDM
from sdm_1d_kid.pympdata_driver.thermo.thermodynamics import Thermodynamics

sys.path.append(os.environ["CLEO_PYTHON_BINDINGS"])
import cleo_python_bindings as cleo


class CleoSDMWrapper:
    """A class wrapping around C++ bindings to CLEO's Superdroplet Model (SDM)
    microphysics scheme (wrapper for compatibility purposes).

    This class wraps around the CLEO SDM microphysics to provide compatibility
    with the python run scripts and tests in this project. It initializes CLEO
    SDM microphysics object and provides wrappers around methods to initialize,
    finalize, and run the microphysics.
    """

    def __init__(
        self,
        config_filename,
        is_motion,
        t_start,
        timestep,
        press,
        temp,
        qvap,
        qcond,
        wvel,
        uvel,
        vvel,
        do_init=True,
    ):
        """Initialize the CleoSDMWrapper.

        This Wrapper only works correctly if addresses of press, temp, qvap, qcond,
        wvel, uvel, and vvel arrays remain unchanged throughout a simulation.
        Undefined behaviour if values are changed by reassigning arrays rather
        than by copying data into the arrays given during wrapper initialisation.
        """
        config = cleo.Config(str(config_filename))
        if do_init:
            cleo.cleo_initialize(config)

        self.microphys = CleoSDM(
            config,
            is_motion,
            t_start,
            timestep,
            press,
            temp,
            qvap,
            qcond,
            wvel,
            uvel,
            vvel,
        )
        self.name = "Wrapper around " + self.microphys.name

        # constants to de-dimensionalise thermodynamics
        self.TEMP0 = 273.15  # Temperature [K]
        self.P0 = 100000.0  # Pressure [Pa]
        self.W0 = 1.0  # Velocity [m/s]

    def initialize(self) -> int:
        """Initialise the microphysics scheme.

        This method calls the microphysics initialisation

        Returns:
            int: 0 upon successful initialisation
        """

        return 0

    def finalize(self) -> int:
        """Finalise the microphysics scheme.

        This method calls the microphysics finalisation.

        Returns:
            int: 0 upon successful finalisation.
        """

        return 0

    def run(self, timestep: float, thermo: Thermodynamics) -> Thermodynamics:
        """Run the microphysics computations.

        This method is a wrapper of the MicrophysicsScheme object's run function to call the
        microphysics computations in a way that's compatible with the test and scripts in this project.

        Args:
            timestep (float):
              Time-step for integration of microphysics (s)
            thermo (Thermodynamics):
              Thermodynamic properties.

        Returns:
            Thermodynamics: Updated thermodynamic properties after microphysics computations.

        """
        # de-dimensionlise variables
        thermo.press /= self.P0
        thermo.temp /= self.TEMP0
        thermo.wvel /= self.W0
        thermo.uvel /= self.W0
        thermo.vvel /= self.W0

        self.microphys.run(timestep)

        # re-dimensionlise variables
        thermo.press *= self.P0
        thermo.temp *= self.TEMP0
        thermo.wvel *= self.W0
        thermo.uvel *= self.W0
        thermo.vvel *= self.W0

        return thermo
