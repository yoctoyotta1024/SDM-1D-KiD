"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: run.py
Project: bulk_scheme
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Run 1-D kid test case with PyMPDATA bulk scheme microphysics
"""


def main(cnfg):
    """runs test of 1-D KiD rainshaft model using bulk scheme for condensation
    extracted from pyMPDATA for the microphysics scheme.

    This function sets up initial conditions and parameters for running a 1-D KiD rainshaft
    test case using the bulk microphysics scheme for condensation from the Shipway and Hill 2012
    pyMPDATA-examples example (via a wrapper). It then runs the test case as specified.
    """
    from pathlib import Path

    import numpy as np

    from sdm_1d_kid.pympdata_driver import perform_1dkid
    from sdm_1d_kid.pympdata_driver.bulk_scheme.bulk_scheme_wrapper import (
        BulkSchemeWrapper,
    )
    from sdm_1d_kid.pympdata_driver.thermo.thermodynamics import Thermodynamics

    ### --- read and check KiD configuration --- ###
    advect_hydrometeors = True  # True == PyMPDATA advection
    kid_config = cnfg["kid_config"]

    run_name = cnfg["run_name"]
    zarrbasedir = Path(cnfg["zarrbasedir"]).resolve()
    quickplotdir = Path(
        zarrbasedir
    ).parent  # path for PyMPDATA quickplot of 1-D KiD run

    if not zarrbasedir.parent.is_dir():
        raise ValueError(
            "Cannot run bulk scheme: directory for outputdata doesn't exist"
        )

    if not quickplotdir.is_dir():
        raise ValueError(
            "Cannot run bulk scheme: directory for quickplots doesn't exist"
        )

    zmin = kid_config["zmin"]
    zmax = kid_config["zmax"]
    zdelta = kid_config["zdelta"]
    assert (zmax - zmin) % zdelta == 0, "z limit is not a multiple of the grid spacing."
    ngbxs = int((zmax - zmin) / zdelta)

    ### --- set initial classes (will be updated by PyMPDATA and bulk scheme) --- ###
    null = np.array([])  # this microphysics test doesn't need winds
    thermo_init = Thermodynamics(
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        np.zeros(ngbxs),
        null,
        null,
        null,
    )

    microphys_scheme = BulkSchemeWrapper()

    ### --- run 1-D KiD with PyMPDATA and bulk scheme microphysics --- ###
    perform_1dkid.main(
        kid_config,
        thermo_init,
        microphys_scheme,
        advect_hydrometeors,
        quickplotdir,
        run_name,
        zarrbasedir=zarrbasedir,
    )
