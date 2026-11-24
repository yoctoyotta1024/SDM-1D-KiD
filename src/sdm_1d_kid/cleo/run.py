"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: run.py
Project: cleo
Created Date: Wednesday 18th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Run 1-D kid test case with Cleo with either only condensation/evaporation
enabled or with precipitation enabled (with_precip =
condensation/evaporation + droplet collisions + droplet fall speeds).
NOTE: This script assumes you have already set the path to your cleo_python_bindings
module as the environment variable ``CLEO_PYTHON_BINDINGS``, e.g.
export CLEO_PYTHON_BINDINGS=/path/to/your/build/_deps/cleo-build/cleo_python_bindings/
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logger = logging.getLogger(__name__)


# %% ------------------------- FUNCTION DEFINITIONS -------------------------- #
def parse_arguments():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument("config_filename", type=Path, help="Path to YAML config file")
    return parser.parse_args()


def main(config_filename):
    from pathlib import Path

    import numpy as np
    import yaml

    from sdm_1d_kid.cleo.helpers.cleo_sdm_wrapper import CleoSDMWrapper
    from sdm_1d_kid.pympdata_driver import perform_1dkid
    from sdm_1d_kid.pympdata_driver.thermo.thermodynamics import Thermodynamics

    ### --- load configuration --- ###
    with open(config_filename) as config_file:
        cnfg = yaml.safe_load(config_file)
    grid_filename = Path(cnfg["inputfiles"]["grid_filename"]).resolve()
    constants_filename = Path(cnfg["inputfiles"]["constants_filename"]).resolve()
    initsupers_filename = Path(cnfg["initsupers"]["initsupers_filename"]).resolve()
    setup_filename = Path(cnfg["outputdata"]["setup_filename"]).resolve()
    zarrbasedir = Path(cnfg["outputdata"]["zarrbasedir"]).resolve()

    if (
        not grid_filename.exists()
        or not constants_filename.exists()
        or not initsupers_filename.exists()
    ):
        raise ValueError(
            "Cannot run Cleo: grid, initsupers and/or constants file(s) don't exist"
        )
    elif (
        not Path(setup_filename).parent.is_dir()
        or not Path(zarrbasedir).parent.is_dir()
    ):
        raise ValueError("Cannot run Cleo: directory for outputdata doesn't exist")

    ### --- read and check KiD configuration --- ###
    t_start = 0.0
    is_motion = True  # True == advect SDs according to wind fields
    advect_hydrometeors = False  # True == PyMPDATA advection (not same as Cleo(!))
    kid_config = cnfg["kid_config"]
    timestep = float(kid_config["timestep"])
    assert float(kid_config["timestep"]) == float(cnfg["timesteps"]["COUPLTSTEP"]), (
        "timesteps don't agree"
    )
    assert float(kid_config["t_end"]) == float(cnfg["timesteps"]["T_END"]), (
        "end times don't agree"
    )

    run_name = cnfg["run_ensemble"]["run_name"]
    quickplotdir = Path(
        zarrbasedir
    ).parent  # path for PyMPDATA quickplot of 1-D KiD run

    if not quickplotdir.is_dir():
        raise ValueError("Cannot run Cleo: directory for quickplots doesn't exist")

    zmin = float(kid_config["zmin"])
    zmax = float(kid_config["zmax"])
    zdelta = float(kid_config["zdelta"])
    assert (zmax - zmin) % zdelta == 0, "z limit is not a multiple of the grid spacing."
    ngbxs = int(np.ceil((zmax - zmin) / zdelta))
    assert ngbxs == int(cnfg["domain"]["ngbxs"]), "ngbxs don't agree"

    ### --- set initial classes (will be updated by PyMPDATA and Cleo) --- ###
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
        np.tile(np.zeros(ngbxs), 2),
        np.tile(np.zeros(ngbxs), 2),
        np.tile(np.zeros(ngbxs), 2),
    )

    microphys_scheme = CleoSDMWrapper(
        config_filename,
        is_motion,
        t_start,
        timestep,
        thermo_init.press,
        thermo_init.temp,
        thermo_init.massmix_ratios["qvap"],
        thermo_init.massmix_ratios["qcond"],
        thermo_init.wvel,
        thermo_init.uvel,
        thermo_init.vvel,
    )

    ### --- run 1-D KiD with PyMPDATA and Cleo --- ###
    perform_1dkid.main(
        kid_config,
        thermo_init,
        microphys_scheme,
        advect_hydrometeors,
        quickplotdir,
        run_name,
    )


# %% --------------------------- RUN PROGRAM --------------------------------- #
if __name__ == "__main__":
    args = parse_arguments()
    main(args.config_filename)
