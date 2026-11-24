"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: run.py
Project: pysdm
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Run 1-D kid test case with PySDM with either only condensation/evaporation
enabled or with precipitation enabled (with_precip =
condensation/evaporation + droplet collisions + droplet fall speeds).
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

    import yaml

    from sdm_1d_kid.pysdm import perform_1dkid

    ### --- load configuration --- ###
    with open(config_filename) as config_file:
        cnfg = yaml.safe_load(config_file)
    outputdir = Path(cnfg["outputdata"]["outputdir"]).resolve()

    if not outputdir.is_dir():
        raise ValueError("Cannot run PySDM: directory for outputdata doesn't exist")

    ### --- read and check KiD configuration --- ###
    sds_init = cnfg["generate_inputfiles"]["superdroplet_initialisation"]
    kid_config = cnfg["kid_config"]
    sdm_settings = cnfg["sdm_settings"]
    outputdata = cnfg["outputdata"]
    run_name = cnfg["run_ensemble"]["run_name"]
    quickplotdir = outputdir.parent

    if not quickplotdir.is_dir():
        raise ValueError("Cannot run PySDM: directory for quickplots doesn't exist")

    ### --- run 1-D KiD with PyMPDATA and PySDM --- ###
    perform_1dkid.main(
        kid_config,
        sds_init,
        sdm_settings,
        outputdata,
        quickplotdir,
        run_name,
    )


# %% --------------------------- RUN PROGRAM --------------------------------- #
if __name__ == "__main__":
    args = parse_arguments()
    main(args.config_filename)
