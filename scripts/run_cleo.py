"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: run_cleo.py
Project: scripts
Created Date: Thursday 1st January 1970
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
This script runs the entire workflow for the Cleo 1-D KiD (generate inputfiles,
run model, post-process output and plot quickplots. It assumes you have already
compiled Cleo's python bindings.
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logging.basicConfig(
    format="%(levelname)s: %(name)s(%(funcName)s): %(message)s", level=logging.INFO
)
logger = logging.getLogger("run_cleo." + __name__)


# %% ---------------------------- PARSE ARGUMENTS ---------------------------- #
def parse_arguments():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument("path2build", type=Path, help="Path to build directory")
    parser.add_argument("ensemble_name", type=str, help="name of ensemble")
    parser.add_argument("ensemble_size", type=int, help="number of runs in ensemble")
    parser.add_argument(
        "--do_inputfiles",
        action="store_true",  # default is False
        help="Generate config, grid, and superdroplet initial condition files",
    )
    parser.add_argument(
        "--do_run_model",
        action="store_true",  # default is False
        help="Run Cleo 1-D KiD",
    )
    parser.add_argument(
        "--do_post_process",
        action="store_true",  # default is False
        help="Post-process output",
    )
    parser.add_argument(
        "--do_plot_results",
        action="store_true",  # default is False
        help="Plot quickplot results",
    )

    return parser.parse_args()


# %% ---------------------- HELPER FUNCTION DEFINITIONS ---------------------- #
def make_directories(path2ensemble, tmppath, sharepath, binpath, savefigpath=None):
    if not path2ensemble.exists():
        logger.warning(
            "path to ensemble doesn't exist, making directory and its parents"
        )

    path2ensemble.mkdir(exist_ok=True)
    tmppath.mkdir(exist_ok=True, parents=False)
    sharepath.mkdir(exist_ok=True, parents=False)
    binpath.mkdir(exist_ok=True, parents=False)

    if savefigpath is not None:
        savefigpath.mkdir(exist_ok=True, parents=False)


# %% ----------------------- MAIN FUNCTION DEFINITIONS ----------------------- #
def inputfiles(path2build, path2ensemble, ensemble_name, ensemble_size):
    from pathlib import Path

    import numpy as np

    from sdm_1d_kid.cleo import generate_inputfiles

    ### --- INPUT PARAMETERS FOR ENSEMBLE --- ###
    gen_config = True
    gen_gbxs = True
    gen_supers = True

    ### directories and filenames for ensemble
    tmppath = path2ensemble / "tmp"
    sharepath = path2ensemble / "share"
    binpath = path2ensemble / "bin"
    savefigpath = binpath
    constants_filename = (
        path2build / "_deps" / "cleo-src" / "libs" / "cleoconstants.hpp"
    )
    grid_filename = sharepath / "dimlessGBxboundaries.dat"

    ### grid_geometry
    zmin = -25  # minimum z-coord of domain [m]
    zmax = 3200  # maximum z-coord of domain [m]
    zdelta = 25  # z-spacing of gridboxes in domain [m]

    ### superdroplet_initialisation
    numconc = 50e6  # number conc. of initial droplet spectrum [m^-3]
    nsupers_pergbx = 256  # no. of super-droplets per gridbox to initialise
    alpha = 0.0  # sampling param: 0 -> const xi, 1 -> xi follows spectrum

    ### timesteps
    CONDTSTEP = 0.1  # time between SD condensation [s]
    COLLTSTEP = 1.25  # time between SD collision [s]
    MOTIONTSTEP = 1.25  # time between SDM motion [s]
    COUPLTSTEP = 1.25  # time between dynamic couplings [s]
    OBSTSTEP = 1.25  # time between SDM observations [s]
    T_END = 3600  # time span of integration from 0s to T_END [s]

    ### python_bindings
    enable_terminal_velocity = (
        True  # true enables terminal velocity in superdroplet motion
    )
    enable_condensation = True  # true enables condensation in microphysics
    enable_collisions = True  # true enables collisions in microphysics

    ### domain
    ngbxs = int(np.ceil((zmax - zmin) / zdelta))  # total number of gridboxes
    initnsupers = ngbxs * nsupers_pergbx  # initial no. of super-droplets to initialise
    maxnsupers = initnsupers * 2  # maximum number of SDs

    ### boundary_conditions
    newnsupers = nsupers_pergbx  # number SDs to add to each gridbox below/above LOWER/UPPER_COORD3LIM
    NUMCONC_a = numconc  # number conc. of 1st droplet lognormal dist [m^-3]

    ### booleans for figures of grid and superdroplet initial conditions
    print_info = False
    show_figures = False
    save_figures = True

    for n in range(ensemble_size):
        ### --- INPUT PARAMETERS FOR EACH RUN --- ###
        ### labels and filenames for single run
        run_name = f"{ensemble_name}_{n}"
        config_filename = tmppath / f"config_{run_name}.yaml"
        initsupers_filename = sharepath / f"dimlessSDsinit_{run_name}.dat"
        setup_filename = binpath / f"setup_{run_name}.txt"
        zarrbasedir = binpath / f"sol_{run_name}.zarr"
        figlabel = f"_{run_name}"

        ### --- GENERATE INPUT FILES --- ###
        path2sdm1dkid = Path(__file__).parent.parent.resolve()
        src_config_filename = (
            path2sdm1dkid / "src" / "extern" / "cleo" / "draft_config.yaml"
        )
        config_params = {
            "ensemble_name": ensemble_name,
            "ensemble_size": ensemble_size,
            "run_name": run_name,
            "constants_filename": str(constants_filename),
            "grid_filename": str(grid_filename),
            "initsupers_filename": str(initsupers_filename),
            "setup_filename": str(setup_filename),
            "zarrbasedir": str(zarrbasedir),
            "zmin": zmin,
            "zmax": zmax,
            "zdelta": zdelta,
            "numconc": numconc,
            "nsupers_pergbx": nsupers_pergbx,
            "alpha": alpha,
            "ngbxs": ngbxs,
            "initnsupers": initnsupers,
            "maxnsupers": maxnsupers,
            "CONDTSTEP": CONDTSTEP,
            "COLLTSTEP": COLLTSTEP,
            "MOTIONTSTEP": MOTIONTSTEP,
            "COUPLTSTEP": COUPLTSTEP,
            "OBSTSTEP": OBSTSTEP,
            "T_END": T_END,
            "newnsupers": newnsupers,
            "NUMCONC_a": NUMCONC_a,
            "enable_terminal_velocity": enable_terminal_velocity,
            "enable_condensation": enable_condensation,
            "enable_collisions": enable_collisions,
        }
        make_directories(
            path2ensemble, tmppath, sharepath, binpath, savefigpath=savefigpath
        )
        generate_inputfiles.main(
            config_filename,
            gen_config=gen_config,
            gen_gbxs=gen_gbxs,
            gen_supers=gen_supers,
            src_config_filename=src_config_filename,
            config_params=config_params,
            print_info=print_info,
            show_figures=show_figures,
            save_figures=save_figures,
            savefigpath=savefigpath,
            figlabel=figlabel,
        )

        gen_gbxs = False  # only need one grid file for whole ensemble of runs


# %% ----------------------------- MAIN PROGRAM ------------------------------ #
def main():
    ### --- SCRIPT INPUT PARAMETERS --- ###
    args = parse_arguments()
    path2build = args.path2build.resolve()
    do_inputfiles = args.do_inputfiles
    do_run_model = args.do_run_model
    do_post_process = args.do_post_process
    do_plot_results = args.do_plot_results
    ensemble_name = args.ensemble_name
    ensemble_size = args.ensemble_size

    path2ensemble = path2build / ensemble_name

    if do_inputfiles:
        logger.info("generating inputfiles")
        inputfiles(path2build, path2ensemble, ensemble_name, ensemble_size)

    if do_run_model:
        logger.info("running Cleo 1-D KiD model")
        # run_model() ###WIP

    if do_post_process:
        logger.info("post-processing Cleo 1-D KiD output")
        # post_process() ###WIP

    if do_plot_results:
        logger.info("plotting results of Cleo 1-D KiD")
        # plot_results() ###WIP

    if (
        not do_inputfiles
        and not do_run_model
        and not do_post_process
        and not do_plot_results
    ):
        logger.warning("run Cleo called without any action")


# %% ----------------------------- RUN PROGRAM ------------------------------- #
if __name__ == "__main__":
    main()
