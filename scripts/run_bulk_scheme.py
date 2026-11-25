"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: run_bulk_scheme.py
Project: scripts
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
This script runs the entire workflow for the PyMPDATA bulk scheme microphysics
1-D KiD (generate configuration, run model, post-process output and plot
quickplots)
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logging.basicConfig(
    format="%(levelname)s: %(name)s(%(funcName)s): %(message)s", level=logging.INFO
)
logger = logging.getLogger("run_bulk_scheme." + __name__)


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
        help="Create directories for model run",
    )
    parser.add_argument(
        "--do_run_model",
        action="store_true",  # default is False
        help="Run bulk scheme 1-D KiD",
    )
    parser.add_argument(
        "--do_postprocess",
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
def make_directories(path2ensemble, binpath, postprocpath):
    if not path2ensemble.exists():
        logger.warning(
            "path to ensemble doesn't exist, making directory and its parents"
        )

    path2ensemble.mkdir(exist_ok=True, parents=True)
    binpath.mkdir(exist_ok=True, parents=False)
    postprocpath.mkdir(exist_ok=True, parents=False)


# %% ----------------------- MAIN FUNCTION DEFINITIONS ----------------------- #
def generate_configuration(path2ensemble, binpath, ensemble_name, ensemble_size):
    ### grid_geometry
    zmin = -25  # minimum z-coord of domain [m]
    zmax = 3200  # maximum z-coord of domain [m]
    zdelta = 25  # z-spacing of gridboxes in domain [m]

    ### timesteps
    timestep = 1.25  # time between dynamic couplings [s]
    t_end = 3600  # time span of integration from 0s to T_END [s]

    ### KiD Dynamics
    wmax_const = 3  # circa. max vertical velocity (is maximum w*rho_dry [kg m^-2 s^-1])
    tscale_const = 600  # timescale of sinusoid [s]
    psurf_hpa = 1000  # surface pressure [hPa]

    ### KiD Dynamics Initial Atmospheric Profile (default is hydrostatic equilibrium)
    is_exner_novapour = (
        False  # use Exner profile without accounting for vapour (Shipway and Hill 2012)
    )
    is_exner_novapour_uniformrho = (
        False  # use Exner also with constant air density (Hill et al. 2023)
    )
    is_approx_drhod_dz = (
        False  # ignore water vapour component in hydrostatic equilibrium profile
    )

    configs = []
    for n in range(ensemble_size):
        run_name = f"{ensemble_name}_{n}"
        zarrbasedir = binpath / f"{run_name}.zarr"

        cnfg = {
            "ensemble_name": ensemble_name,
            "ensemble_size": ensemble_size,
            "run_name": run_name,
            "zarrbasedir": zarrbasedir,
            "kid_config": {
                "wmax_const": wmax_const,
                "tscale_const": tscale_const,
                "is_exner_novapour": is_exner_novapour,
                "is_exner_novapour_uniformrho": is_exner_novapour_uniformrho,
                "is_approx_drhod_dz": is_approx_drhod_dz,
                "psurf_hpa": psurf_hpa,
                "zmin": zmin,
                "zmax": zmax,
                "zdelta": zdelta,
                "timestep": timestep,
                "t_end": t_end,
            },
        }
        configs.append(cnfg)

    return configs


def run_model(configs):
    from sdm_1d_kid.pympdata_driver.bulk_scheme import run

    for cnfg in configs:
        run.main(cnfg)


def postprocess(configs, binpath, postprocpath):
    import glob
    import os

    from sdm_1d_kid.pympdata_driver.bulk_scheme import postprocess

    run_datasets = sorted(glob.glob(os.path.join(binpath, "*.zarr")))
    if not configs or not run_datasets:
        raise ValueError(
            f"Cannot postprocess bulk scheme: no config given and/or output data found in search of {binpath}"
        )
    else:
        postprocess.main(configs[0], run_datasets, postprocpath=postprocpath)


# %% ----------------------------- MAIN PROGRAM ------------------------------ #
def main():
    ### --- SCRIPT INPUT PARAMETERS --- ###
    args = parse_arguments()
    path2build = args.path2build.resolve()
    do_inputfiles = args.do_inputfiles
    do_run_model = args.do_run_model
    do_postprocess = args.do_postprocess
    do_plot_results = args.do_plot_results
    ensemble_name = args.ensemble_name
    ensemble_size = args.ensemble_size

    ### directories for specific ensemble
    path2ensemble = path2build / "bulk_scheme" / ensemble_name
    binpath = path2ensemble / "bin"
    postprocpath = path2ensemble / "postproc"
    configs = generate_configuration(
        path2ensemble, binpath, ensemble_name, ensemble_size
    )

    if do_inputfiles:
        logger.info("generating inputfiles")
        make_directories(path2ensemble, binpath, postprocpath)
        logger.info("generating inputfiles -> completed")

    if do_run_model:
        logger.info("running bulk scheme 1-D KiD model")
        run_model(configs)
        logger.info("running bulk scheme 1-D KiD model -> completed")

    if do_postprocess:
        logger.info("post-processing bulk scheme 1-D KiD output")
        postprocess(configs, binpath, postprocpath)
        logger.info("post-processing bulk scheme 1-D KiD output -> completed")

    if do_plot_results:
        logger.info("plotting results of bulk scheme 1-D KiD")
        # plot_results() ###WIP
        logger.info("plotting results of bulk scheme 1-D KiD -> completed")

    if (
        not do_inputfiles
        and not do_run_model
        and not do_postprocess
        and not do_plot_results
    ):
        logger.warning("run bulk scheme called without any action")


# %% ----------------------------- RUN PROGRAM ------------------------------- #
if __name__ == "__main__":
    main()
