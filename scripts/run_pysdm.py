"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: run_pysdm.py
Project: scripts
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
This script runs the entire workflow for the PySDM 1-D KiD (generate input,
run model, post-process output and plot quickplots).
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logging.basicConfig(
    format="%(levelname)s: %(name)s(%(funcName)s): %(message)s", level=logging.INFO
)
logger = logging.getLogger("run_pysdm." + __name__)


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
        help="Run PySDM 1-D KiD",
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
def make_directories(path2ensemble, tmppath, binpath, postprocpath):
    if not path2ensemble.exists():
        logger.warning(
            "path to ensemble doesn't exist, making directory and its parents"
        )

    path2ensemble.mkdir(exist_ok=True, parents=True)
    tmppath.mkdir(exist_ok=True, parents=False)
    binpath.mkdir(exist_ok=True, parents=False)
    postprocpath.mkdir(exist_ok=True, parents=False)


# %% ----------------------- MAIN FUNCTION DEFINITIONS ----------------------- #
def inputfiles(tmppath, binpath, ensemble_name, ensemble_size):
    from pathlib import Path

    import numpy as np
    from run_cleo import generate_random_alpha_seeds

    from sdm_1d_kid.pysdm import generate_inputfiles

    ### grid_geometry
    zmin = 0  # minimum z-coord of domain [m]
    zmax = 3200  # maximum z-coord of domain [m]
    zdelta = 25  # z-spacing of gridboxes in domain [m]

    ### superdroplet_initialisation
    ### *(!)Note:* if changing geomean or geosig, be mindful of MINRADIUS/MAXRADIUS in config.yaml
    numconc = 50e6  # number conc. of initial droplet spectrum [m^-3]
    geomean = 4.0e-8  # geometric mean radius of initial droplet spectrum [m]
    geosig = 1.4  # geometric standard deviation of initial droplet spectrum
    nsupers_pergbx = 256  # no. of super-droplets per gridbox to initialise
    alpha = 1.0  # sampling param: 0 -> const xi, 1 -> xi follows spectrum
    alpha_seeds = generate_random_alpha_seeds(
        ensemble_size
    )  # seeds for random number generator for alpha sampling (seed = an integer)

    ### SDM Runtime Parameter
    kappa = 0.9  # Aerosol hygroscopicity factor
    is_precip = True  # enable/disable precipitation scheme
    n_iters = 1  # number of iterations to run
    # collision_kernel = Hydrodynamic(),  ###WIP

    ### timesteps
    timestep = 1.25  # time between SD collision [s]
    t_end = 3600  # time span of integration from 0s to T_END [s]

    ### KiD Dynamics
    wmax_const = 3  # circa. max vertical velocity (is maximum w*rho_dry [kg m^-2 s^-1])
    tscale_const = 600  # timescale of sinusoid [s]
    psurf_hpa = 1000  # surface pressure [hPa]

    ### KiD Dynamics Initial Atmospheric Profile (default is hydrostatic equilibrium)
    is_approx_drhod_dz = (
        False  # ignore water vapour component in hydrostatic equilibrium profile
    )

    ### Output Parameters ###
    output_timestep = 1.25  # time between SDM observations [s]
    cloud_water_radius_range = (
        2e-6,
        50e-6,
    )  # range of droplet sizes to identify as "cloud" water [m]
    rain_water_radius_range = (
        50e-6,
        np.inf,
    )  # range of droplet sizes to identify as "rain" water [m]

    for n in range(ensemble_size):
        ### --- INPUT PARAMETERS FOR EACH RUN --- ###
        ### labels and filenames for single run
        run_name = f"{ensemble_name}_{n}"
        config_filename = tmppath / f"config_{run_name}.yaml"
        outputdir = binpath / f"{run_name}"

        ### --- GENERATE INPUT FILES --- ###
        path2sdm1dkid = Path(__file__).parent.parent.resolve()
        src_config_filename = (
            path2sdm1dkid / "src" / "extern" / "pysdm" / "draft_config.yaml"
        )

        config_params = {
            "ensemble_name": ensemble_name,
            "ensemble_size": ensemble_size,
            "run_name": run_name,
            "geomean": geomean,
            "geosig": geosig,
            "numconc": numconc,
            "nsupers_pergbx": nsupers_pergbx,
            "alpha": alpha,
            "alpha_seed": alpha_seeds[n],
            "zmin": zmin,
            "zmax": zmax,
            "zdelta": zdelta,
            "wmax_const": wmax_const,
            "tscale_const": tscale_const,
            "psurf_hpa": psurf_hpa,
            "is_approx_drhod_dz": is_approx_drhod_dz,
            "timestep": timestep,
            "t_end": t_end,
            "is_precip": is_precip,
            "kappa": kappa,
            "n_iters": n_iters,
            "outputdir": str(outputdir),
            "output_timestep": output_timestep,
            "cloud_water_radius_range": cloud_water_radius_range,
            "rain_water_radius_range": rain_water_radius_range,
        }

        generate_inputfiles.main(
            config_filename=config_filename,
            src_config_filename=src_config_filename,
            config_params=config_params,
        )


def run_model(path2build):
    print("###WIP")


def postprocess(binpath, postprocpath):
    print("###WIP")


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
    path2ensemble = path2build / "pysdm" / ensemble_name
    tmppath = path2ensemble / "tmp"
    binpath = path2ensemble / "bin"
    postprocpath = path2ensemble / "postproc"

    if do_inputfiles:
        logger.info("generating inputfiles")
        make_directories(path2ensemble, tmppath, binpath, postprocpath)
        inputfiles(tmppath, binpath, ensemble_name, ensemble_size)

    if do_run_model:
        logger.info("running PySDM 1-D KiD model")
        run_model(path2build)

    if do_postprocess:
        logger.info("post-processing PySDM 1-D KiD output")
        postprocess(binpath, postprocpath)

    if do_plot_results:
        logger.info("plotting results of PySDM 1-D KiD")
        # plot_results() ###WIP

    if (
        not do_inputfiles
        and not do_run_model
        and not do_postprocess
        and not do_plot_results
    ):
        logger.warning("run PySDM called without any action")


# %% ----------------------------- RUN PROGRAM ------------------------------- #
if __name__ == "__main__":
    main()
