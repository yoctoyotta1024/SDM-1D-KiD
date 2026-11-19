"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: generate_inputfiles.py
Project: cleo
Created Date: Wednesday 18th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script to generate a real configuration file (by copying and editing a draft one), the initial
superdroplets data ("initsupers_filename"), and the grid data ("grid_filename") which Cleo
requires as model input data, using cleopy modules.
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logger = logging.getLogger(__name__)


# %% ------------------------- FUNCTION DEFINITIONS -------------------------- #
def parse_arguments():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config_filename", type=Path, help="Path to destination (real) YAML config file"
    )
    parser.add_argument(
        "src_config_filename",
        type=Path,
        default=None,
        help="Path to source (draft) YAML config file",
    )
    parser.add_argument(
        "config_params",
        type=dict,
        default=None,
        help="Dictionary of configuration parameters to set",
    )
    parser.add_argument(
        "figlabel",
        type=str,
        default="",
        help="Label to save initialiation figures with (if save_figures is True)",
    )
    parser.add_argument(
        "--gen_config",
        action="store_true",  # default is False
        help="Generate YAML configuration file",
    )
    parser.add_argument(
        "--gen_gbxs",
        action="store_true",  # default is False
        help="Generate gridbox boundaries binary file conditions",
    )
    parser.add_argument(
        "--gen_supers",
        action="store_true",  # default is False
        help="Generate initial superdroplet conditions binary file",
    )
    parser.add_argument(
        "--print_info",
        action="store_true",  # default is False
        help="Print details of grid and superdroplet initialisation",
    )
    parser.add_argument(
        "--show_figures",
        action="store_true",  # default is False
        help="Show initialiation figures",
    )
    parser.add_argument(
        "--save_figures",
        action="store_true",  # default is False
        help="Save initialiation figures in savefigpath",
    )
    parser.add_argument(
        "--savefigpath",
        type=Path,
        default=None,
        help="Directory to save initialiation figures in (if save_figures is True)",
    )
    return parser.parse_args()


def generate_configuration_file(src_config_filename, config_filename, config_params):
    import shutil

    from cleopy import editconfigfile

    ### --- copy src_config_filename into config_filename --- ##
    config_filename.unlink(missing_ok=True)  # delete any existing config
    shutil.copy(src_config_filename, config_filename)

    ### --- edit selected parameters of config_filename --- ###
    editconfigfile.edit_config_params(config_filename, config_params)


def generate_gridfile(
    config_filename, print_info, show_figures, save_figures, savefigpath, figlabel
):
    from pathlib import Path

    import numpy as np
    import yaml
    from cleopy import geninitconds

    ### --- load configuration --- ###
    with open(config_filename) as config_file:
        cnfg = yaml.safe_load(config_file)
    constants_filename = cnfg["inputfiles"]["constants_filename"]
    grid_filename = cnfg["inputfiles"]["grid_filename"]
    grid_geom = cnfg["generate_inputfiles"]["grid_geometry"]

    if not Path(constants_filename).exists():
        raise ValueError(
            "Cannot generate grid file: Cleo constants .hpp file doesn't exist"
        )
    if not Path(grid_filename).parent.is_dir():
        raise ValueError(
            "Cannot generate grid file: parent directory for grid file doesn't exist"
        )
    if save_figures and not Path(savefigpath).is_dir():
        raise ValueError(
            "Cannot generate grid file: directory for figures doesn't exist"
        )

    ### --- input parameters for z, x, and y coords of gridbox boundaries --- ###
    zgrid = [
        float(grid_geom["zmin"]),
        float(grid_geom["zmax"]),
        float(grid_geom["zdelta"]),
    ]
    xgrid = np.asarray([0, 1])
    ygrid = np.asarray([0, 1])

    ### --- binary file generation --- ###
    geninitconds.generate_gridbox_boundaries(
        grid_filename,
        zgrid,
        xgrid,
        ygrid,
        constants_filename,
        isprintinfo=print_info,
        isfigures=[show_figures, save_figures],
        savefigpath=savefigpath,
        savelabel=figlabel,
    )


def generate_superdroplet_initial_conditions_file(
    config_filename,
    print_info,
    show_figures,
    save_figures,
    savefigpath,
    figlabel,
):
    from pathlib import Path

    import yaml
    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import attrsgen, crdgens, dryrgens
    from PySDM.initialisation import spectra

    from sdm_1d_kid.cleo.helpers import alpha_sampling_wrapper as asw
    # from sdm_1d_kid.cleo.helpers import get_initial_pressure_profile ###WIP

    ### --- load configuration --- ###
    with open(config_filename) as config_file:
        cnfg = yaml.safe_load(config_file)
    constants_filename = cnfg["inputfiles"]["constants_filename"]
    grid_filename = cnfg["inputfiles"]["grid_filename"]
    initsupers_filename = cnfg["initsupers"]["initsupers_filename"]
    sds_init = cnfg["generate_inputfiles"]["superdroplet_initialisation"]

    if not Path(constants_filename).exists():
        raise ValueError(
            "Cannot generate superdroplet initial conditions file: Cleo constants .hpp file doesn't exist"
        )
    if not Path(grid_filename).exists():
        raise ValueError(
            "Cannot generate superdroplet initial conditions file: grid file doesn't exist"
        )
    if save_figures and not Path(savefigpath).is_dir():
        raise ValueError(
            "Cannot generate superdroplet initial conditions file: directory for figures doesn't exist"
        )

    ### --- Choice of Droplet Radius Probability Distribution and Radii Generator --- ###
    geomean = float(sds_init["geomean"])
    geosig = float(sds_init["geosig"])
    numconc = float(sds_init["numconc"])
    spectrum = spectra.Lognormal(norm_factor=1.0, m_mode=geomean, s_geom=geosig)
    nsupers_pergbx = int(sds_init["nsupers_pergbx"])
    alpha = float(sds_init["alpha"])  ###WIP, issue with randomness?
    default_cdf_range = (0.00001, 0.99999)
    rspan = spectrum.percentiles(
        default_cdf_range
    )  # min and max range of radii to sample [m]
    radiigen, xiprobdist = asw.AlphaSamplingWrapper(spectrum, alpha, rspan)
    numconc_tolerance = (
        0.001  # 0.1% tolerance on actual numconc not being equal to input numconc
    )
    xi_by_pressure = False  ###WIP

    ### --- Initialise Superdroplets According to Pressure profile --- ###
    # if xi_by_pressure:
    #     is_exner_novapour = False  # (!) Settings here MUST match kid_dynamics.py (!)
    #     is_exner_novapour_uniformrho = False  # (!) Settings MUST match kid_dynamics.py (!)
    #     p_surf = 1000  # [hPa] PSURF (!) Settings here MUST match kid_dynamics.py (!)
    #     z_min = -25  # [m] (!) Settings here MUST match those in run_cleo_1dkid.py (!)
    #     z_max = 3200  # [m] (!) Settings here MUST match those in run_cleo_1dkid.py (!)
    #     z_delta = 25  # [m] (!) Settings here MUST match those in run_cleo_1dkid.py (!)
    #     press, press_ref = get_initial_pressure_profile(
    #         grid_filename=grid_filename,
    #         is_exner_novapour=is_exner_novapour,
    #         is_exner_novapour_uniformrho=is_exner_novapour_uniformrho,
    #         p_surf=p_surf,
    #         z_min=z_min,
    #         z_max=z_max,
    #         z_delta=z_delta,
    #     )
    # else:
    #     press, press_ref = None, 0.0
    press, press_ref = None, 0.0

    ### --- Choice of Droplet Dry Radii (i.e. mass of solute) Generator --- ###
    dryr_sf = 1.0  # scale factor for dry radii [m]
    dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)  # dryradii are 1/sf of radii [m]

    ### --- Choice of Droplet Spatial Coordinate Generators --- ###
    coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 range randomly or not
    coord1gen = None  # do not generate superdroplet coord1s
    coord2gen = None  # do not generate superdroplet coord2s

    ### --- Extra Settings for Plots of Initial Superdroplets --- ###
    if show_figures or save_figures:
        gbxs2plt = [
            0,
            64,
            128,
        ]  # indexes of gridboxes of SDs to plot (nb. "all" can be very slow)

    ### --- binary file generation --- ###
    initattrsgen = attrsgen.AttrsGenerator(
        radiigen,
        dryradiigen,
        xiprobdist,
        coord3gen,
        coord1gen,
        coord2gen,
        xi_by_pressure=xi_by_pressure,
        press=press,
        press_ref=press_ref,
    )
    geninitconds.generate_initial_superdroplet_conditions(
        initattrsgen,
        initsupers_filename,
        config_filename,
        constants_filename,
        grid_filename,
        nsupers_pergbx,
        numconc,
        numconc_tolerance=numconc_tolerance,
        isprintinfo=print_info,
        isfigures=[show_figures, save_figures],
        savefigpath=savefigpath,
        gbxs2plt=gbxs2plt,
        savelabel=figlabel,
    )


# %% -------------------------------- MAIN ----------------------------------- #
def main(
    config_filename,
    gen_config=False,
    gen_gbxs=False,
    gen_supers=False,
    src_config_filename=None,
    config_params=None,
    print_info=False,
    show_figures=False,
    save_figures=False,
    savefigpath=None,
    figlabel="",
):
    if gen_config:
        if src_config_filename is None or config_filename is None:
            raise ValueError(
                "Cannot generate config file: source and/or destination config filename(s) missing"
            )
        elif not src_config_filename.exists():
            raise ValueError(
                "Cannot generate config file: source config file doesn't exist"
            )
        else:
            logger.info("generating configuration file")
            generate_configuration_file(
                src_config_filename, config_filename, config_params
            )

    if gen_gbxs:
        if not config_filename.exists():
            raise ValueError("Cannot generate grid file: config doesn't exist")
        else:
            logger.info("generation grid file")
            generate_gridfile(
                config_filename,
                print_info,
                show_figures,
                save_figures,
                savefigpath,
                figlabel,
            )

    if gen_supers:
        if not config_filename.exists():
            raise ValueError(
                "Cannot generate superdroplet initial conditions file: config doesn't exist"
            )
        else:
            logger.info("generating superdroplet initial conditions file ")
            generate_superdroplet_initial_conditions_file(
                config_filename,
                print_info,
                show_figures,
                save_figures,
                savefigpath,
                figlabel,
            )


# %% --------------------------- RUN PROGRAM --------------------------------- #
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.config_filename.resolve(),
        gen_config=args.gen_config,
        gen_gbxs=args.gen_gbxs,
        gen_supers=args.gen_supers,
        src_config_filename=args.src_config_filename.resolve(),
        config_params=args.config_params,
        print_info=args.print_info,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
        savefigpath=args.savefigpath,
        figlabel=args.figlabel,
    )
