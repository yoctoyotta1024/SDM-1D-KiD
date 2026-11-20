"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: get_initial_pressure_profile.py
Project: helpers
Created Date: Friday 20th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Method to get dictionary of {gbxindex: pressure} to set initial pressure profile
when superdroplet xi (number concentration) initialisation requires it
(see generate_inputfiles.py)
"""


def get_initial_pressure_profile(grid_filename, kid_config):
    import numpy as np
    from cleopy.gbxboundariesbinary_src import read_gbxboundaries as rgrid
    from PyMPDATA_examples.Shipway_and_Hill_2012 import si

    from sdm_1d_kid.pympdata_driver.settings import Settings

    ### add SI units to config variables
    p_surf = kid_config["psurf_hpa"] * si.hPa
    z_min = kid_config["zmin"] * si.m
    z_max = kid_config["zmax"] * si.m
    dz = kid_config["zdelta"] * si.m
    is_exner_novapour = kid_config["is_exner_novapour"]
    is_exner_novapour_uniformrho = kid_config["is_exner_novapour_uniformrho"]
    is_approx_drhod_dz = kid_config["is_approx_drhod_dz"]

    ### dummy settings just to get presure profile
    settings = Settings(
        dt=None,
        dz=dz,
        wmax_const=None,
        tscale_const=None,
        t_max=None,
        z_min=z_min,
        z_max=z_max,
        p_surf=p_surf,
        is_exner_novapour=is_exner_novapour,
        is_exner_novapour_uniformrho=is_exner_novapour_uniformrho,
        is_approx_drhod_dz=is_approx_drhod_dz,
    )
    zfull = np.arange(z_min + dz / 2, z_max + dz / 2, dz)
    press_prof = settings.press0(zfull)

    gbxbounds = rgrid.read_dimless_gbxboundaries_binary(grid_filename, isprint=False)

    assert len(press_prof) == len(gbxbounds.keys()), (
        "number of gbxindexes and number of pressure values are not equal"
    )

    return dict(zip(gbxbounds.keys(), press_prof)), press_prof[0]
