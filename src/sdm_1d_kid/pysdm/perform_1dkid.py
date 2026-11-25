"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: perform_1dkid.py
Project: pysdm
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
interface called to run the 1-D KiD with specifically PySDM as the microphysics
scheme, and then plot a quickplot of the moisture field results.
"""


def main(
    kid_config,
    sds_init,
    sdm_settings,
    outputdata,
    figpath,
    run_name,
):
    from pathlib import Path

    import numpy as np

    import sdm_1d_kid.pysdm.helpers.quickplot_moisture as qpm
    from sdm_1d_kid.pysdm.helpers.pysdm_kid_dynamics import PySDMKiDDynamics

    print("--- Running PySDM 1-D KiD Rainshaft Model ---")
    kid = PySDMKiDDynamics(kid_config, sds_init, sdm_settings, outputdata)
    output = kid.run().products
    for variable, values in output.items():
        assert isinstance(values, np.ndarray)
        outputdir = Path(outputdata["outputdir"]).resolve()
        file = outputdir / f"{variable}.npy".replace(" ", "_")
        np.save(file, values)
    print("-----------------------------------------------")

    print("--- Quick-plotting Results ---")
    assert figpath.exists(), "The specified figpath does not exist."
    assert run_name, "The run_name cannot be empty."
    qpm.quickplot_moisture(
        outputdata,
        kid_config,
        sdm_settings,
        figpath,
        run_name,
    )
    print("------------------------------")

    return output
