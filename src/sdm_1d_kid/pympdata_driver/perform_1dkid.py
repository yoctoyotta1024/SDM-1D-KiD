"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: perform_1dkid.py
Project: pympdata_driver
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
interface called to run the 1-D KiD with a generic microphysics scheme
and then plot a quickplot of the moisture field results.
"""


def main(
    z_min,
    z_max,
    z_delta,
    time_end,
    timestep,
    thermo_init,
    microphys_scheme,
    advect_hydrometeors,
    figpath,
    run_name,
):
    """
    Run test case for a 1-D KiD rainshaft model.

    This function runs a 1-D KiD rainshaft model with a specified microphysics scheme and
    KiD dynamics given the initial thermodynamics. The data is then saved/plotted in the
    figpath directory using the run_name as a label.

    Args:
        z_min (float): Lower limit of 1-D column (max half-cell) (m).
        z_max (float): Upper limit of 1-D column (max half-cell) (m).
        z_delta (float): Grid spacing of 1-D column (m).
        time_end (float): End time for the simulation (s).
        timestep (float): Timestep for the simulation (s).
        thermo_init (Thermodynamics): Initial thermodynamics.
        microphys_scheme: Microphysics scheme to use in test run.
        figpath (str): Path to the directory where data/plots will be saved.
        run_name (str): Name of the test run (used for labeling output).

    Raises:
        AssertionError: If the specified figpath does not exist or if run_name is empty.

    Returns:
        None

    """
    from pathlib import Path

    import sdm_1d_kid.pympdata_driver.quickplot as qp

    print("--- Running PYMPDTA 1-D KiD Rainshaft Model ---")
    out = run_pympdata_as_1dkid(
        z_min,
        z_max,
        z_delta,
        time_end,
        timestep,
        thermo_init,
        microphys_scheme,
        advect_hydrometeors,
    )
    print("-----------------------------------------------")

    print("--- Quick-plotting Results ---")
    assert Path(figpath).exists(), "The specified figpath does not exist."
    assert run_name, "The run_name cannot be empty."
    qp.quickplot_moisture(out, z_min, z_max, z_delta, figpath, run_name)
    print("------------------------------")

    return out


def run_pympdata_as_1dkid(
    z_min,
    z_max,
    z_delta,
    time_end,
    timestep,
    thermo,
    microphys_scheme,
    advect_hydrometeors,
):
    """Run 1-D KiD rainshaft model with a specified microphysics scheme and KiD dynamics.

    This function runs a 1-D KiD rainshaft model with the given initial
    thermodynamic conditions, and microphysics scheme from time to time_end with
    a constant timestep based on the Shipway and Hill (2012) setup.

    Parameters:
        z_min (float):
          Lower limit of 1-D column (m).
        z_max (float):
          Upper limit of 1-D column (m).
        z_delta (float):
          Grid spacing od 1-D column (m).
        time_end (float):
          End time for the simulation (s).
        timestep (float):
          Timestep for the simulation (s).
        thermo (Thermodynamics):
          Initial thermodynamic conditions.
        microphys_scheme:
          Microphysics scheme to use.

    Returns:
          OutputThermodynamics: Output containing thermodynamic data from the model run.
    """
    from sdm_1d_kid.pympdata_driver.kid_dynamics import KiDDynamics
    from sdm_1d_kid.pympdata_driver.thermo.output_thermodynamics import (
        OutputThermodynamics,
    )

    ### type of dynamics rainshaft will undergo
    kid_dynamics = KiDDynamics(
        z_min,
        z_max,
        z_delta,
        timestep,
        time_end,
        advect_hydrometeors=advect_hydrometeors,
    )

    ### run dynamics + microphysics from time to time_end
    microphys_scheme.initialize()

    ### data to output during model run
    ntime = int(time_end / timestep) + 1
    nz = len(kid_dynamics.zhalf) - 1
    shape = (ntime, nz)
    out = OutputThermodynamics(shape, zhalf=kid_dynamics.zhalf)

    time = 0.0
    thermo = kid_dynamics.set_thermo(time, thermo)
    out.output_thermodynamics(time, thermo)
    while time < time_end:
        thermo = kid_dynamics.run(time, timestep, thermo)
        thermo = microphys_scheme.run(timestep, thermo)
        kid_dynamics.set_advectees(thermo)

        time += timestep

        out.output_thermodynamics(time, thermo)

    microphys_scheme.finalize()

    out.finalize()
    return out
