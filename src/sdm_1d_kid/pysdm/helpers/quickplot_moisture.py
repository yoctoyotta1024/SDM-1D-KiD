"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: quickplot_moisture.py
Project: helpers
Created Date: Monday 23rd November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""


def quickplot_moisture(
    outputdata,
    kid_config,
    sdm_settings,
    figpath,
    run_name,
):
    from pathlib import Path

    import numpy as np

    from sdm_1d_kid.pympdata_driver.quickplot import plot_utilities

    ### load relevant variables
    outputdir = Path(outputdata["outputdir"]).resolve()
    dataset = {
        var: np.load(outputdir / f"{var}.npy")
        for var in ("t", "z", "cloud_water_mixing_ratio", "rain_water_mixing_ratio")
    }
    dataset["q_cond"] = (
        dataset["cloud_water_mixing_ratio"] + dataset["rain_water_mixing_ratio"]
    )

    ### choose plotting styles
    if sdm_settings["is_precip"]:
        line = {7: ":", 9: "--", 20: "-", 30: "-.", 60: ":"}
        color = {
            7: "black",
            9: "darkblue",
            20: "mediumorchid",
            30: "crimson",
            60: "darkred",
        }
    else:
        line = {5: ":", 7: "--", 9: "-", 11: "-.", 60: ":"}
        color = {5: "grey", 7: "black", 9: "darkblue", 11: "violet", 60: "darkred"}

    ### plot and save quicklook
    figname = run_name + "_moisture.png"
    fig = plot_with_iterations(
        var="q_cond",
        qlabel="q_cond",
        output=dataset,
        line=line,
        colors=color,
    )

    plot_utilities.save_figure(fig, figpath, figname)


def plot_with_iterations(
    var,
    qlabel,
    output,
    line,
    colors,
    vmin=None,
    vmax=None,
    cmin=None,
    cmax=None,
):
    import matplotlib.pyplot as plt
    import numpy as np
    from PySDM.physics import convert_to, si

    line = line or {15: ":", 20: "--", 25: "-", 30: "-."}
    colors = colors or {15: "k", 20: "tab:blue", 25: "tab:orange", 30: "tab:green"}

    dt = output["t"][1] - output["t"][0]
    dz = output["z"][1] - output["z"][0]
    tgrid = np.concatenate(((output["t"][0] - dt / 2,), output["t"] + dt / 2))
    zgrid = np.concatenate(((output["z"][0] - dz / 2,), output["z"] + dz / 2))
    convert_to(zgrid, si.km)

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(25, 5)
    ax1 = fig.add_subplot(gs[:-1, 0:4])
    if cmin is not None and cmax is not None:
        mesh = ax1.pcolormesh(
            tgrid, zgrid, output[var], cmap="BuPu", vmin=cmin, vmax=cmax
        )
    else:
        mesh = ax1.pcolormesh(tgrid, zgrid, output[var], cmap="BuPu")

    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("z [km]")
    ax1.set_ylim(0, None)

    cbar = fig.colorbar(mesh, fraction=0.05, location="top")
    cbar.set_label(qlabel)

    ax2 = fig.add_subplot(gs[:-1, 4:], sharey=ax1)
    ax2.set_xlabel(qlabel)

    last_t = 0
    for i, t in enumerate(output["t"]):
        x, z = output[var][:, i], output["z"].copy()
        convert_to(z, si.km)
        for line_t, line_s in line.items():
            if last_t < line_t * si.min <= t:
                ax2.plot(x, z, color=colors[line_t])
                if vmin is not None and vmax is not None:
                    ax1.axvline(
                        t, ymin=vmin, ymax=vmax, color=colors[line_t], ls=line_s
                    )
                else:
                    ax1.axvline(t, color=colors[line_t], ls=line_s)
        last_t = t

    return fig
