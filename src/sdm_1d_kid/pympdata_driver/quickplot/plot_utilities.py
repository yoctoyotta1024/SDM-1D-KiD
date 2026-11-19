"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: plot_utilities.py
Project: quickplot
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Helpful functions for plotting
"""


def save_figure(fig, figpath, figname):
    """
    Save a Matplotlib figure as a PNG file with high resolution and tight bounding box.

    Args:
        fig (matplotlib.figure.Figure): The Matplotlib figure to be saved.
        figpath (Path): The directory where the figure will be saved.
        figname (str): The name of the PNG file to save in figpath directory.

    Returns:
        None

    """
    from pathlib import Path

    filename = Path(figpath) / figname
    fig.savefig(
        filename,
        dpi=400,
        bbox_inches="tight",
        facecolor="w",
        format="png",
    )
    print("Figure .png saved as: " + str(figpath) + "/" + figname)
