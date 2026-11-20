"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: postprocess.py
Project: bulk_scheme
Created Date: Friday 20th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Post-process ensemble of (one or more) runs of 1-D kid test case with the
PyMPDATA bulk scheme into single dataset. The steps, variable names and units follow
``sdm_1d_kid.cleo.postprocess`` so that both datasets can be compared with the same
diagnostics. The thermodynamic diagnostics, dimension ordering and units corrections
are imported from there so that they are calculated identically for both models.
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logger = logging.getLogger(__name__)


# %% ------------------------- FUNCTION DEFINITIONS -------------------------- #
def parse_arguments():
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config",
        type=json.loads,
        help="Dictionary of configuration for bulk scheme KiD run",
    )
    parser.add_argument(
        "run_datasets", nargs="+", type=Path, help="Paths to .zarr datasets"
    )
    parser.add_argument(
        "--postprocpath",
        type=Path,
        default=None,
        help="Path to directory to write post-processed ensemble Zarr dataset in",
    )
    return parser.parse_args()


# %% ------------------------ OPEN ENSEMBLE OF RUNS -------------------------- #
def open_ensemble(run_datasets):
    """Open bulk scheme datasets of each ensemble member as one dataset with an "ensemble"
    dimension.

    All members must have identical output times (join="exact" raises an error otherwise).

    Args:
        run_datasets (list[str]): Paths to the "*.zarr" dataset of each ensemble member.

    Returns:
        xarray.Dataset with dimensions (ensemble, zfull, time) and coordinate "zhalf".
    """
    from pathlib import Path

    import xarray as xr

    ds = xr.open_mfdataset(
        run_datasets,
        engine="zarr",
        combine="nested",
        concat_dim="ensemble",
        consolidated=False,
        data_vars="all",
        coords="minimal",
        compat="override",
        join="exact",
    )

    members = [Path(d).stem for d in run_datasets]
    return ds.assign_coords(
        ensemble=("ensemble", members, {"long_name": "ensemble member (run name)"})
    ).assign(
        sources=(
            "ensemble",
            run_datasets,
            {"long_name": "path to dataset of each ensemble member"},
        )
    )


# %% ------------------------- POST-PROCESSING STEPS ------------------------- #
def add_height(ds):
    """Replace "zfull" with height of gridbox centres and add gridbox bounds (from "zhalf").

    The coordinates are named as in the post-processed Cleo dataset. Unlike Cleo, the bulk
    scheme has no gridbox volume (it is not needed for any of its diagnostics).

    Args:
        ds (xarray.Dataset): Dataset with "zfull" dimension and "zhalf" coordinate [m].
    """
    zhalf = ds.zhalf.values
    assert len(zhalf) == ds.sizes["zfull"] + 1, "gridbox bounds do not match heights"

    ds = ds.rename({"zfull": "height"}).drop_vars("zhalf")
    return ds.assign_coords(
        height=ds.height.assign_attrs(units="m", long_name="height of gridbox centre"),
        height_lower=(
            "height",
            zhalf[:-1],
            {"units": "m", "long_name": "height of gridbox lower boundary"},
        ),
        height_upper=(
            "height",
            zhalf[1:],
            {"units": "m", "long_name": "height of gridbox upper boundary"},
        ),
    )


def drop_ice(ds):
    """Drop the mass mixing ratios of ice-phase hydrometeors (cloud ice, snow and graupel).

    The 1-D KiD test case here is warm (no ice microphysics), so these are always zero and are
    not in Cleo's output.

    Args:
        ds (xarray.Dataset): Dataset with (some of) "qice", "qsnow" and "qgrau".
    """
    return ds.drop_vars(["qice", "qsnow", "qgrau"], errors="ignore")


def convert_units(ds):
    """Convert pressure and mass mixing ratios to the units of Cleo's output.

    Pressure [Pa] -> [hPa] and mass mixing ratios [kg/kg] -> [g/kg]; all other attributes
    are kept. Also writes the units of the dry air density in the same style as Cleo's
    post-processed variables ("kg m-3" -> "kg m^-3").

    Args:
        ds (xarray.Dataset): Dataset with "press" [Pa], "rhod" [kg m-3] and mass mixing
          ratios ("qvap", "qcond", ...) [kg/kg].
    """
    converted = {
        "press": (ds.press / 100).assign_attrs({**ds.press.attrs, "units": "hPa"})
    }
    for name in ("qvap", "qcond", "qice", "qrain", "qsnow", "qgrau"):
        if name in ds:
            converted[name] = (ds[name] * 1000).assign_attrs(
                {**ds[name].attrs, "units": "g/kg"}
            )
    converted["rhod"] = ds.rhod.assign_attrs(
        units="kg m^-3", long_name="dry air density (prescribed by KiD driver)"
    )
    return ds.assign(converted)


def add_liquid_water(ds):
    """Add liquid water content and liquid water path.

    The bulk scheme's mass mixing ratios are per unit mass of dry air (as in PyMPDATA's
    Shipway and Hill 2012 example), so LWC = rhod * qcond. LWP is the finite-volume sum of
    LWC * dz over gridboxes above the surface (z >= 0), as for Cleo.

    Args:
        ds (xarray.Dataset): Dataset with "rhod" [kg m^-3], "qcond" [g/kg] and gridbox
          bounds (see add_height).
    """
    lwc = (ds.rhod * ds.qcond).assign_attrs(
        units="g m^-3",
        long_name="liquid water content",
    )  # [kg m^-3] * [g/kg] = [g m^-3]
    dz = ds.height_upper - ds.height_lower
    lwp = ((lwc * dz).where(ds.height_lower >= 0.0).sum("height") / 1000).assign_attrs(
        units="kg m^-2",
        long_name="liquid water path (above surface)",
    )

    return ds.assign(lwc=lwc, lwp=lwp)


def add_config_metadata(ds, config):
    """Copy the run configuration (e.g. from run_bulk_scheme.py) into the global attributes.

    Values which are not plain numbers or strings (e.g. paths) are converted to strings.

    Args:
        ds (xarray.Dataset): Dataset to add the metadata to.
        config (dict): Configuration of (any) one run of the ensemble.
    """
    import json

    return ds.assign_attrs(bulk_scheme_config=json.dumps(config, default=str, indent=2))


def check_dry_density(ds, rtol=1e-3):
    """Log whether the dry air density recomputed by MetPy ("rho_dry") matches
    the density prescribed by the KiD driver ("rhod").

    The two are calculated differently (the driver's formulae and constants vs. MetPy's
    from p, T and qvap), so they need not agree exactly; a large difference indicates an
    inconsistency between the driver's thermodynamics and the post-processing (e.g. units).
    """
    import numpy as np

    reldiff = float((np.abs(ds.rho_dry - ds.rhod) / ds.rhod).max())
    if reldiff > rtol:
        logger.warning(
            f"rho_dry (MetPy) differs from rhod (KiD driver) by up to {reldiff:.2e} "
            "(relative to rhod)"
        )
        return False
    else:
        logger.info(f"rho_dry (MetPy) matches rhod (KiD driver) ({reldiff:.2e})")
        return True


# %% ----------------------- CREATE ENSEMBLE DATASET ------------------------- #
def create_ensemble_dataset(run_datasets):
    """Combine bulk scheme output of an ensemble of (one or more) runs into one
    post-processed dataset.

    Args:
        run_datasets (list): Paths to the ".zarr" dataset of each ensemble member.

    Returns:
        xarray.Dataset with dimensions (ensemble, height, time).
    """
    from sdm_1d_kid.cleo.postprocess import (
        add_thermodynamics,
        metadata_corrections,
        order_dimensions,
    )

    if not run_datasets:
        raise ValueError("Cannot postprocess bulk scheme runs: no run datasets given")
    run_datasets = sorted(str(d) for d in run_datasets)

    return (
        open_ensemble(run_datasets)
        .pipe(add_height)
        .pipe(drop_ice)
        .pipe(convert_units)
        .pipe(add_liquid_water)
        .pipe(add_thermodynamics)
        .pipe(order_dimensions)
        .pipe(metadata_corrections)
    )


# %% -------------------------------- MAIN ----------------------------------- #
def main(
    cnfg,
    run_datasets,
    postprocpath=None,
):
    """Post-process an ensemble of bulk scheme runs and (optionally) write it to
    'postprocpath/<ensemble_name>.zarr'.

    Args:
        config (dict): Configuration of (any) one run, stored as metadata.
        run_datasets (list): Paths to the ".zarr" dataset of each ensemble member.
        postprocpath (pathlib.Path, optional): Directory to write the dataset in.
    """
    ensemble_name = cnfg["ensemble_name"]
    ds = create_ensemble_dataset(run_datasets).pipe(add_config_metadata, cnfg)

    is_rhod_equal = check_dry_density(ds, rtol=1e-4)
    if is_rhod_equal:
        ds = ds.drop_vars(["rhod"], errors="ignore")

    if postprocpath is not None:
        if not postprocpath.is_dir():
            raise ValueError(
                f"Cannot write post-processed dataset for bulk scheme: parent directory doesn't exist, {postprocpath}"
            )
        filename = postprocpath / f"{ensemble_name}.zarr"
        ds.to_zarr(filename, mode="w", zarr_format=2)
        logger.info(f"post-processed ensemble dataset written to: {filename}")

    return ds


# %% --------------------------- RUN PROGRAM --------------------------------- #
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.config,
        args.run_datasets,
        postprocpath=args.postprocpath,
    )
