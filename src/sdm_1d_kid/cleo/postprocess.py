"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: postprocess.py
Project: cleo
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Post-process ensemble of (one or more) runs of 1-D kid test case with Cleo
into single dataset.
"""

# %% ---------------------------- GLOBAL IMPORTS ----------------------------- #
import logging

logger = logging.getLogger(__name__)


# %% ------------------------- FUNCTION DEFINITIONS -------------------------- #
def parse_arguments():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument("config_filename", type=Path, help="Path to .yaml config file")
    parser.add_argument("setup_filename", type=Path, help="Path to .txt setup file")
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
    """Open CLEO datasets of each ensemble member as one dataset with an "ensemble" dimension.

    Superdroplet (ragged array) data are dropped because their length differs between members.
    All members must have identical output times (join="exact" raises an error otherwise).

    Args:
        run_datasets (list[str]): Paths to the "*.zarr" dataset of each ensemble member.

    Returns:
        xarray.Dataset with dimensions (ensemble, time, gbxindex).
    """
    from pathlib import Path

    import xarray as xr

    def drop_superdroplets(ds):
        return ds.drop_dims("superdroplets", errors="ignore").drop_vars(
            "raggedcount", errors="ignore"
        )

    ds = xr.open_mfdataset(
        run_datasets,
        engine="zarr",
        combine="nested",
        concat_dim="ensemble",
        preprocess=drop_superdroplets,
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
def add_height(ds, gbxs):
    """Replace gridbox index with height of gridbox centres and add gridbox bounds and volume.

    Args:
        ds (xarray.Dataset): Dataset with "gbxindex" dimension.
        gbxs (dict): CLEO gridboxes (from cleopy's pygbxsdat.get_gridboxes).
    """
    zhalf = gbxs["zhalf"]
    volume = gbxs["gbxvols"][0, 0, :]
    assert len(volume) == len(gbxs["zfull"]), "gridbox volumes do not match heights"

    ds = ds.rename_dims({"gbxindex": "height"}).drop_vars("gbxindex")
    return ds.assign_coords(
        height=(
            "height",
            gbxs["zfull"],
            {"units": "m", "long_name": "height of gridbox centre"},
        ),
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
    ).assign(volume=("height", volume, {"units": "m^3", "long_name": "gridbox volume"}))


def add_liquid_water(ds):
    """Add droplet number concentration, liquid water content and liquid water path.

    LWP is the finite-volume sum of LWC * dz over gridboxes above the surface (z >= 0), so the
    superdroplet reservoir gridbox below the surface is excluded.

    Args:
        ds (xarray.Dataset): Dataset with CLEO mass moments ("massmom0" [], "massmom1" [g])
          and gridbox volume and bounds (see add_height).
    """
    numconc = (ds.massmom0 / ds.volume / 1e6).assign_attrs(
        units="cm^-3", long_name="droplet number concentration"
    )
    lwc = (ds.massmom1 / ds.volume).assign_attrs(
        units="g m^-3", long_name="liquid water content"
    )
    dz = ds.height_upper - ds.height_lower
    lwp = ((lwc * dz).where(ds.height_lower >= 0.0).sum("height") / 1000).assign_attrs(
        units="kg m^-2", long_name="liquid water path (above surface)"
    )

    return ds.assign(numconc=numconc, lwc=lwc, lwp=lwp)


def add_surface_precipitation(ds, obststep, rolling_window=None):
    """Add surface precipitation rate, accumulation and (optionally) rolling mean rate.

    CLEO's "precip" is the depth of water [m] moving downwards through the lower boundary of
    each gridbox per output interval, so surface precipitation is taken from the gridbox whose
    lower boundary is z = 0.

    Args:
        ds (xarray.Dataset): Dataset with "precip" and gridbox bounds (see add_height).
        obststep (float): Time between outputs [s].
        rolling_window (int, optional): Width of centred rolling mean window [number of output
          steps]. No rolling mean if None. Note: values within half a window of the
          start/end of the run are NaN.
    """
    import numpy as np

    if "precip" not in ds:
        logger.warning("no 'precip' in dataset, skipping surface precipitation")
        return ds

    is_surface = np.isclose(ds.height_lower.values, 0.0)
    assert is_surface.sum() == 1, (
        "expected exactly one gridbox with lower boundary at z = 0"
    )
    depth = (
        ds.precip.isel(height=int(np.argmax(is_surface))) * 1000
    )  # [mm per output interval]

    rate = (depth / (obststep / 3600)).assign_attrs(
        units="mm hr^-1", long_name="surface precipitation rate"
    )
    cumulative = depth.cumsum("time").assign_attrs(
        units="mm", long_name="accumulated surface precipitation"
    )
    ds = ds.assign(surfprecip_rate=rate, surfprecip_cumulative=cumulative)

    if rolling_window is not None:
        rolling = rate.rolling(time=rolling_window, center=True).mean()
        ds = ds.assign(
            surfprecip_rolling=rolling.assign_attrs(
                units="mm hr^-1",
                long_name="centred rolling mean of surface precipitation rate",
                rolling_window_nsteps=rolling_window,
            )
        )

    return ds


# ---- MetPy calculations on numpy blocks (applied lazily with xarray.apply_ufunc) ---- #
def _vapour_pressure_hpa(press_hpa, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.vapor_pressure(press_hpa * units.hPa, qvap_kgkg * units("kg/kg"))
        .to("hPa")
        .magnitude
    )


def _theta_k(press_hpa, temp_k):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.potential_temperature(press_hpa * units.hPa, temp_k * units.kelvin)
        .to("K")
        .magnitude
    )


def _theta_virtual_k(press_hpa, temp_k, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.virtual_potential_temperature(
            press_hpa * units.hPa, temp_k * units.kelvin, qvap_kgkg * units("kg/kg")
        )
        .to("K")
        .magnitude
    )


def _density_kgm3(press_hpa, temp_k, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.density(
            press_hpa * units.hPa, temp_k * units.kelvin, qvap_kgkg * units("kg/kg")
        )
        .to("kg/m^3")
        .magnitude
    )


def _relative_humidity_percent(press_hpa, temp_k, qvap_kgkg):
    from metpy import calc as mtpy_calc
    from metpy.units import units

    return (
        mtpy_calc.relative_humidity_from_mixing_ratio(
            press_hpa * units.hPa, temp_k * units.kelvin, qvap_kgkg * units("kg/kg")
        )
        .to("percent")
        .magnitude
    )


def _lazy(func, *arrays):
    """Apply numpy function to (dask-backed) DataArrays block-wise without loading them."""
    import xarray as xr

    return xr.apply_ufunc(func, *arrays, dask="parallelized", output_dtypes=[float])


def add_thermodynamics(ds):
    """Add vapour and dry-air pressures, densities, potential temperatures and relative humidity.

    Densities and virtual potential temperature include water vapour but not condensate.

    Args:
        ds (xarray.Dataset): Dataset with CLEO thermodynamics ("press" [hPa], "temp" [K]
          and "qvap" [g/kg]).
    """
    qvap = ds.qvap / 1000  # [kg/kg]
    press_vapour = _lazy(_vapour_pressure_hpa, ds.press, qvap)
    rho = _lazy(_density_kgm3, ds.press, ds.temp, qvap)

    return ds.assign(
        relh=_lazy(_relative_humidity_percent, ds.press, ds.temp, qvap).assign_attrs(
            units="%", long_name="relative humidity"
        ),
        theta=_lazy(_theta_k, ds.press, ds.temp).assign_attrs(
            units="K", long_name="potential temperature"
        ),
        theta_virtual=_lazy(_theta_virtual_k, ds.press, ds.temp, qvap).assign_attrs(
            units="K", long_name="virtual potential temperature"
        ),
        rho=rho.assign_attrs(
            units="kg m^-3", long_name="density of moist air (dry air + vapour)"
        ),
        rho_dry=(rho / (1 + qvap)).assign_attrs(
            units="kg m^-3", long_name="density of dry air"
        ),
        press_vapour=press_vapour.assign_attrs(
            units="hPa", long_name="partial pressure of water vapour"
        ),
        press_dry=(ds.press - press_vapour).assign_attrs(
            units="hPa", long_name="partial pressure of dry air"
        ),
    )


def order_dimensions(ds):
    """Order dimensions of every variable as (ensemble, height, time)."""
    return ds.transpose("ensemble", "height", "time", ..., missing_dims="ignore")


def add_setup_metadata(ds, setup_filename):
    """Copy the contents of CLEO's setup .txt file into the dataset's global attributes."""
    from pathlib import Path

    return ds.assign_attrs(
        cleo_setup=Path(setup_filename).read_text(),
        cleo_setup_filename=str(setup_filename),
    )


# %% ----------------------- CREATE ENSEMBLE DATASET ------------------------- #
def create_ensemble_dataset(config, gbxs, run_datasets, precip_rolling_window=None):
    """Combine CLEO output of an ensemble of (one or more) runs into one post-processed dataset.

    Args:
        config (dict): CLEO configuration (from setup .txt file).
        gbxs (dict): CLEO gridboxes (from grid file).
        run_datasets (list): Paths to the ".zarr" dataset of each ensemble member.
        precip_rolling_window (int, optional): Width of rolling mean of surface
          precipitation rate [number of output steps]; no rolling mean if None.

    Returns:
        xarray.Dataset with dimensions (ensemble, height, time).
    """
    if not run_datasets:
        raise ValueError("Cannot postprocess Cleo runs: no run datasets given")
    run_datasets = sorted(str(d) for d in run_datasets)

    return (
        open_ensemble(run_datasets)
        .pipe(add_height, gbxs)
        .pipe(add_liquid_water)
        .pipe(add_surface_precipitation, config["OBSTSTEP"], precip_rolling_window)
        .pipe(add_thermodynamics)
        .pipe(order_dimensions)
    )


# %% -------------------------------- MAIN ----------------------------------- #
def main(
    config_filename,
    setup_filename,
    run_datasets,
    postprocpath=None,
):
    import yaml
    from cleopy.sdmout_src import pygbxsdat, pysetuptxt

    with open(config_filename) as config_file:
        cnfg = yaml.safe_load(config_file)
    grid_filename = cnfg["inputfiles"]["grid_filename"]
    ensemble_name = cnfg["run_ensemble"]["ensemble_name"]

    config = pysetuptxt.get_config(setup_filename, nattrs=3, isprint=False)
    consts = pysetuptxt.get_consts(setup_filename, isprint=False)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=False)
    precip_rolling_window = 100  # [number of OBSTSTEP timesteps, 1 timestep~1.25s]

    ds = create_ensemble_dataset(
        config, gbxs, run_datasets, precip_rolling_window=precip_rolling_window
    ).pipe(add_setup_metadata, setup_filename)

    if postprocpath is not None:
        if not postprocpath.is_dir():
            raise ValueError(
                f"Cannot write post-processed dataset for Cleo: parent directory doesn't exist, {postprocpath}"
            )
        filename = postprocpath / f"{ensemble_name}.zarr"
        ds.to_zarr(filename, mode="w", zarr_format=2)
        logger.info(f"post-processed ensemble dataset written to: {filename}")

    return ds


# %% --------------------------- RUN PROGRAM --------------------------------- #
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.config_filename,
        args.setup_filename,
        args.run_datasets,
        postprocpath=args.postprocpath,
    )
