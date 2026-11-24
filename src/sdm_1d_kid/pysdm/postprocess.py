"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: postprocess.py
Project: pysdm
Created Date: Tuesday 24th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Post-process ensemble of (one or more) runs of 1-D kid test case with PySDM
into a single dataset.

The post-processed dataset uses the same variable names, units, dimensions and
diagnostics as the post-processed CLEO dataset (see sdm_1d_kid.cleo.postprocess)
wherever PySDM's output allows, and it also post-processes extra output variables
from PySDM which Cleo doesn't have.
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
    parser.add_argument(
        "run_datasets",
        nargs="+",
        type=Path,
        help="Paths to directories of .npy output of each run",
    )
    parser.add_argument(
        "--postprocpath",
        type=Path,
        default=None,
        help="Path to directory to write post-processed ensemble Zarr dataset in",
    )
    return parser.parse_args()


# %% ------------------------- PYSDM OUTPUT FORMAT --------------------------- #
def _radius_bin_coords(radius_edges, dim, description):
    """Coordinates of radius bins (as in CLEO's post-processing) along dimension 'dim'."""
    import numpy as np

    lower, upper = radius_edges[:-1], radius_edges[1:]
    return {
        dim: (
            dim,
            np.sqrt(lower * upper),
            {
                "units": "m",
                "long_name": "radius of bin centre (geometric mean of bounds)",
                "description": description,
            },
        ),
        f"{dim}_lower": (
            dim,
            lower,
            {
                "units": "m",
                "long_name": "lower bound of radius bin",
                "description": description,
            },
        ),
        f"{dim}_upper": (
            dim,
            upper,
            {
                "units": "m",
                "long_name": "upper bound of radius bin",
                "description": description,
            },
        ),
        f"{dim}_dlnr": (
            dim,
            np.log(upper / lower),
            {"units": "", "long_name": "width of radius bin in ln(radius)"},
        ),
    }


def _spectrum_per_bin(spectrum, radius_edges, itimes):
    """Number concentration in each radius bin [cm^-3] with dims (time_sds, height, bin).

    PySDM's ParticleSizeSpectrumPerVolume is the number concentration per unit radius,
    dN/dr [m^-3 m^-1], with shape (height, bin, time of spectra), so it is multiplied
    by the bin width dr [m] to give the number concentration in each bin.
    """
    import numpy as np

    dr = np.diff(radius_edges)  # [m]
    spectrum = np.asarray(spectrum[:, :, itimes])  # only loads sub-sampled times
    numconc = spectrum * dr[np.newaxis, :, np.newaxis] / 1e6  # [cm^-3]
    return np.transpose(numconc, (2, 0, 1)).astype(np.float32)


# %% ------------------------ OPEN ENSEMBLE OF RUNS -------------------------- #
def open_run(rundir, sddiags_every):
    """Open PySDM .npy output of a single run as a dataset with dimensions (height, time)
    and, for size distributions, (time_sds, height, radius_bin) and
    (time_sds, height, radius_bin_dry).

    Variables are converted to CLEO's post-processed names and units where there is an
    equivalent. Products that are not in the output directory (e.g. condensation or
    collision products when those processes are disabled) are skipped.

    Args:
        rundir (pathlib.Path): Directory of .npy output files of the run.
        sddiags_every (int): Size distributions at every 'sddiags_every'th spectrum time.
    """
    from pathlib import Path

    import numpy as np
    import xarray as xr

    rundir = Path(rundir)

    def load(name, mmap=False):
        file = rundir / f"{name}.npy"
        if not file.is_file():
            return None
        return np.load(file, mmap_mode="r" if mmap else None)

    time = load("t")  # [s]
    height = load("z")  # [m] gridbox centres (including reservoir below surface)
    dims = ("height", "time")

    # (raw PySDM name, post-processed name, factor to convert units, attributes)
    products = (
        ("water_vapour_mixing_ratio", "qvap", 1000, {"units": "g/kg"}),
        ("p", "press", 1 / 100, {"units": "hPa"}),
        ("T", "temp", 1, {"units": "K"}),
        ("rhod", "rhod", 1, {"units": "kg m^-3"}),
        ("thd", "thd", 1, {"units": "K"}),
        ("RH", "rh", 1, {"units": "%"}),
        ("cloud_water_mixing_ratio", "qcloud", 1, {"units": "g/kg"}),
        ("rain_water_mixing_ratio", "qrain", 1, {"units": "g/kg"}),
        ("liquid_water_mixing_ratio", "qcond", 1, {"units": "g/kg"}),
        ("na", "numconc_aerosol", 1e-6, {"units": "cm^-3"}),
        ("nc", "numconc_cloud", 1e-6, {"units": "cm^-3"}),
        ("nr", "numconc_rain", 1e-6, {"units": "cm^-3"}),
        ("super_droplet_count_per_gridbox", "nsupers", 1, {"units": ""}),
        ("peak_saturation", "peak_saturation", 1, {"units": ""}),
    )
    data_vars = {}
    for raw_name, name, factor, attrs in products:
        values = load(raw_name)
        if values is None:
            logger.info(f"no '{raw_name}' in {rundir}, skipping it")
            continue
        assert values.shape == (len(height), len(time)), (
            f"unexpected shape of {raw_name}"
        )
        data_vars[name] = (dims, values * factor, attrs)

    surfprecip = load("surface_precipitation")  # [m s^-1] mean over previous timestep
    if surfprecip is not None:
        data_vars["surfprecip_rate"] = (
            "time",
            surfprecip * 1000 * 3600,  # [mm hr^-1]
            {"units": "mm hr^-1"},
        )

    ds = xr.Dataset(
        data_vars,
        coords={
            "time": ("time", time, {"units": "s", "long_name": "time"}),
            "height": ("height", height),
        },
    )

    # size distributions at every 'sddiags_every'th time of PySDM's spectra
    spectra_times = load("t_spectra")  # [s] NaN where no spectrum was saved
    if spectra_times is None or np.isnan(spectra_times).all():
        logger.info(
            f"no spectra saved in {rundir} (do_save_spectrum: False?), skipping them; "
            "time_sds is every 'sddiags_every'th output time instead"
        )
        return ds.assign_coords(
            time_sds=(
                "time_sds",
                time[::sddiags_every],
                {"units": "s", "long_name": "time of superdroplet diagnostics"},
            )
        )

    itimes = np.arange(0, len(spectra_times), sddiags_every)
    if np.isnan(spectra_times[itimes]).any():
        raise ValueError(
            f"Cannot postprocess spectra in {rundir}: no spectrum saved at "
            f"{np.count_nonzero(np.isnan(spectra_times[itimes]))} sub-sampled times"
        )
    ds = ds.assign_coords(
        time_sds=(
            "time_sds",
            spectra_times[itimes],
            {"units": "s", "long_name": "time of superdroplet diagnostics"},
        )
    )

    spectra = (
        (
            "wet_spectrum",
            "numconc_distrib",
            "radius_bin",
            "radius of sphere with the particle's water volume (excludes solute volume)",
        ),
        (
            "dry_spectrum",
            "numconc_distrib_dry",
            "radius_bin_dry",
            "dry (solute) radius of the particle",
        ),
    )
    for raw_name, name, dim, description in spectra:
        spectrum = load(raw_name, mmap=True)
        if spectrum is None:
            logger.info(f"no '{raw_name}' in {rundir}, skipping it")
            continue
        edges = load(f"{raw_name}_radius_bins_edges")  # [m]
        if edges is None:
            raise ValueError(
                f"Cannot postprocess '{raw_name}': no radius bin edges in {rundir}"
            )
        assert spectrum.shape == (len(height), len(edges) - 1, len(spectra_times)), (
            f"shape of {raw_name} doesn't match radius bins and times of spectra"
        )
        ds = ds.assign_coords(_radius_bin_coords(edges, dim, description)).assign(
            {
                name: (
                    ("time_sds", "height", dim),
                    _spectrum_per_bin(spectrum, edges, itimes),
                    {"units": "cm^-3"},
                )
            }
        )

    return ds


def open_ensemble(run_datasets, sddiags_every):
    """Open PySDM output of each ensemble member as one dataset with an "ensemble" dimension.

    All members must have identical times and heights (join="exact" raises an error
    otherwise).

    Args:
        run_datasets (list[str]): Paths to the .npy output directory of each ensemble member.
        sddiags_every (int): Size distributions at every 'sddiags_every'th spectrum time.
    """
    from pathlib import Path

    import xarray as xr

    ds = xr.concat(
        [open_run(d, sddiags_every) for d in run_datasets],
        dim="ensemble",
        data_vars="all",
        coords="minimal",
        compat="override",
        join="exact",
    )

    members = [Path(d).name for d in run_datasets]
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
def _at_time_sds(da, time_sds):
    """Select 'da' at the times of the size distributions (allowing for float round-off
    between PySDM's output times and the times of its spectra)."""
    tolerance = 1e-6 * float(da.time.max())  # [s]
    return da.sel(time=time_sds, method="nearest", tolerance=tolerance).drop_vars(
        "time"
    )


def add_height(ds, dz):
    """Add attributes to height of gridbox centres and add gridbox bounds.

    Args:
        ds (xarray.Dataset): Dataset with "height" coordinate of gridbox centres [m].
        dz (float): Gridbox depth [m].
    """
    import numpy as np

    assert np.allclose(np.diff(ds.height.values), dz), "expected uniform gridbox depth"
    return ds.assign_coords(
        height=ds.height.assign_attrs(units="m", long_name="height of gridbox centre"),
        height_lower=(
            "height",
            ds.height.values - dz / 2,
            {"units": "m", "long_name": "height of gridbox lower boundary"},
        ),
        height_upper=(
            "height",
            ds.height.values + dz / 2,
            {"units": "m", "long_name": "height of gridbox upper boundary"},
        ),
    )


def add_liquid_water(ds):
    """Add droplet number concentration, liquid water mixing ratio, liquid water content and
    liquid water path, and the total number of superdroplets.

    PySDM's liquid water is the water mass of all particles (cloud, rain and haze), excluding
    the solute mass, as in CLEO's "qcond". (A difference is in CLEO's "lwc" and "lwp" are
    calculated from massmom1, so they include the solute mass.) PySDM's mixing ratio is per mass
    of dry air, so LWC = qcond * rhod. LWP is the finite-volume sum of LWC * dz over gridboxes
    above the surface (z >= 0), so the superdroplet reservoir below the surface is excluded.

    Args:
        ds (xarray.Dataset): Dataset with "numconc_aerosol", "numconc_cloud", "numconc_rain"
          [cm^-3], "qcond" [g/kg] and "rhod" [kg m^-3] and gridbox bounds (see add_height).
    """
    numconc = (ds.numconc_aerosol + ds.numconc_cloud + ds.numconc_rain).assign_attrs(
        units="cm^-3", long_name="droplet number concentration"
    )
    lwc = (ds.qcond * ds.rhod).assign_attrs(  # [g/kg] * [kg m^-3] = [g m^-3]
        units="g m^-3",
        long_name="liquid water content",
        description="(excludes solute mass)",
    )
    dz = ds.height_upper - ds.height_lower
    lwp = ((lwc * dz).where(ds.height_lower >= 0.0).sum("height") / 1000).assign_attrs(
        units="kg m^-2",
        long_name="liquid water path (above surface)",
        description="(excludes solute mass)",
    )
    totnsupers = ds.nsupers.sum("height").assign_attrs(
        units="", long_name="total number of superdroplets in domain"
    )

    return ds.assign(numconc=numconc, lwc=lwc, lwp=lwp, totnsupers=totnsupers)


def add_incloud_numconc(ds, cloud_water_radius_range, lwc_threshold=1e-5):
    """Add in-cloud droplet number concentration at the times of the size distributions.

    As in CLEO's post-processing (and Hill et al. 2023): droplets with radius >= 1 micron in
    gridboxes with LWC > 1e-5 kg m^-3 (zero outside cloud). PySDM's number concentration of
    droplets with radius >= 1 micron is the cloud plus rain number concentration.

    *NOTE:* this assumes the setting of the lower bound for cloud droplet radii to
    1 micron in PySDM's runs (see ``cloud_water_radius_range``).

    Args:
        ds (xarray.Dataset): Dataset with "numconc_cloud", "numconc_rain" [cm^-3] and "lwc"
          [g m^-3] (see add_liquid_water) and "time_sds" coordinate.
        cloud_water_radius_range (tuple): [min, max) radius of cloud water [m] (from config);
          its lower bound must be 1 micron.
        lwc_threshold (float): Minimum LWC of a cloudy gridbox [kg m^-3].
    """
    import numpy as np

    cloud_rmin = float(cloud_water_radius_range[0])
    assert np.isclose(cloud_rmin, 1e-6), (
        f"numconc_d assumes cloud droplets start at radius 1micron, not {cloud_rmin}m"
    )
    numconc_d = (ds.numconc_cloud + ds.numconc_rain).where(
        ds.lwc / 1000 > lwc_threshold, 0.0
    )
    numconc_d = _at_time_sds(numconc_d, ds.time_sds)
    return ds.assign(
        numconc_d=numconc_d.assign_attrs(
            units="cm^-3",
            long_name="in-cloud droplet number concentration",
            description="(radius >= 1 micron, LWC > 1e-5 kg m^-3; "
            "zero outside cloud; LWC calculation excludes solute mass)",
        )
    )


def add_surface_precipitation(ds, timestep, rolling_window=None):
    """Add surface precipitation accumulation and (optionally) rolling mean rate.

    PySDM's surface precipitation rate at each output time is the mean over the preceding
    timestep of the volume of water falling through the surface (z = 0) per unit area.

    Args:
        ds (xarray.Dataset): Dataset with "surfprecip_rate" [mm hr^-1].
        timestep (float): Time between outputs (== PySDM timestep) [s].
        rolling_window (int, optional): Width of centred rolling mean window [number of output
          steps]. No rolling mean if None. Note: values within half a window of the
          start/end of the run are NaN.
    """
    import numpy as np

    if "surfprecip_rate" not in ds:
        logger.warning(
            "no 'surfprecip_rate' in dataset, skipping surface precipitation"
        )
        return ds

    assert np.allclose(np.diff(ds.time.values), timestep), (
        "expected PySDM output at every timestep"
    )

    rate = ds.surfprecip_rate.assign_attrs(
        units="mm hr^-1", long_name="surface precipitation rate"
    )
    cumulative = (
        (rate * timestep / 3600)
        .cumsum("time")
        .assign_attrs(units="mm", long_name="accumulated surface precipitation")
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


def add_var_metadata(ds, cloud_water_radius_range, rain_water_radius_range):
    """Add a "long_name" attribute (and a "description" for some variables) to each of
    PySDM's output variables.

    Existing attributes are never overwritten; only missing attributes are added.
    Variables that are not in the dataset are skipped.

    Args:
        ds (xarray.Dataset): Dataset of PySDM output.
        cloud_water_radius_range (tuple): [min, max) radius of cloud water [m].
        rain_water_radius_range (tuple): [min, max) radius of rain water [m].
    """

    def _range(radius_range):
        rmin, rmax = (float(r) * 1e6 for r in radius_range)
        return f"{rmin:g} <= radius < {rmax:g} microns (radius excludes solute volume)"

    cloud, rain = _range(cloud_water_radius_range), _range(rain_water_radius_range)
    aerosol = _range((0.0, cloud_water_radius_range[0]))
    attrs = {
        "qvap": {"long_name": "water vapour mass mixing ratio"},
        "press": {"long_name": "air pressure"},
        "temp": {"long_name": "air temperature"},
        "rhod": {
            "long_name": "density of dry air",
            "description": "(PySDM's dry air density used in the model)",
        },
        "thd": {"long_name": "dry air potential temperature"},
        "rh": {
            "long_name": "relative humidity in PySDM",
            "description": "relative humidity PySDM uses in the model",
        },
        "qcloud": {
            "long_name": "cloud water mass mixing ratio",
            "description": f"({cloud}; excludes solute mass)",
        },
        "qrain": {
            "long_name": "rain water mass mixing ratio",
            "description": f"({rain}; excludes solute mass)",
        },
        "qcond": {
            "long_name": "liquid water mass mixing ratio",
            "description": "liquid water only (excludes solute mass in calculation)",
        },
        "numconc_aerosol": {
            "long_name": "aerosol (unactivated) particle number concentration",
            "description": f"({aerosol})",
        },
        "numconc_cloud": {
            "long_name": "cloud droplet number concentration",
            "description": f"({cloud})",
        },
        "numconc_rain": {
            "long_name": "raindrop number concentration",
            "description": f"({rain})",
        },
        "mean_radius": {
            "long_name": "number weighted mean radius of all particles",
            "description": "(radius excludes solute volume)",
        },
        "effective_radius_cloud": {
            "long_name": "effective radius of cloud droplets",
            "description": f"({cloud})",
        },
        "nsupers": {"long_name": "number of superdroplets in gridbox"},
        "vterm_rain": {
            "long_name": "volume weighted mean terminal velocity of raindrops",
            "description": f"({rain})",
        },
        "peak_saturation": {
            "long_name": "peak saturation ratio during condensation substeps",
            "description": "(maximum since previous output time)",
        },
        "activating_rate": {
            "long_name": "rate of activation events per mass of dry air",
        },
        "deactivating_rate": {
            "long_name": "rate of deactivation events per mass of dry air",
        },
        "ripening_rate": {
            "long_name": "rate of ripening events per mass of dry air",
        },
        "collision_rate": {
            "long_name": "rate of collision events per mass of dry air",
        },
        "collision_deficit_rate": {
            "long_name": "rate of collision deficit per mass of dry air",
            "description": "(collisions that multiplicities could not represent)",
        },
        "coalescence_rate": {
            "long_name": "rate of coalescence events per mass of dry air",
        },
        "numconc_distrib": {
            "long_name": "number concentration in each radius bin",
            "description": "(radius of water volume, i.e. excludes solute volume)",
        },
        "numconc_distrib_dry": {
            "long_name": "number concentration in each dry radius bin",
        },
    }
    return ds.assign(
        {
            # existing attributes come last so they take precedence over the new ones
            name: ds[name].assign_attrs({**var_attrs, **ds[name].attrs})
            for name, var_attrs in attrs.items()
            if name in ds
        }
    )


def add_config_metadata(ds, config_filename):
    """Copy the contents of the PySDM .yaml config file into the dataset's global attributes."""
    from pathlib import Path

    return ds.assign_attrs(
        pysdm_config=Path(config_filename).read_text(),
        pysdm_config_filename=str(config_filename),
    )


def check_radius_ranges(cloud_water_radius_range, rain_water_radius_range):
    """Assert that cloud and rain radius ranges join up and rain is unbounded above.

    Post-processing assumes the radius ranges [0, r_cloud_min), [r_cloud_min, r_rain_min)
    and [r_rain_min, inf) cover all particles exactly once, so that e.g.
    numconc = numconc_aerosol + numconc_cloud + numconc_rain and
    numconc_d = numconc_cloud + numconc_rain (radius >= r_cloud_min).
    """
    import numpy as np

    cloud = [float(r) for r in cloud_water_radius_range]  # [m] (YAML may give strings)
    rain = [float(r) for r in rain_water_radius_range]  # [m]

    assert 0.0 < cloud[0] < cloud[1], f"invalid cloud water radius range {cloud}"
    assert np.isclose(rain[0], cloud[1]), (
        f"rain water radius range must start where cloud water ends ({cloud[1]} m), "
        f"not at {rain[0]} m"
    )
    assert np.isinf(rain[1]) and rain[1] > 0, (
        f"rain water radius range must be unbounded above, not end at {rain[1]} m"
    )


def drop_repeated_variables(ds):
    from sdm_1d_kid import common as cmn

    source2 = "PySDM"
    rtol = 1e-4

    is_rhod_equal = cmn.check_dry_density(ds, logger, source2=source2, rtol=rtol)
    if is_rhod_equal:
        ds = ds.drop_vars(["rhod"], errors="ignore")

    is_rh_equal = cmn.check_relative_humidity(ds, logger, source2=source2, rtol=rtol)
    if is_rh_equal:
        ds = ds.drop_vars(["rh"], errors="ignore")

    is_thd_equal = cmn.check_dry_potential_temperature(
        ds, logger, source2=source2, rtol=rtol
    )
    if is_thd_equal:
        ds = ds.drop_vars(["thd"], errors="ignore")

    return ds


def check_size_distributions(ds, rtol=1e-3):
    """Log whether the sum over radius bins reproduces PySDM's number concentration.

    Number concentration summed over all radius bins should equal "numconc" at the
    sub-sampled times, unless particles lie outside the radius bins. The maximum difference
    is relative to the mean of the non-zero values; differences of order 1e-7 to 1e-5 are
    expected from float32 output.
    """
    import numpy as np

    distrib = "numconc_distrib"
    if distrib not in ds:
        return
    binned = ds[distrib].sum("radius_bin")
    reference = _at_time_sds(ds.numconc, ds.time_sds)
    scale = float(np.abs(reference).where(reference != 0.0).mean())
    if np.isnan(scale):
        logger.info(f"numconc is zero everywhere, no check of {distrib}")
        return
    maxdiff = float(np.abs(binned - reference).max()) / scale
    if maxdiff > rtol:
        logger.warning(
            f"sum over bins of {distrib} differs from numconc by up to {maxdiff:.2e} "
            "(relative to non-zero mean of numconc)"
        )
    else:
        logger.info(f"sum over bins of {distrib} matches numconc ({maxdiff:.2e})")


# %% ----------------------- CREATE ENSEMBLE DATASET ------------------------- #
def create_ensemble_dataset(
    config,
    run_datasets,
    sddiags_every,
    precip_rolling_window=None,
):
    """Combine PySDM output of an ensemble of (one or more) runs into one post-processed
    dataset.

    Args:
        config (dict): PySDM configuration (from .yaml config file).
        run_datasets (list): Paths to the .npy output directory of each ensemble member.
        sddiags_every (int): Size distributions at every 'sddiags_every'th spectrum time.
        precip_rolling_window (int, optional): Width of rolling mean of surface
          precipitation rate [number of output steps]; no rolling mean if None.

    Returns:
        xarray.Dataset with dimensions (ensemble, height, time) and, for superdroplet
        diagnostics, (ensemble, height, time_sds) and (ensemble, height, time_sds,
        radius_bin) (or radius_bin_dry).
    """
    from sdm_1d_kid.cleo import postprocess as cleo_postprocess

    if not run_datasets:
        raise ValueError("Cannot postprocess PySDM runs: no run datasets given")
    run_datasets = [str(d) for d in run_datasets]

    kid_config, outputdata = config["kid_config"], config["outputdata"]
    check_radius_ranges(
        outputdata["cloud_water_radius_range"], outputdata["rain_water_radius_range"]
    )
    return (
        open_ensemble(run_datasets, sddiags_every)
        .pipe(add_height, float(kid_config["zdelta"]))
        .pipe(add_liquid_water)
        .pipe(add_incloud_numconc, outputdata["cloud_water_radius_range"])
        .pipe(
            add_surface_precipitation,
            float(kid_config["timestep"]),
            precip_rolling_window,
        )
        .pipe(cleo_postprocess.add_thermodynamics)
        .pipe(cleo_postprocess.order_dimensions)
        .pipe(
            add_var_metadata,
            outputdata["cloud_water_radius_range"],
            outputdata["rain_water_radius_range"],
        )
        .pipe(cleo_postprocess.metadata_corrections)
    )


# %% -------------------------------- MAIN ----------------------------------- #
def main(
    config_filename,
    run_datasets,
    postprocpath=None,
):
    from pathlib import Path

    import yaml

    ### --- PARAMETERS FOR SURFACE PRECIPITATION --- ###
    precip_rolling_window = 100  # [number of output steps]

    ### --- PARAMETERS FOR SUPERDROPLET DIAGNOSTICS --- ###
    sddiags_every = 10  # size distributions at every 'sddiags_every'th spectrum time

    with open(config_filename) as config_file:
        config = yaml.safe_load(config_file)
    ensemble_name = config["run_ensemble"]["ensemble_name"]

    ds = create_ensemble_dataset(
        config,
        run_datasets,
        sddiags_every,
        precip_rolling_window=precip_rolling_window,
    ).pipe(add_config_metadata, config_filename)

    check_size_distributions(ds, rtol=1e-4)

    ds = drop_repeated_variables(ds)

    if postprocpath is not None:
        postprocpath = Path(postprocpath)
        if not postprocpath.is_dir():
            raise ValueError(
                f"Cannot write post-processed dataset for PySDM: parent directory doesn't exist, {postprocpath}"
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
        args.run_datasets,
        postprocpath=args.postprocpath,
    )
