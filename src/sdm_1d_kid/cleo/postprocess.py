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
into a single dataset.
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


# %% -------------------- SUPERDROPLET DIAGNOSTICS -------------------------- #
def _gridbox_sum(gbxindex, weights, ngbxs):
    """Sum 'weights' of superdroplets in each gridbox."""
    import numpy as np

    return np.bincount(gbxindex, weights=weights, minlength=ngbxs)


def _bin_sum(flat, is_binned, weights, ngbxs, nbins):
    """Sum 'weights' of binned superdroplets into (gridbox, radius bin)."""
    import numpy as np

    return np.bincount(
        flat, weights=weights[is_binned], minlength=ngbxs * nbins
    ).reshape(ngbxs, nbins)


def _superdroplet_diagnostics(
    gbxindex, radius, xi, msol, lwc, ngbxs, consts, radius_edges, dry_radius_edges
):
    """Diagnostics of superdroplets in each gridbox at a single output time.

    Args:
        gbxindex (numpy.ndarray): Gridbox index of each superdroplet (in-domain only).
        radius (numpy.ndarray): Wet radius of each superdroplet [m].
        xi (numpy.ndarray): Multiplicity of each superdroplet.
        msol (numpy.ndarray): Solute mass of each superdroplet [g].
        lwc (numpy.ndarray): Liquid water content of each gridbox [g m^-3] (from massmom1).
        ngbxs (int): Number of gridboxes.
        consts (dict): CLEO constants (from setup .txt file).
        radius_edges (numpy.ndarray): Edges of (wet) radius bins [m].
        dry_radius_edges (numpy.ndarray): Edges of dry radius bins [m].

    Returns:
        dict of arrays with (ngbxs,) or, for size distributions, (ngbxs, nbins) or
        (ngbxs, nbins_dry), with quantities summed over superdroplets (i.e. not yet divided
        by gridbox volume):
          - "xi_incloud": number of droplets with radius > 1 micron in cloudy gridboxes
          - "dvol": mass-weighted mean volume diameter [m]
          - "dvol_sigma": standard deviation of diameter about dvol [m]
          - "numconc", "mass", "reflectivity": per (wet) radius bin [number, g, mm^6]
          - "numconc_dry", "msol_mass": per dry radius bin [number, g]
          - "nsupers_not_binned": number of superdroplets outside the (wet) radius bins
          - "nsupers_not_binned_dry": number of superdroplets outside the dry radius bins
    """
    import numpy as np

    rho_l, rho_sol = consts["RHO_L"], consts["RHO_SOL"]  # [kg m^-3]

    # mass of water in each superdroplet [kg] (volume of droplet minus volume of solute)
    m_water = rho_l * (4.0 / 3.0 * np.pi * radius**3 - msol / 1000 / rho_sol)
    m_water = np.where(m_water <= 0.0, 0.0, m_water)
    diameter = (6.0 / np.pi / rho_l * m_water) ** (
        1 / 3
    )  # volume diameter of water [m]

    # in-cloud droplet number: droplets with radius > 1 micron in gridboxes with
    # LWC > 1e-5 kg m^-3 (zero outside cloud), as in Hill et al. (2023)
    is_cloudy = lwc / 1000 > 1e-5
    xi_incloud = _gridbox_sum(gbxindex, np.where(radius > 1e-6, xi, 0.0), ngbxs)
    xi_incloud = np.where(is_cloudy, xi_incloud, 0.0)

    # mass-weighted mean volume diameter of all superdroplets,
    # D_vol = sum(xi m D) / sum(xi m) (NaN if no water)
    sum_xi_m = _gridbox_sum(gbxindex, xi * m_water, ngbxs)
    sum_xi_m_d = _gridbox_sum(gbxindex, xi * m_water * diameter, ngbxs)
    dvol = np.full(ngbxs, np.nan)
    np.divide(sum_xi_m_d, sum_xi_m, out=dvol, where=sum_xi_m > 0.0)

    # standard deviation of diameter of all superdroplets about D_vol
    # (number weighted, with Bessel's correction)
    sum_xi = _gridbox_sum(gbxindex, xi, ngbxs)
    sum_xi_dev2 = _gridbox_sum(gbxindex, xi * (diameter - dvol[gbxindex]) ** 2, ngbxs)
    dvol_sigma = np.full(ngbxs, np.nan)
    np.divide(sum_xi_dev2, sum_xi - 1, out=dvol_sigma, where=sum_xi > 1)
    dvol_sigma = np.sqrt(dvol_sigma)

    # size distributions binned by wet radius
    nbins = len(radius_edges) - 1
    ibin = np.searchsorted(radius_edges, radius, side="right") - 1
    is_binned = (ibin >= 0) & (ibin < nbins)
    flat = gbxindex[is_binned] * nbins + ibin[is_binned]

    # size distributions binned by dry radius, i.e. radius of a sphere of the solute,
    # r_dry = (3 msol / (4 pi rho_sol))^(1/3) (as PySDM's dry radius)
    dry_radius = (3.0 * msol / 1000 / (4.0 * np.pi * rho_sol)) ** (1 / 3)  # [m]
    nbins_dry = len(dry_radius_edges) - 1
    ibin_dry = np.searchsorted(dry_radius_edges, dry_radius, side="right") - 1
    is_binned_dry = (ibin_dry >= 0) & (ibin_dry < nbins_dry)
    flat_dry = gbxindex[is_binned_dry] * nbins_dry + ibin_dry[is_binned_dry]

    # droplet mass [g] as in CLEO's SuperdropAttrs::mass(), i.e. including the solute
    mass = 4.0 / 3.0 * np.pi * rho_l * 1000 * radius**3 + (1.0 - rho_l / rho_sol) * msol

    return {
        "xi_incloud": xi_incloud,
        "dvol": dvol,
        "dvol_sigma": dvol_sigma,
        "numconc": _bin_sum(flat, is_binned, xi, ngbxs, nbins),
        "mass": _bin_sum(flat, is_binned, xi * mass, ngbxs, nbins),
        "reflectivity": _bin_sum(
            flat, is_binned, xi * (2 * radius * 1000) ** 6, ngbxs, nbins
        ),
        "numconc_dry": _bin_sum(flat_dry, is_binned_dry, xi, ngbxs, nbins_dry),
        "msol_mass": _bin_sum(flat_dry, is_binned_dry, xi * msol, ngbxs, nbins_dry),
        "nsupers_not_binned": np.count_nonzero(~is_binned),
        "nsupers_not_binned_dry": np.count_nonzero(~is_binned_dry),
    }


def _radius_bin_coords(edges, dim, long_name_prefix="", description=""):
    """Coordinates of radius bins along dimension 'dim': geometric mean of bin bounds
    ('dim'), bounds ('dim'_lower, 'dim'_upper) and width in ln(radius) ('dim'_dlnr)."""
    import numpy as np

    lower, upper = edges[:-1], edges[1:]
    return {
        dim: (
            dim,
            np.sqrt(lower * upper),
            {
                "units": "m",
                "long_name": f"{long_name_prefix}radius of bin centre",
                "description": "centre is geometric mean of bounds" + description,
            },
        ),
        f"{dim}_lower": (
            dim,
            lower,
            {"units": "m", "long_name": f"lower bound of {long_name_prefix}radius bin"},
        ),
        f"{dim}_upper": (
            dim,
            upper,
            {"units": "m", "long_name": f"upper bound of {long_name_prefix}radius bin"},
        ),
        f"{dim}_dlnr": (
            dim,
            np.log(upper / lower),
            {
                "units": "",
                "long_name": f"width of {long_name_prefix}radius bin in ln(radius)",
                "description": "width = ln|upper| - ln|lower|",
            },
        ),
    }


def process_superdroplets(
    ds, consts, gbxvols, radius_edges, dry_radius_edges, sddiags_every
):
    """Superdroplet diagnostics and size distributions of a single run.

    All superdroplet diagnostics are calculated at every 'sddiags_every'th output time,
    i.e. at the times of the "time_sds" coordinate. Returns 'ds' with the following variables
    added:
      with dimensions (time_sds, gbxindex):
        - numconc_d: in-cloud droplet number concentration [cm^-3], counting droplets with
          radius > 1 micron in gridboxes with LWC > 1e-5 kg m^-3 (zero outside cloud)
        - dvol: mass-weighted mean volume diameter of all superdroplets [um]
        - dvol_sigma: standard deviation of diameter of all superdroplets about dvol [um]
      with dimensions (time_sds, gbxindex, radius_bin):
        - numconc_distrib: number concentration [cm^-3] (sum over bins = massmom0 / volume)
        - mass_distrib: droplet mass concentration [g m^-3] (sum over bins = massmom1 / volume)
        - reflectivity_distrib: radar reflectivity factor [mm^6 m^-3]
      with dimensions (time_sds, gbxindex, radius_bin_dry):
        - numconc_distrib_dry: number concentration [cm^-3] (sum over bins = massmom0 / volume)
        - msol_mass_distrib: solute mass concentration [g m^-3]

    Droplet water mass is m_w = rho_l * (4/3 pi r^3 - msol / rho_sol) and its volume diameter
    is D = (6 m_w / (pi rho_l))^(1/3). Droplet mass in mass_distrib is calculated as in CLEO's
    SuperdropAttrs::mass() (i.e. including the solute), so that the sum over all bins equals
    CLEO's massmom1. The radar reflectivity factor is Z = sum(xi * (2r)^6) / volume, which
    assumes liquid spheres much smaller than the radar wavelength (Rayleigh regime). The dry
    radius is the radius of a sphere of the solute, r_dry = (3 msol / (4 pi rho_sol))^(1/3).

    Superdroplets out of the domain are excluded from all diagnostics and superdroplets
    outside the (wet or dry) radius bins are excluded from those size distributions (a
    warning is logged if there are any).

    Args:
        ds (xarray.Dataset): CLEO output of a single run with superdroplet ragged arrays.
        consts (dict): CLEO constants (from setup .txt file).
        gbxvols (numpy.ndarray): Volume of each gridbox [m^3], ordered by gbxindex.
        radius_edges (numpy.ndarray): Edges of (wet) radius bins [m], monotonically increasing.
        dry_radius_edges (numpy.ndarray): Edges of dry radius bins [m], monotonically
          increasing.
        sddiags_every (int): Superdroplet diagnostics at every 'sddiags_every'th output time.
    """
    import numpy as np

    ngbxs = ds.sizes["gbxindex"]
    assert np.array_equal(ds.gbxindex.values, np.arange(ngbxs)), (
        "expected gbxindex to be 0, 1, ..., ngbxs-1"
    )
    nbins, nbins_dry = len(radius_edges) - 1, len(dry_radius_edges) - 1
    vol = np.asarray(gbxvols)[np.newaxis, :]  # [m^3]
    lwc = ds.massmom1.values / vol  # [g m^-3]

    # start/end of each output time in the contiguous ragged arrays
    offsets = np.concatenate([[0], np.cumsum(ds.raggedcount.values.astype(np.int64))])

    itimes = np.arange(0, ds.sizes["time"], sddiags_every)
    shape = (len(itimes), ngbxs)
    shape_bins, shape_bins_dry = (*shape, nbins), (*shape, nbins_dry)
    diags = {name: np.zeros(shape) for name in ("xi_incloud", "dvol", "dvol_sigma")}
    diags.update(
        {name: np.zeros(shape_bins) for name in ("numconc", "mass", "reflectivity")}
    )
    diags.update(
        {name: np.zeros(shape_bins_dry) for name in ("numconc_dry", "msol_mass")}
    )
    nsupers_out_of_domain, nsupers_not_binned, nsupers_not_binned_dry = 0, 0, 0

    for n, it in enumerate(itimes):
        supers = (
            ds[["sdgbxindex", "radius", "xi", "msol"]]
            .isel(superdroplets=slice(offsets[it], offsets[it + 1]))
            .load()
        )
        gbxindex = supers.sdgbxindex.values.astype(np.int64)
        in_domain = (gbxindex >= 0) & (gbxindex < ngbxs)
        nsupers_out_of_domain += np.count_nonzero(~in_domain)

        diags_n = _superdroplet_diagnostics(
            gbxindex[in_domain],
            supers.radius.values[in_domain].astype(np.float64) * 1e-6,  # [m]
            supers.xi.values[in_domain].astype(np.float64),
            supers.msol.values[in_domain].astype(np.float64),  # [g]
            lwc[it],
            ngbxs,
            consts,
            radius_edges,
            dry_radius_edges,
        )
        for name, value in diags.items():
            value[n] = diags_n[name]
        nsupers_not_binned += diags_n["nsupers_not_binned"]
        nsupers_not_binned_dry += diags_n["nsupers_not_binned_dry"]

    if nsupers_out_of_domain:
        logger.warning(
            f"{nsupers_out_of_domain} superdroplet records (summed over sub-sampled times) "
            "are out of the domain and are excluded from the superdroplet diagnostics"
        )
    if nsupers_not_binned:
        logger.warning(
            f"{nsupers_not_binned} superdroplet records (summed over sub-sampled times) are "
            "outside the radius bins and are excluded from the size distributions"
        )
    if nsupers_not_binned_dry:
        logger.warning(
            f"{nsupers_not_binned_dry} superdroplet records (summed over sub-sampled times) "
            "are outside the dry radius bins and are excluded from the dry size distributions"
        )

    vol_bins = vol[:, :, np.newaxis]  # [m^3]
    dims = ("time_sds", "gbxindex")
    dims_bins = ("time_sds", "gbxindex", "radius_bin")
    dims_bins_dry = ("time_sds", "gbxindex", "radius_bin_dry")
    coords = {
        "time_sds": (
            "time_sds",
            ds.time.values[itimes],
            {
                "units": ds.time.attrs.get("units", "s"),
                "long_name": "time of superdroplet diagnostics",
            },
        ),
        **_radius_bin_coords(
            radius_edges, "radius_bin", description="; radius includes solute volume"
        ),
        **_radius_bin_coords(
            dry_radius_edges,
            "radius_bin_dry",
            "dry ",
            description="; dry radius of the particle (solute)",
        ),
    }

    return ds.assign_coords(coords).assign(
        numconc_d=(
            dims,
            diags["xi_incloud"] / vol / 1e6,
            {
                "units": "cm^-3",
                "long_name": "in-cloud droplet number concentration",
                "description": "(radius > 1 micron, LWC > 1e-5 kg m^-3; "
                "zero outside cloud; LWC calculation includes droplet solute mass)",
            },
        ),
        dvol=(
            dims,
            diags["dvol"] * 1e6,
            {
                "units": "um",
                "long_name": "mass weighted mean volume diameter",
                "description": "(excludes solute mass in calculation)",
            },
        ),
        dvol_sigma=(
            dims,
            diags["dvol_sigma"] * 1e6,
            {
                "units": "um",
                "long_name": "standard deviation of volume diameter about dvol",
                "description": "(excludes solute mass in calculation)",
            },
        ),
        numconc_distrib=(
            dims_bins,
            (diags["numconc"] / vol_bins / 1e6).astype(np.float32),
            {
                "units": "cm^-3",
                "long_name": "number concentration in each radius bin",
                "description": "(calculation bins by radius which includes solute volume)",
            },
        ),
        mass_distrib=(
            dims_bins,
            (diags["mass"] / vol_bins).astype(np.float32),
            {
                "units": "g m^-3",
                "long_name": "droplet mass concentration in each radius bin",
                "description": "(includes droplet solute mass)",
            },
        ),
        reflectivity_distrib=(
            dims_bins,
            (diags["reflectivity"] / vol_bins).astype(np.float32),
            {
                "units": "mm^6 m^-3",
                "long_name": "radar reflectivity factor in each radius bin (Rayleigh regime)",
                "description": "(radius calculation includes droplet solute mass)",
            },
        ),
        numconc_distrib_dry=(
            dims_bins_dry,
            (diags["numconc_dry"] / vol_bins / 1e6).astype(np.float32),
            {
                "units": "cm^-3",
                "long_name": "number concentration in each dry radius bin",
            },
        ),
        msol_mass_distrib=(
            dims_bins_dry,
            (diags["msol_mass"] / vol_bins).astype(np.float32),
            {
                "units": "g m^-3",
                "long_name": "solute mass concentration in each dry radius bin",
            },
        ),
    )


# %% ------------------------ OPEN ENSEMBLE OF RUNS -------------------------- #
def open_ensemble(
    run_datasets, consts, gbxs, radius_edges, dry_radius_edges, sddiags_every
):
    """Open CLEO datasets of each ensemble member as one dataset with an "ensemble" dimension.

    Superdroplet (ragged array) data are processed into diagnostics and size distributions
    (see process_superdroplets) and then dropped, because their length differs between
    members.
    All members must have identical output times (join="exact" raises an error otherwise).

    Args:
        run_datasets (list[str]): Paths to the "*.zarr" dataset of each ensemble member.
        consts (dict): CLEO constants (from setup .txt file).
        gbxs (dict): CLEO gridboxes (from cleopy's pygbxsdat.get_gridboxes).
        radius_edges (numpy.ndarray): Edges of (wet) radius bins [m].
        dry_radius_edges (numpy.ndarray): Edges of dry radius bins [m].
        sddiags_every (int): Superdroplet diagnostics at every 'sddiags_every'th output time.

    Returns:
        xarray.Dataset with dimensions (ensemble, time, gbxindex) and, for superdroplet
        diagnostics, (ensemble, time_sds, gbxindex) and for superdroplet distributions
        (ensemble, time_sds, gbxindex, radius_bin) or (ensemble, time_sds, gbxindex,
        radius_bin_dry).
    """
    from pathlib import Path

    import xarray as xr

    gbxvols = gbxs["gbxvols"][0, 0, :]

    def preprocess(ds):
        ds = process_superdroplets(
            ds, consts, gbxvols, radius_edges, dry_radius_edges, sddiags_every
        )
        return ds.drop_dims("superdroplets", errors="ignore").drop_vars(
            "raggedcount", errors="ignore"
        )

    ds = xr.open_mfdataset(
        run_datasets,
        engine="zarr",
        combine="nested",
        concat_dim="ensemble",
        preprocess=preprocess,
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
        units="g m^-3",
        long_name="liquid water content",
        description="(includes droplet solute mass)",
    )
    dz = ds.height_upper - ds.height_lower
    lwp = ((lwc * dz).where(ds.height_lower >= 0.0).sum("height") / 1000).assign_attrs(
        units="kg m^-2",
        long_name="liquid water path (above surface)",
        description="(includes droplet solute mass)",
    )

    return ds.assign(numconc=numconc, lwc=lwc, lwp=lwp)


def add_radar_reflectivity(ds):
    """Add total radar reflectivity (sum over radius bins) in dBZ.

    Gridboxes with zero reflectivity (no droplets) are NaN.

    Args:
        ds (xarray.Dataset): Dataset with "reflectivity_distrib" (see process_superdroplets).
    """
    import numpy as np

    reflectivity = ds.reflectivity_distrib.sum("radius_bin")
    dbz = 10 * np.log10(
        reflectivity.where(reflectivity > 0.0)
    )  # Z relative to 1 mm^6 m^-3
    return ds.assign(
        reflectivity_dbz=dbz.assign_attrs(
            units="dBZ",
            long_name="radar reflectivity (Rayleigh regime)",
            description="(radius calculation includes droplet solute mass)",
        )
    )


def add_surface_precipitation(ds, obststep, rolling_window=None):
    """Convert precipitation to [mm] and add surface precipitation rate, accumulation and
    (optionally) rolling mean rate.

    CLEO's "precip" is the depth of water [m] moving downwards through the lower boundary of
    each gridbox per output interval. It is converted to [mm] and surface precipitation is
    taken from the gridbox whose lower boundary is z = 0.

    Args:
        ds (xarray.Dataset): Dataset with "precip" [m] and gridbox bounds (see add_height).
        obststep (float): Time between outputs [s].
        rolling_window (int, optional): Width of centred rolling mean window [number of output
          steps]. No rolling mean if None. Note: values within half a window of the
          start/end of the run are NaN.
    """
    import numpy as np

    if "precip" not in ds:
        logger.warning("no 'precip' in dataset, skipping surface precipitation")
        return ds

    precip = (ds.precip * 1000).assign_attrs(
        units="mm",
        long_name="depth of water falling through gridbox lower boundary per output interval",
    )
    ds = ds.assign(precip=precip)

    is_surface = np.isclose(ds.height_lower.values, 0.0)
    assert is_surface.sum() == 1, (
        "expected exactly one gridbox with lower boundary at z = 0"
    )
    depth = ds.precip.isel(
        height=int(np.argmax(is_surface))
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


def add_var_metadata(ds):
    """Add a "long_name" attribute (and a "description" for some variables) to each of CLEO's
    output variables.

    Existing attributes (e.g. "units", or a "long_name" already written by CLEO) are never
    overwritten; only missing attributes are added. Raindrops are droplets with
    radius >= 40 microns (see CLEO's libs/observers/massmoments_observer.hpp). Variables that
    are not in the dataset are skipped.
    """
    attrs = {
        "massmom0": {
            "long_name": "0th mass moment of droplet distribution",
            "description": "number of droplets in gridbox",
        },
        "massmom0_raindrops": {
            "long_name": "0th mass moment of raindrop distribution",
            "description": "number of raindrops in gridbox (radius >= 40 microns)",
        },
        "massmom1": {
            "long_name": "1st mass moment of droplet distribution",
            "description": "mass of droplets in gridbox (includes droplet solute mass)",
        },
        "massmom1_raindrops": {
            "long_name": "1st mass moment of raindrop distribution",
            "description": "mass of raindrops in gridbox (radius >= 40 microns; "
            "includes droplet solute mass)",
        },
        "massmom2": {
            "long_name": "2nd mass moment of droplet distribution",
            "description": "proxy for radar reflectivity (Rayleigh regime; "
            "includes droplet solute mass)",
        },
        "massmom2_raindrops": {
            "long_name": "2nd mass moment of raindrop distribution",
            "description": "proxy for radar reflectivity of raindrops "
            "(Rayleigh regime; radius >= 40 microns; includes droplet solute mass)",
        },
        "nsupers": {"long_name": "number of superdroplets in gridbox"},
        "totnsupers": {"long_name": "total number of superdroplets in domain"},
        "press": {"long_name": "air pressure"},
        "qcond": {
            "long_name": "liquid water mass mixing ratio",
            "description": "liquid water only (excludes solute mass in calculation)",
        },
        "qvap": {"long_name": "water vapour mass mixing ratio"},
        "temp": {"long_name": "air temperature"},
        "uvel": {"long_name": "eastward (x) wind velocity"},
        "vvel": {"long_name": "northward (y) wind velocity"},
        "wvel": {"long_name": "vertical (z) wind velocity"},
    }
    return ds.assign(
        {
            # existing attributes come last so they take precedence over the new ones
            name: ds[name].assign_attrs({**var_attrs, **ds[name].attrs})
            for name, var_attrs in attrs.items()
            if name in ds
        }
    )


def add_setup_metadata(ds, setup_filename):
    """Copy the contents of CLEO's setup .txt file into the dataset's global attributes."""
    from pathlib import Path

    return ds.assign_attrs(
        cleo_setup=Path(setup_filename).read_text(),
        cleo_setup_filename=str(setup_filename),
    )


def check_size_distributions(ds, rtol=1e-3):
    """Log whether the sums over radius bins reproduce CLEO's mass moments.

    Number concentration and mass summed over all radius bins should equal "numconc" and
    "lwc" (from massmom0 and massmom1) at the sub-sampled times, unless superdroplets lie
    outside the radius bins. Likewise number concentration summed over all dry radius bins
    should equal "numconc". The maximum difference is relative to the mean of the
    non-zero values; differences of order 1e-7 to 1e-5 are expected from float32 output.
    """
    import numpy as np

    checks = (
        ("numconc_distrib", "numconc", "radius_bin"),
        ("mass_distrib", "lwc", "radius_bin"),
        ("numconc_distrib_dry", "numconc", "radius_bin_dry"),
    )
    for distrib, total, bin_dim in checks:
        binned = ds[distrib].sum(bin_dim)
        reference = ds[total].sel(time=ds.time_sds)
        scale = float(np.abs(reference).where(reference != 0.0).mean())
        if np.isnan(scale):
            logger.info(f"{total} is zero everywhere, no check of {distrib}")
            continue
        maxdiff = float(np.abs(binned - reference).max()) / scale
        if maxdiff > rtol:
            logger.warning(
                f"sum over bins of {distrib} differs from {total} by up to {maxdiff:.2e} "
                f"(relative to non-zero mean of {total})"
            )
        else:
            logger.info(f"sum over bins of {distrib} matches {total} ({maxdiff:.2e})")


# %% ----------------------- CREATE ENSEMBLE DATASET ------------------------- #
def create_ensemble_dataset(
    config,
    consts,
    gbxs,
    run_datasets,
    radius_edges,
    dry_radius_edges,
    sddiags_every,
    precip_rolling_window=None,
):
    """Combine CLEO output of an ensemble of (one or more) runs into one post-processed dataset.

    Args:
        config (dict): CLEO configuration (from setup .txt file).
        consts (dict): CLEO constants (from setup .txt file).
        gbxs (dict): CLEO gridboxes (from grid file).
        run_datasets (list): Paths to the ".zarr" dataset of each ensemble member.
        radius_edges (numpy.ndarray): Edges of (wet) radius bins for size distributions [m].
        dry_radius_edges (numpy.ndarray): Edges of dry radius bins for size distributions
          [m].
        sddiags_every (int): Superdroplet diagnostics at every 'sddiags_every'th output time.
        precip_rolling_window (int, optional): Width of rolling mean of surface
          precipitation rate [number of output steps]; no rolling mean if None.

    Returns:
        xarray.Dataset with dimensions (ensemble, height, time) and, for superdroplet
        diagnostics, (ensemble, height, time_sds) and (ensemble, height, time_sds,
        radius_bin) (or radius_bin_dry).
    """
    import sdm_1d_kid.common as cmn

    if not run_datasets:
        raise ValueError("Cannot postprocess Cleo runs: no run datasets given")
    run_datasets = sorted(str(d) for d in run_datasets)

    return (
        open_ensemble(
            run_datasets, consts, gbxs, radius_edges, dry_radius_edges, sddiags_every
        )
        .pipe(add_height, gbxs)
        .pipe(add_liquid_water)
        .pipe(add_radar_reflectivity)
        .pipe(add_surface_precipitation, config["OBSTSTEP"], precip_rolling_window)
        .pipe(cmn.add_thermodynamics)
        .pipe(cmn.order_dimensions)
        .pipe(add_var_metadata)
        .pipe(cmn.metadata_corrections)
    )


# %% -------------------------------- MAIN ----------------------------------- #
def main(
    config_filename,
    setup_filename,
    run_datasets,
    postprocpath=None,
):
    from pathlib import Path

    import numpy as np
    import yaml
    from cleopy.sdmout_src import pygbxsdat, pysetuptxt

    ### --- PARAMETERS FOR SURFACE PRECIPITATION --- ###
    precip_rolling_window = 100  # [number of output steps]

    ### --- PARAMETERS FOR SUPERDROPLET DIAGNOSTICS --- ###
    sddiags_every = (
        10  # all superdroplet diagnostics at every 'sddiags_every'th output time
    )
    radius_min = 5.0e-9  # lower edge of smallest radius bin [m]
    radius_max = 5.0e-3  # upper edge of largest radius bin [m]
    nbins = 128  # number of (log-spaced) radius bins
    dry_radius_min = 1.0e-9  # lower edge of smallest dry radius bin [m]
    dry_radius_max = 1.0e-4  # upper edge of largest dry radius bin [m]
    nbins_dry = 128  # number of (log-spaced) dry radius bins

    with open(config_filename) as config_file:
        cnfg = yaml.safe_load(config_file)
    grid_filename = Path(cnfg["inputfiles"]["grid_filename"]).resolve()
    ensemble_name = cnfg["run_ensemble"]["ensemble_name"]

    config = pysetuptxt.get_config(setup_filename, nattrs=3, isprint=False)
    consts = pysetuptxt.get_consts(setup_filename, isprint=False)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=False)
    radius_edges = np.geomspace(radius_min, radius_max, nbins + 1)  # [m]
    dry_radius_edges = np.geomspace(
        dry_radius_min, dry_radius_max, nbins_dry + 1
    )  # [m]

    ds = create_ensemble_dataset(
        config,
        consts,
        gbxs,
        run_datasets,
        radius_edges,
        dry_radius_edges,
        sddiags_every,
        precip_rolling_window=precip_rolling_window,
    ).pipe(add_setup_metadata, setup_filename)

    check_size_distributions(ds, rtol=1e-4)

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
