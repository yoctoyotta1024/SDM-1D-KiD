"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: products.py
Project: helpers
Created Date: Wednesday 25th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Helper file for extra PySDM products, so that PySDM's output can be post-processed
into the same superdroplet diagnostics as CLEO's (mass weighted mean volume diameter,
droplet mass and radar reflectivity distributions).
"""


def make_radius_moment_product(rank):
    """Product class for sum over all particles in each gridbox of xi * r^rank [m^rank],
    where r is PySDM's wet radius (of the particle's water volume)."""
    from PySDM.products.size_spectral.arbitrary_moment import (
        make_arbitrary_moment_product,
    )

    return make_arbitrary_moment_product(
        rank=rank,
        attr="radius",
        attr_unit="m",
        skip_division_by_m0=True,  # sum (not mean) over particles
        skip_division_by_dv=True,  # per gridbox (not per volume of air)
    )


def make_sum_per_radius_bin_product(attr, rank, unit):
    """Product class for sum over particles in each wet radius bin of xi * attr^rank per
    volume of air, i.e. sum(xi * attr^rank) / dv.

    Particles are binned by wet radius (of the particle's water volume), with the bins
    filtered by "volume" as in PySDM's ParticleSizeSpectrumPerVolume. 'unit' must be the
    SI unit of attr^rank per m^3 (PySDM requires each product class to have a default SI
    unit, hence the factory).
    """
    import numpy as np
    from PySDM.products.impl import SpectrumMomentProduct, register_product

    @register_product()
    class SumPerRadiusBin(SpectrumMomentProduct):
        def __init__(self, radius_bins_edges, name=None, unit=unit):
            super().__init__(name=name, unit=unit, attr_unit="m^3")
            self.radius_bins_edges = radius_bins_edges

        def register(self, particulator):
            particulator.request_attribute("volume")
            particulator.request_attribute(attr)
            volume_bins_edges = particulator.formulae.trivia.volume(
                self.radius_bins_edges
            )
            self.attr_bins_edges = particulator.backend.Storage.from_ndarray(
                volume_bins_edges
            )
            super().register(particulator)
            self.shape = (*particulator.mesh.grid, len(self.attr_bins_edges) - 1)

        def _impl(self, **kwargs):
            vals = np.empty(
                [self.particulator.mesh.n_cell, len(self.attr_bins_edges) - 1]
            )
            self._recalculate_spectrum_moment(
                attr=attr, rank=rank, filter_attr="volume"
            )
            for i in range(vals.shape[1]):
                # mean of attr^rank in bin (weighted by multiplicity) ...
                self._download_spectrum_moment_to_buffer(rank=rank, bin_number=i)
                vals[:, i] = self.buffer.ravel()
                # ... times sum of multiplicities in bin
                self._download_spectrum_moment_to_buffer(rank=0, bin_number=i)
                vals[:, i] *= self.buffer.ravel()
            return vals / self.particulator.mesh.dv

    return SumPerRadiusBin


# sum over particles in each gridbox of xi * r^k for k = 0, 1, 2, 3, 4 (for D_vol and sigma)
RadiusMoments = {rank: make_radius_moment_product(rank) for rank in range(5)}

# sums over particles in each wet radius bin per volume of air
WaterMassPerRadiusBin = make_sum_per_radius_bin_product(
    attr="water mass", rank=1, unit="kg / m**3"
)
DryVolumePerRadiusBin = make_sum_per_radius_bin_product(
    attr="dry volume", rank=1, unit="m**3 / m**3"
)
RadiusSixthPowerPerRadiusBin = make_sum_per_radius_bin_product(
    attr="radius", rank=6, unit="m**6 / m**3"
)
