"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: kid_dynamics.py
Project: pympdata_driver
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
class for driving KiD rainshaft test case
"""

import numpy as np
from PyMPDATA_examples.Shipway_and_Hill_2012 import si

import sdm_1d_kid.pympdata_driver.common as cmn
from sdm_1d_kid.pympdata_driver.mpdata import MPDATA
from sdm_1d_kid.pympdata_driver.settings import Settings


class KiDDynamics:
    """A class for driving the KiD rainshaft test case, based on Shipway and Hill 2012.

    This class is a wrapper around the MPDATA driver of the Shipway and Hill 2012 example
    in the PyMPDATA-examples library.

    See https://github.com/open-atmos/PyMPDATA/tree/main/examples/PyMPDATA_examples/Shipway_and_Hill_2012
    for the original source code.
    """

    def __init__(self, kid_config, advect_hydrometeors=True):
        """Initialize the KiDDynamics object.

        Args:
            kid_config (dict): configuration parameters for settings of 1-D KiD.
        """
        ### read in config
        # ``wmax_const`` is similar to 'w1' of equation (6) in Shipway and Hill (2012), but accounting for density
        # ``tscale_const`` is  't1' of equation (6) in Shipway and Hill (2012)
        wmax_const = kid_config[
            "wmax_const"
        ]  # circa. max vertical velocity (is maximum w*rho_dry [kg m^-2 s^-1])
        tscale_const = kid_config["tscale_const"]  # timescale of sinusoid [s]
        is_exner_novapour = kid_config["is_exner_novapour"]
        is_exner_novapour_uniformrho = kid_config["is_exner_novapour_uniformrho"]
        is_approx_drhod_dz = kid_config["is_approx_drhod_dz"]

        ### add SI units to config variables
        wmax_const = wmax_const * si.m / si.s * si.kg / si.m**3
        tscale_const = tscale_const * si.s
        p_surf = kid_config["psurf_hpa"] * si.hPa
        t_max = kid_config["t_end"] * si.s
        dt = kid_config["timestep"] * si.s
        z_min = kid_config["zmin"] * si.m
        z_max = kid_config["zmax"] * si.m
        dz = kid_config["zdelta"] * si.m

        # *NOTE:*(!)must match above SI units^
        self.written_kid_config_units = {
            "zmin": "m",
            "zmax": "m",
            "zdelta": "m",
            "timestep": "s",
            "t_end": "s",
            "wmax_const": "kg m^-2 s^-1 (peak rho_d*w)",
            "tscale_const": "s",
            "psurf_hpa": "hPa",
        }

        ### initialise options, settings and MPDATA
        self.settings = Settings(
            dt=dt,
            dz=dz,
            wmax_const=wmax_const,
            tscale_const=tscale_const,
            t_max=t_max,
            z_min=z_min,
            z_max=z_max,
            p_surf=p_surf,
            is_exner_novapour=is_exner_novapour,
            is_exner_novapour_uniformrho=is_exner_novapour_uniformrho,
            is_approx_drhod_dz=is_approx_drhod_dz,
        )

        options = cmn.get_options()
        self.mpdata = MPDATA(
            nz=self.settings.nz,
            dt=self.settings.dt,
            qv_of_zZ_at_t0=lambda zZ: self.settings.qv0(
                zZ * self.settings.dz + self.settings.z_min
            ),
            g_factor_of_zZ=lambda zZ: self.settings.rhod(
                zZ * self.settings.dz + self.settings.z_min
            ),
            options=options,
        )

        self.advect_hydrometeors = advect_hydrometeors

        assert self.settings.nz == int((z_max - z_min) / dz)
        assert self.settings.z_min == z_min
        assert self.settings.z_max == z_max
        assert self.settings.dz == dz
        self.zfull = np.arange(z_min + dz / 2, z_max + dz / 2, dz)
        self.zhalf = np.arange(z_min, z_max + dz, dz)
        assert np.all((self.zhalf[1:] + self.zhalf[:-1]) / 2 == self.zfull)

        ### print configuration
        key = f"nr={self.mpdata.nr}, dz={self.settings.dz}, dt={self.settings.dt}, options={options}"
        print(f"Simulating {self.settings.nt} timesteps using {key}")
        print(
            self._format_kid_config(
                kid_config, nz=self.settings.nz, nt=self.settings.nt
            )
        )

    def _format_kid_config(self, kid_config, nz=None, nt=None):
        """Return KiD configuration as aligned 'name = value [units]' lines for printing."""
        width = max(len(key) for key in kid_config)
        lines = ["KiDDynamics configuration:"]
        for key, value in kid_config.items():
            units = self.written_kid_config_units.get(key)
            lines.append(
                f"  {key:<{width}} = {value}" + (f" [{units}]" if units else "")
            )
        if nz is not None and nt is not None:
            lines.append(f"  {'(derived)':<{width}} = {nz} gridboxes, {nt} timesteps")
        return "\n".join(lines)

    def set_thermo(self, time, thermo):
        """Set thermodynamics from the dynamics solver.

        Args:
            thermo (Thermodynamics): Object representing the thermodynamic state.

        Returns:
            Thermodynamics: Updated thermodynamic state.
        """
        thermo.massmix_ratios["qvap"][:] = self.mpdata["qvap"].advectee.get()
        thermo.massmix_ratios["qcond"][:] = self.mpdata["qcond"].advectee.get()
        thermo.massmix_ratios["qice"][:] = self.mpdata["qice"].advectee.get()
        thermo.massmix_ratios["qrain"][:] = self.mpdata["qrain"].advectee.get()
        thermo.massmix_ratios["qsnow"][:] = self.mpdata["qsnow"].advectee.get()
        thermo.massmix_ratios["qgrau"][:] = self.mpdata["qgrau"].advectee.get()
        thermo.rhod[:] = self.settings.rhod(self.zfull)
        thermo.temp[:] = self.settings.temp(self.zfull)
        thermo.press[:] = self.settings.press(self.zfull, thermo.massmix_ratios["qvap"])

        if thermo.wvel.size != 0:
            z_half_reps = np.repeat(
                np.arange(
                    self.settings.z_min,
                    self.settings.z_max + self.settings.dz,
                    self.settings.dz,
                ),
                2,
            )[1:-1]
            wmagnitude = self.settings.rhod_w(
                time + self.settings.dt / 2
            ) / self.settings.rhod(z_half_reps)
            assert np.ones_like(thermo.wvel).shape == wmagnitude.shape
            thermo.wvel[:] = wmagnitude

        return thermo

    def run(self, time, timestep, thermo):
        """
        Run the 1-D KiD motion computations.

        This method integrates the equations from time to time+timestep
        for the dynamics of a 1-D KiD rainshaft.

        Args:
            time (float): Current time [s].
            timestep (float): Size of timestep for the simulation [s].
            thermo (Thermodynamics): Object representing the thermodynamic state.

        Returns:
            Thermodynamics: Updated thermodynamic state.
        """
        assert timestep == self.settings.dt, "Timestep must match initialised value."
        assert time % timestep == 0, "Time not a multiple of the timestep."

        GC = (
            self.settings.rhod_w(time + self.settings.dt / 2)
            * self.settings.dt
            / self.settings.dz
        )  # courant_number * velocity * density (i.e. momentum * dt/dz)
        advector_0 = np.ones_like(self.settings.nz_vec) * GC

        if self.advect_hydrometeors:
            for field in self.mpdata.fields:
                self.mpdata[field].advector.get_component(0)[:] = advector_0
                self.mpdata[field].advance(1)
        else:
            field = "qvap"
            self.mpdata[field].advector.get_component(0)[:] = advector_0
            self.mpdata[field].advance(1)

        thermo = self.set_thermo(time, thermo)
        return thermo

    def set_advectees(self, thermo):
        """
        Set the values of variables that are advected 1-D KiD motion computations.

        This method sets the values of advected quantities in KiD dynamics ("advectees")
        to match the values given by the thermo data structure.

        It is neccessary to call this function after changes to thermo if the data in thermo
        was created from a copy not a reference to the advectees.

        Args:
            thermo (Thermodynamics): Object representing the thermodynamic state.
        """
        self.mpdata["qvap"].advectee.get()[:] = thermo.massmix_ratios["qvap"]
        self.mpdata["qcond"].advectee.get()[:] = thermo.massmix_ratios["qcond"]
        self.mpdata["qice"].advectee.get()[:] = thermo.massmix_ratios["qice"]
        self.mpdata["qrain"].advectee.get()[:] = thermo.massmix_ratios["qrain"]
        self.mpdata["qsnow"].advectee.get()[:] = thermo.massmix_ratios["qsnow"]
        self.mpdata["qgrau"].advectee.get()[:] = thermo.massmix_ratios["qgrau"]
