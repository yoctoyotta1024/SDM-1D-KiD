"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- SDM-1D-KiD -----
File: output_thermodynamics.py
Project: thermo
Created Date: Thursday 19th November 2026
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
class for storing thermodynamics output during model run
"""

import numpy as np

from sdm_1d_kid.pympdata_driver.thermo.thermodynamics import Thermodynamics


class OutputVariable:
    """Class to output a variable with some of its metadata.

    Attributes:
        name (string):
          Name of variable.
        units (string):
          Units of variable.
        values (np.ndarray):
          Values of variable; must be able to be appended to.
    """

    def __init__(self, name, units, shape):
        """Initialise OutputVariable instance

        Parameters:
          name (string):
            Name of variable.
          units (string):
            Units of variable.
          shape (array of integers)
            Shape of final array. Note that values are added by write function
            as slices along first dimension (e.g. first dimension is time)
        """

        self.name = name
        self.units = units
        self.values = np.zeros(shape)
        self._i = 0

    def write(self, val):
        """Copy values 'val' to the array in slice along its first dimension.

        Parameters:
          val (any): Value(s) to append.
        """
        assert np.shape(self.values)[0] > self._i, "no space left in array"
        assert np.shape(val) == np.shape(self.values)[1:], (
            "values to add must be complete slice along 0th dimension"
        )
        self.values[self._i] = val
        self._i += 1

    def set(self, val):
        """Copy values 'val' to the entire array.

        Parameters:
          val (any): Value(s) of array.
        """
        self.values = val
        self._i = np.nan  # destroy ability of write() function

    def finalize(self):
        """Convert values to a NumPy array."""
        print(self.name + " shape: " + str(self.values.shape))


class OutputThermodynamics:
    """Class is method and store for thermodynamic variables output during model timestep.

    Thermodynamic variables include pressure, temperature, moist air density and the specific
    content (mass mixing ratio) of vapour and condensates.

    Attributes:
        time (OutputVariable):
          time (s).
        temp (OutputVariable):
          Temperature (K).
        rhod (OutputVariable):
          Density of dry air (kg/m3).
        press (OutputVariable):
          Pressure (Pa).
        qvap (OutputVariable):
          Specific water vapor content (kg/kg).
        qcond (OutputVariable):
          Specific cloud water content (kg/kg).
        qice (OutputVariable):
          Specific cloud ice content (kg/kg).
        qrain (OutputVariable):
          Specific rain content (kg/kg).
        qsnow (OutputVariable):
          Specific snow content kg/kg).
        qgrau (OutputVariable):
          Specific graupel content (kg/kg).
    """

    def __init__(self, shape, zhalf=None, xhalf=None, yhalf=None, zarrbasedir=None):
        """Initialize an OutputThermodynamics object.

        Args:
            shape (tuple): Shape of the output variables, time is first dimension.
            zhalf (np.ndarray, optional): Half-level z-coordinates (m). Defaults to None.
            xhalf (np.ndarray, optional): Half-level x-coordinates (m). Defaults to None.
            yhalf (np.ndarray, optional): Half-level y-coordinates (m). Defaults to None.
        """

        self.time = OutputVariable("time", "s", [shape[0]])
        self.temp = OutputVariable("temp", "K", shape)
        self.rhod = OutputVariable("rhod", "kg m-3", shape)
        self.press = OutputVariable("press", "Pa", shape)
        self.qvap = OutputVariable("qvap", "kg/kg", shape)
        self.qcond = OutputVariable("qcond", "kg/kg", shape)
        self.qice = OutputVariable("qice", "kg/kg", shape)
        self.qrain = OutputVariable("qrain", "kg/kg", shape)
        self.qsnow = OutputVariable("qsnow", "kg/kg", shape)
        self.qgrau = OutputVariable("qgrau", "kg/kg", shape)

        if zhalf is not None:
            self.zhalf = OutputVariable("zhalf", "m", [len(zhalf)])
            self.zhalf.set(zhalf)
        else:
            self.zhalf = None

        if xhalf is not None:
            self.xhalf = OutputVariable("xhalf", "m", [len(xhalf)])
            self.xhalf.set(xhalf)
        else:
            self.xhalf = None

        if yhalf is not None:
            self.yhalf = OutputVariable("yhalf", "m", [len(yhalf)])
            self.yhalf.set(yhalf)
        else:
            self.yhalf = None

        self.zarrbasedir = zarrbasedir

    def _time_varying_2d_variables(self):
        """Returns list of output variables with dimensions (time, z)."""
        return [
            self.temp,
            self.rhod,
            self.press,
            self.qvap,
            self.qcond,
            self.qice,
            self.qrain,
            self.qsnow,
            self.qgrau,
        ]

    @staticmethod
    def _global_attrs():
        """Metadata For Xarray Dataset."""
        from datetime import UTC, datetime

        return {
            "title": "1-D KiD with PyMPDATA bulk scheme",
            "created": datetime.now(UTC).isoformat(),
        }

    def to_xarray(self):
        """Return the output as an xarray Dataset with dimensions (zfull, time).

        Only the time slices which have actually been written are included. Variables are
        defined at the centres of gridboxes (full levels, "z"); the gridbox boundaries
        (half levels) are included as the coordinate "zhalf" if given.
        """
        import xarray as xr

        nt = self.time._i  # number of output times actually written
        long_names = {
            "temp": "temperature",
            "rhod": "dry air density",
            "press": "pressure",
            "qvap": "water vapour mass mixing ratio",
            "qcond": "cloud condensate mass mixing ratio",
            "qice": "cloud ice mass mixing ratio",
            "qrain": "rain mass mixing ratio",
            "qsnow": "snow mass mixing ratio",
            "qgrau": "graupel mass mixing ratio",
        }

        coords = {
            "time": ("time", self.time.values[:nt], {"units": "s", "long_name": "time"})
        }
        if self.zhalf is not None:
            zhalf = np.asarray(self.zhalf.values)
            coords["zhalf"] = (
                "zhalf",
                zhalf,
                {
                    "units": "m",
                    "long_name": "height of gridbox boundaries (half levels)",
                },
            )
            coords["zfull"] = (
                "zfull",
                0.5 * (zhalf[1:] + zhalf[:-1]),
                {"units": "m", "long_name": "height of gridbox centres (full levels)"},
            )
        for name, var in (("xhalf", self.xhalf), ("yhalf", self.yhalf)):
            if var is not None:
                coords[name] = (name, np.asarray(var.values), {"units": var.units})

        data_vars = {
            var.name: (
                ("zfull", "time"),
                var.values[:nt].T,
                {"units": var.units, "long_name": long_names[var.name]},
            )
            for var in self._time_varying_2d_variables()
        }

        return xr.Dataset(
            data_vars=data_vars, coords=coords, attrs=self._global_attrs()
        )

    def output_thermodynamics(self, time: float, thermo: Thermodynamics):
        """output thermodynamics from thermo to each variable in thermodynamics output.

        This method writes time and thermodynamic variables from thermo such as temperature, density,
        pressure, and specific mass mixing ratios to the respective output variables.

        Parameters:
            time (float):
              The time at which the thermodynamic variables are output (s).
            thermo (Thermodynamics):
              An instance of the Thermodynamics class containing the thermodynamic variables to output.

        Returns:
            None
        """
        self.time.write(time)
        self.temp.write(thermo.temp)
        self.rhod.write(thermo.rhod)
        self.press.write(thermo.press)
        self.qvap.write(thermo.massmix_ratios["qvap"])
        self.qcond.write(thermo.massmix_ratios["qcond"])
        self.qice.write(thermo.massmix_ratios["qice"])
        self.qrain.write(thermo.massmix_ratios["qrain"])
        self.qsnow.write(thermo.massmix_ratios["qsnow"])
        self.qgrau.write(thermo.massmix_ratios["qgrau"])

    def __call__(self, time: float, thermo: Thermodynamics):
        """Invoke the object as a function to call the `output_thermodynamics` method.

        Parameters:
            time (float):
              The time at which the thermodynamic variables are output (s).
            thermo (Thermodynamics):
              An instance of the Thermodynamics class containing the thermodynamic variables to output.

        Returns:
            None
        """
        self.output_thermodynamics(time, thermo)

    def __getitem__(self, key):
        """
        Get the value of an attribute by its name using bracket notation.

        Parameters:
            key (str): The name of the attribute to access.

        Returns:
            object: The value of the attribute if it exists, otherwise None.
        """
        return getattr(self, key, None)

    def finalize(self):
        """Finalize the thermodynamics output.

        This method finalizes the output of thermodynamic variables by calling the `finalize` method
        of each output variable, e.g. to ensure that all variables are properly formatted
        for further use or analysis.

        Returns:
            None
        """
        for var in [self.time, *self._time_varying_2d_variables()]:
            var.finalize()
        for var in (self.zhalf, self.xhalf, self.yhalf):
            if var is not None:
                var.finalize()

        if self.zarrbasedir is not None:
            ds = self.to_xarray()
            ds.to_zarr(self.zarrbasedir, mode="w", zarr_format=2)
            print(f"output thermodynamics written to: {self.zarrbasedir}")
