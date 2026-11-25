Running the 1-D KiD with PyMPDATA Bulk Scheme
===============================================

The 1-D KiD can be run with a simple bulk scheme condensation saturation
adjustment from the PyMPDATA library with PyMPDATA as the thermodynamic driver.

Running The PyMPDATA Bulk Scheme
----------------------------------

The script ``scripts/run_bulk_scheme.py`` can be used to run the entire workflow
for an ensemble of runs of the bulk scheme 1-D KiD. For each member of the
ensemble, it consists of 4 stages:

1. :ref:`generating inputfiles<do_inputfiles_bulk>`
2. :ref:`running, i.e. time-stepping, the bulk scheme 1-D KiD model<do_run_model_bulk>`
3. :ref:`post-processing the bulk scheme 1-D KiD output<do_postprocess_bulk>`
4. :ref:`plotting "quicklooks" of the post-processed results<do_plot_results_bulk>`

Each stage can be turned on (the default is off) by calling ``run_bulk_scheme.py``
with the associated flag for that stage, e.g. to run all the stages for
an ensemble called "your_ensemble" consisting of 3 runs you would do:

.. code-block:: console

  $ uv run python ./scripts/run_bulk_scheme.py /path/to/your/build "your_ensemble" 3 \
      --do_inputfiles --do_run_model --do_postprocess --do_plot_results

To make things potentially easier, you can also choose to call this script through
the bash script ``run_sdm_1d_kid.sh``:

.. code-block:: console

  $ ./scripts/bash/run_sdm_1d_kid.sh [machine] [model] [run_flags] [path2sdm1dkid] [path2build] [python]


E.g. on Levante:

.. code-block:: console

  $ ./scripts/bash/run_sdm_1d_kid.sh levante bulk \
       "--do_inputfiles --do_run_model --do_postprocess --do_plot_results" \
       ./ ./build $(uv run which python)

The following sections describe each stage in a bit more detail.


.. _do_inputfiles_bulk:

Generating Inputfiles
~~~~~~~~~~~~~~~~~~~~~~~

To generate the input files for each run
in your ensemble, call:

.. code-block:: console

  $  uv run python ./scripts/run_bulk_scheme.py /path/to/your/build "your_ensemble" 3 --do_inputfiles


This will actually just make some directories so that you can run the model.
The actually configuration of the model for each run is set in the
``--do_run_model`` stage, through a python dictionary.


.. _do_run_model_bulk:

Running, i.e. Time-stepping, the Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running the 1-D KiD for the bulk scheme consists of two steps:

* generate the configuration (see ``generate_configuration`` function)
* run the model (see ``run_model`` function)

The function ``generate_configuration`` in ``run_bulk_scheme.py`` contains all
the settings you may want to configure yourself. The configuration for each run
is then passed through to the ``run_model`` function and so no configuration
file is ever written (unlike for Cleo).

By default, for an ensemble called "your_ensemble", this script will generate
the output (a Zarr directory and quickplot of the moisture evolution) for each
run in ``/path/to/your/build/bulk_scheme/your_ensemble/bin/``.

*NOTE:* The bulk scheme is deterministic, so all members of an ensemble are
identical. It therefore doesn't make sense to run more than one member, however,
you will still want to post-process the raw output as if it was an
ensemble (see :ref:`do_postprocess_bulk`).


.. _do_postprocess_bulk:

Post-Processing the Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the ``--do_postprocess`` stage to combine the ``.zarr`` output of all the runs in
your ensemble into a single post-processed dataset (see
``src/sdm_1d_kid/pympdata_driver/bulk_scheme/postprocess.py`` for details). By default,
this dataset is written in ``/path/to/your/build/bulk_scheme/your_ensemble/postproc/``
(see ``run_bulk_scheme.py``) as ``your_ensemble.zarr``.

The post-processed dataset has dimensions ``(ensemble, height, time)`` and uses the same
variable names, units and thermodynamic diagnostics (e.g. ``relh``, ``theta``,
``rho_dry``) as the post-processed Cleo dataset, so that the two can be compared
directly.

*NOTE:* The coordinates and variables in the bulk scheme dataset are an exact
subset of the post-processed ensemble dataset from Cleo, except it may also
contain the variables ``qrain`` and ``rhod``, which the Cleo dataset does not have.
Also the bulk scheme's ``lwc`` and ``lwp`` contain no solute (unlike Cleo's).

.. _interp_datasets:

Interpolating Between Datasets
##################################

The bulk scheme output may not be on the same timestep as Cleo's outputs,
so to compare the two datasets you may need to match their times, e.g.

.. code-block:: python

  import numpy as np
  import xarray as xr

  ds_bulk = xr.open_zarr("/path/to/your/build/bulk_scheme/your_ensemble/postproc/your_ensemble.zarr")
  ds_cleo = xr.open_zarr("/path/to/your/build/cleo/your_ensemble/postproc/your_ensemble.zarr")

  if np.all(np.isin(ds_cleo.time.values, ds_bulk.time.values)):
      ds_bulk_matched = ds_bulk.sel(time=ds_cleo.time)  # exactly matching times
  else:
      ds_bulk_matched = ds_bulk.interp(time=ds_cleo.time)  # linear interpolation

If Cleo's output times are a subset of the bulk scheme's (i.e. ``OBSTSTEP`` is a
multiple of the bulk scheme's timestep), ``sel`` selects the actual model values;
otherwise ``interp`` interpolates linearly in time (NaN outside the bulk run's time
range). To compare with Cleo's superdroplet diagnostics, select their times instead,
e.g. ``ds_bulk_matched.sel(time=ds_cleo.time_sds)``.


.. _do_plot_results_bulk:

Plotting the "Quicklooks"
~~~~~~~~~~~~~~~~~~~~~~~~~~~
###WIP
