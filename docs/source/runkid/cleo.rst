Running the 1-D KiD with Cleo
==============================

The 1-D KiD with Cleo uses the PyMPDATA library as the thermodynamic driver. To
run Cleo's C++ code alongside PyMPDATA, we first make the Cleo's python bindings
(python module called ``cleo_python_bindings``) using pybind11.

Creating The Python Bindings
----------------------------

Cleo's libraries are fetched and built via the ``CMakeLists.txt`` in the
directory ``src/extern/cleo/``.

To build Cleo's python bindings you can simply execute the bash script
``compile_cleo_python_bindings.sh`` in ``scripts/bash/`` as follows:

.. code-block:: console

  $ ./scripts/bash/compile_cleo_python_bindings.sh [machine] [path2sdm1dkid] [path2build] [python] [yacyaxtroot]

* ``machine`` is the computer architecture you are using. Currently supported is
  either ``levante`` or ``vanilla``. The "vanilla" architecture is supposed to run
  on any generic computer, e.g. your MacBook or Ubuntu Linux machine.

* ``path2sdm1dkid`` is the (absolute) path to your SDM-1D-KiD repository,
  e.g. ``${HOME}/SDM-1D-KiD``.

* ``path2build`` is the (absolute) path for where you want to build the python
  bindings, e.g. ``/path/to/SDM-1D-KiD/build``.

* ``python`` is the (absolute) path to your python interpreter. *hint:* with
  ``uv`` you can set this to ``$(uv run which python)``.

* ``yacyaxtroot`` is the (absolute) path to the parent directory of your
  installations of YAC and YAXT, i.e. ``/path/to/yacyaxt/`` which contains
  a directory called ``yac`` and a directory called ``yaxt`` containing each
  installation respectively. (See :doc:`../usage/requirements` for more details).

E.g. on Levante:

.. code-block:: console

  $ ./scripts/bash/compile_cleo_python_bindings.sh levante ./ ./build \
      $(uv run which python) /work/mh0731/m300950/yacyaxt/gcc


After making the bindings, the ``cleo_python_bindings`` python module can be used
just like any other ordinary python module. By default, the bindings are built
in ``/path/to/your/build/_deps/cleo-build/cleo_python_bindings/``.


Running Cleo
-------------

The script ``scripts/run_cleo.py`` can be used to run the entire workflow for
an ensemble of runs of the Cleo 1-D KiD. It assumes you have already
compiled Cleo's python bindings and for each member of the ensemble,
it consists of 4 stages:

1. :ref:`generating inputfiles<do_inputfiles_cleo>`
2. :ref:`running, i.e. time-stepping, the Cleo 1-D KiD model<do_run_model_cleo>`
3. :ref:`post-processing the Cleo 1-D KiD output<do_postprocess_cleo>`
4. :ref:`plotting "quicklooks" of the post-processed results<do_plot_results_cleo>`

Each stage can be turned on (the default is off) by calling ``run_cleo.py``
with the associated flag for that stage, e.g. to run all the stages for
an ensemble called "your_ensemble" consisting of 3 runs you would do:

.. code-block:: console

  $ uv run python ./scripts/run_cleo.py /path/to/your/build "your_ensemble" 3 \
      --do_inputfiles --do_run_model --do_postprocess --do_plot_results

To make things potentially easier, you can also choose to call this script through
the bash script ``run_sdm_1d_kid.sh``:

.. code-block:: console

  $ ./scripts/bash/run_sdm_1d_kid.sh [machine] [model] [run_flags] [path2sdm1dkid] [path2build] [python]


E.g. on Levante:

.. code-block:: console

  $ ./scripts/bash/run_sdm_1d_kid.sh levante cleo \
       "--do_inputfiles --do_run_model --do_postprocess --do_plot_results" \
       ./ ./build $(uv run which python)

The following sections describe each stage in a bit more detail.


.. _do_inputfiles_cleo:

Generating Inputfiles
~~~~~~~~~~~~~~~~~~~~~~~

To generate the input files for Cleo for each run in your ensemble, call:

.. code-block:: console

  $  uv run python ./scripts/run_cleo.py /path/to/your/build "your_ensemble" 3 --do_inputfiles


This will generate three files for each run:

* the configuration file,
* the grid file,
* the superdroplet initial conditions file.

The function ``inputfiles`` in ``run_cleo.py`` contains all the settings
you may want to configure yourself. You can also have a look at the source
it calls, ``src/sdm_1d_kid/cleo/generate_inputfiles.py``, for more details.

By default, for an ensemble called "your_ensemble", this script will generate the
configuration file for your setup in ``/path/to/your/build/cleo/your_ensemble/tmp/``,
and the grid and superdroplet initial conditions files in
``/path/to/your/build/cleo/your_ensemble/share/``. It will also plot some figures
about the grid and superdroplet initial conditions in
``/path/to/your/build/cleo/your_ensemble/bin/``.


.. _do_run_model_cleo:

Running, i.e. Time-stepping, the Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the Cleo 1-D KiD with the ``--do_run_model`` flag, either via ``run_sdm_1d_kid.sh``
or directly via ``run_cleo.py``. This will output a ``.txt`` setup file
and a ``.zarr`` output data directory, as well as a quickplot of the moisture
evolution for each run. By default, these files will be written in
``/path/to/your/build/cleo/your_ensemble/bin/`` (as per the settings in your
YAML config files, generated by the ``inputfiles`` function in ``run_cleo.py``).


.. _do_postprocess_cleo:

Post-Processing the Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the ``--do_postprocess`` stage to combine the ``.zarr`` output of all the runs in
your ensemble into a single post-processed dataset (see
``src/sdm_1d_kid/cleo/postprocess.py`` for details). By default, this post-processed
dataset is written in ``/path/to/your/build/cleo/your_ensemble/postproc/``
(see ``run_cleo.py``) as ``your_ensemble.zarr``.

The dataset has dimensions ``(ensemble, height, time)``, plus
``(ensemble, height, time_sds)`` for the superdroplet diagnostics and
``(ensemble, height, time_sds, radius_bin)`` or ``(ensemble, height, time_sds, radius_bin_dry)``
for the superdroplet (dry) size distributions. The radius bins are geometrically spaced,
i.e. equally spaced in logarithmic space. The superdroplet diagnostics are calculated
every ``sddiags_every``-th output time, i.e. at the times of the ``time_sds`` coordinate.
The parameters for these diagnostics, including the radius bins
(``radius_min``, ``radius_max`` and ``nbins``, and likewise for dry radii),
can be edited in the ``main`` function of ``postprocess.py``.


.. _normalising_distribs:

Normalising the Droplet Size Distributions
###########################################

The size distributions (``numconc_distrib``, ``mass_distrib``, ``reflectivity_distrib``,
``msol_mass_distrib``, etc.) are *not* normalised by the bin width. To normalise them,
divide by the bin-width coordinates, e.g.

.. code-block:: python

  import xarray as xr

  ds = xr.open_zarr("/path/to/your/build/cleo/your_ensemble/postproc/your_ensemble.zarr")

  dndlnr = ds.numconc_distrib / ds.radius_bin_dlnr  # number concentration per unit ln(radius), dN/dlnr [cm^-3]
  dmdlnr = ds.mass_distrib / ds.radius_bin_dlnr  # droplet mass concentration per unit ln(radius), dM/dlnr [g m^-3]
  dzdlnr = ds.reflectivity_distrib / ds.radius_bin_dlnr  # reflectivity per unit ln(radius), dZ/dlnr [mm^6 m^-3]

or per log10(r) instead of ln(r):

.. code-block:: python

  import numpy as np
  [...]
  dlogr = ds.radius_bin_dlnr / np.log(10)
  dndlog10r = ds.numconc_distrib / dlogr  # number concentration per unit log10(radius), dN/dlogr [cm^-3]

or per unit radius (linear instead of ln(r)):

.. code-block:: python

  [...]
  dr = (ds.radius_bin_upper - ds.radius_bin_lower) * 1e6  # [um]
  dndr = ds.numconc_distrib / dr  # number concentration per unit radius, dN/dr [cm^-3 um^-1]

The radius bins are logarithmically spaced, so ``radius_bin_dlnr`` (and ``dlogr``)
are the same for every bin. Normalising by linear radius (``dN/dr``) is different because
the linear bin widths grow with r, so it will change the shape of the distribution a lot.

Sanity check: integrating back the normalised distributions should recover the totals
(to float32 precision), e.g.
``(dndlnr * ds.radius_bin_dlnr).sum("radius_bin")`` ≈ ``ds.numconc.sel(time=ds.time_sds)``.

.. _do_plot_results_cleo:

Plotting the "Quicklooks"
~~~~~~~~~~~~~~~~~~~~~~~~~~~
###WIP


Comparison to the PySDM 1-D KiD
--------------------------------
You may want to have a look at the similarities and differences between the PySDM
and Cleo KiD setups here: :doc:`cleo_pysdm_comparison`.
