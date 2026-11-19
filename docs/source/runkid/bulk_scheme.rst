Running the 1-D KiD with PyMPDATA Bulk Scheme
===============================================

The 1-D KiD can be run with a simple bulk scheme condensation saturation
adjustment from the PyMPDATA libray with PyMPDATA as the thermodynamic driver.

Running The PyMPDATA Bulk Scheme
----------------------------------

The script ``scripts/run_bulk_scheme.py`` can be used to run the entire workflow
for an ensemble of runs of the bulk scheme 1-D KiD. For each member of the
ensemble, it consists of 4 stages:

1. :ref:`generating inputfiles<do_inputfiles_bulk>`
2. :ref:`running, i.e. time-stepping, the bulk scheme 1-D KiD model<do_run_model_bulk>`
3. :ref:`post-processing the bulk scheme 1-D KiD output<do_post_process_bulk>`
4. :ref:`plotting "quicklooks" of the post-processed results<do_plot_results_bulk>`

Each stage can be turned on (the default is off) by calling ``run_bulk_scheme.py``
with the associated flag for that stage, e.g. to run all the stages for
an ensemble called "your_ensemble" consisting of 3 runs you would do:

.. code-block:: console

  $ uv run python ./scripts/run_bulk_scheme.py /path/to/your/build "your_ensemble" 3 \
      --do_inputfiles --do_run_model --do_post_process --do_plot_results

To make things potentially easier, you can also choose to call this script through
the bash script ``run_sdm_1d_kid.sh``:

.. code-block:: console

  $ ./scripts/bash/run_sdm_1d_kid.sh [machine] [model] [path2sdm1dkid] [path2build] [python] [run_flags]


E.g. on Levante:

.. code-block:: console

  $ ./scripts/bash/run_sdm_1d_kid.sh levante bulk \
       "--do_inputfiles --do_run_model --do_post_process --do_plot_results" \
       ./ ./build $(uv run which python)

The following sections describe each stage in a bit more detail.


.. _do_inputfiles_bulk:

Generating Inputfiles
~~~~~~~~~~~~~~~~~~~~~~~

To generate the input files  for each run
in your ensemble, call:

.. code-block:: console

  $  uv run python ./scripts/run_bulk_scheme.py /path/to/your/build "your_ensemble" 3 --do_inputfiles


This will actually just make some directories so that you can run the model.
The actualy configuration of the model for each run is set in the
``--do_run_model`` stage, through a python dictionary.


.. _do_run_model_bulk:

Running, i.e. Time-stepping, the Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running the 1-D KiD for the bulk scheme consits of two steps:
* generate the configuration (see ``generate_configuration`` function)
* run the model (see ``run_model`` function)

The function ``generate_configuration`` in ``run_bulk_scheme.py`` contains all
the settings you may want to configure yourself. The configuration for each run
is then passed through to the ``run_model`` function and so no configuration
file is ever written  (unlike for Cleo).

By default, for an ensemble called "your_ensemble", this script will generate
the output (a Zarr direcotry and quickplot of the moisture evolution) for each
run in ``/path/to/your/build/bulk_scheme/your_ensemble/bin/``.


.. _do_post_process_bulk:

Post-Processing the Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~
###WIP

.. _do_plot_results_bulk:

Plotting the "Quicklooks"
~~~~~~~~~~~~~~~~~~~~~~~~~~~
###WIP
