Running the 1-D KiD with Cleo
==============================

The 1-D KiD with Cleo uses the PYMPDATA libray as the thermodynamic driver. To
run Cleo's C++ code alongside PYMPDATA, we first make the Cleo's python bindings
(Python module called ``cleo_python_bindings``) using pybind11.

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

* ``yacyaxtroot`` is the (absolute) path to the parent directory of your
  installations of YAC and YAXT, i.e. ``/path/to/yacyaxt/`` which contains
  a directory called ``yac`` and a directory called ``yaxt`` containing each
  installation respectively. (See :doc:`../usage/requirements` for more details).

After making the bindings, the ``cleo_python_bindings`` Python module can used just
like any other ordinary Python module. By default, the bindings are built in
``/path/to/your/build/_deps/cleo-build/cleo_python_bindings/``.


Running Cleo
-------------

The script ``scripts/run_cleo.py`` can be used to run the entire workflow for
the Cleo 1-D KiD. It assumes you have already compiled Cleo's python bindings
and it consists of 4 stages:

1. :ref:`generating inputfiles<do_inputfiles>`
2. :ref:`running, i.e. time-stepping, the Cleo 1-D KiD model<do_run_model>`
3. :ref:`post-processing the Cleo 1-D KiD output<do_post_process>`
4. :ref:`plotting "quicklooks" of the post-processed results<do_plot_results>`

Each stage can be turned on (the default is off) by calling ``run_cleo.py``
with the associated flag for that stage, e.g. to run all the stages
you would do:

.. code-block:: console

  $  uv run python ./scripts/run_cleo.py /path/to/your/build \
      --do_inputfiles --do_run_model --do_post_process --do_plot_results

The following sections describe each stage in a bit more detail.


.. _do_inputfiles:

Generating Inputfiles
~~~~~~~~~~~~~~~~~~~~~~~

To generate the input files to run Cleo call:

.. code-block:: console

  $  uv run python ./scripts/run_cleo.py /path/to/your/build --do_inputfiles


This will generate three files:
* the configuration file,
* the grid file,
* the superdroplet initial conditions file.

The function ``inputfiles`` in ``run_cleo.py`` contains all the settings
you may want to configure yourself. You can also have a look at the source
it calls, ``src/sdm_1d_kid/cleo/generate_inputfiles.py``, for more details.

By default, this script will generate the configuration file for your
setup in ``/path/to/your/build/tmp/``, and the grid and superdroplet
initial conditions files in ``/path/to/your/build/share/``. It will also plot
some figures about the grid and superdroplet initial conditions in
``/path/to/your/build/bin/``.


.. _do_run_model:

Running, i.e. Time-stepping, the Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. _do_post_process:

Post-Processing the Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. _do_plot_results:

Plotting the "Quicklooks"
~~~~~~~~~~~~~~~~~~~~~~~~~~~
