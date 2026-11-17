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
