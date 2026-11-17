.. _getstart:

Getting Started
===============

Clone the GitHub repository:

.. code-block:: console

  $ git clone https://github.com/yoctoyotta1024/SDM-1D-KiD.git

We provide a ``pyproject.toml`` file for creating the python version and all the pip-installable
packages you need to use this repository (e.g to run the scripts, to build the documentation,
or to run pre-commit etc.). We provide a ``environment.yaml`` file for all non-pip-installable
packages available from conda-forge.

Everything will be much easier for you if you install the python version and package depenencies as
per our ``pyproject.toml`` file. You can acheive this quickly and simply using
`uv <https://docs.astral.sh/uv/>`_. *hint*: optionally, you can install ``uv`` within a
Conda/Micromamba environment, which also installs our other non-pip-installable packages (e.g. doxygen)
using our ``environment.yaml`` file, e.g.

.. code-block:: console

  $ mamba create --file=environment.yaml
  $ mamba activate sdm1dkid

Using ``uv`` to setup python, simply do:

.. code-block:: console

  $ uv sync

Sometimes it is useful to know where uv has installed your python, in which case you can run:

.. code-block:: console

  $ uv python find

Finally, install our pre-commit hooks:

.. code-block:: console

  $ uv run pre-commit install

That's it, you're done!

We kindly ask that you also :ref:`contact us <contact>` or
`open a new issue <https://github.com/yoctoyotta1024/SDM-1D-KiD/issues/new>`_
on our GitHub repository if you discover this environment is missing any packages.

Using PyBind11 on ARM64 (Macbooks with M1 Processors) instead of x86_64 Architectures
-------------------------------------------------------------------------------------
If you are using Python on an ARM64 architecture (e.g. you have a Mac with an M1 silicon chip) then
you need to ensure the python libraries have been installed for ARM64 and *not* x86_64
architectures. You can check this by entering either ``file [your_python_library].dylib`` or
``file [your_python_library].so``, which will tell you which architecture your python is installed
for.

If you don't know where your python libraries are installed, have a look at ``which python`` to find
the path to your python interpreter, e.g. something like ``/path/to/python/bin/python``. Then
look in ``/path/to/python/lib/`` for something like ``libpython[X.Y].dylib`` or ``libpython[X.Y].so``
where X.Y is your python version (which you can get from ``python --version``).

If your python libraries are not ARM64, then you will have problems... One solution is to create a
conda enviroment which specifically uses conda-forge's osx-arm64 directory to find packages, e.g. via

.. code-block:: console

  $ CONDA_SUBDIR=osx-arm64 conda create -n [my_arm64_env] numpy -c conda-forge
  $ conda activate [name_my_arm64_env]
  $ conda config --env --set subdir osx-arm64

where ``[name_my_arm64_env]`` is the name of your new enviroment.
