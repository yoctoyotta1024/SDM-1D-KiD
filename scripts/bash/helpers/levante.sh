#!/bin/bash

# Necessary Levante packages for compiling Cleo's python bindings
levante_gcc="gcc/11.2.0-gcc-11.2.0" # bcn7mbu # module load
levante_gcc_openmpi="openmpi/4.1.2-gcc-11.2.0" # module load
levante_gcc_cmake="cmake@3.26.3%gcc@=11.2.0/fuvwuhz" # spack load
levante_gcc_netcdf_yac="netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0" # module load
levante_gcc_openblas_yac="openblas@0.3.18%gcc@=11.2.0" # spack load

load_build_packages() {
  module load ${levante_gcc} ${levante_gcc_openmpi} ${levante_gcc_netcdf_yac} # for Cleo and YAC
  module list
  spack load ${levante_gcc_cmake} ${levante_gcc_openblas_yac} # for Cleo and YAC
  spack find --loaded -p
  export cxx_compiler="$(command -v mpic++)"
  export c_compiler="$(command -v mpicc)"
  export cxx_flags="-Werror -Wall -Wextra -pedantic -Wno-unused-parameter -O3 -mfma"
}

load_runtime_packages() {
  local path2sdm1dkid="$1"
  local path2build="$2"

  module load ${levante_gcc} ${levante_gcc_openmpi} ${levante_gcc_netcdf_yac}
  module list
  echo "export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64:$LD_LIBRARY_PATH"
  export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64:$LD_LIBRARY_PATH

  # (optional) prepend to python path to make import searches faster
  site_packages="$(uv run python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")"
  path2cleo_python_bindings="${path2build}/_deps/cleo-build/cleo_python_bindings"
  echo "export PYTHONPATH=${site_packages}:${path2cleo_python_bindings}:${path2sdm1dkid}:${PYTHONPATH}"
  export PYTHONPATH=${site_packages}:${path2cleo_python_bindings}:${path2sdm1dkid}:${PYTHONPATH}
}
