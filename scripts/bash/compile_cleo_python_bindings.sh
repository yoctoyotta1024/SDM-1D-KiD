#!/bin/bash

### script to compile Cleo's python bindings
### either for Levante HPC or a generic "vanilla" computer

### ------------------ input parameters ---------------- ###
### ------- You may want to edit these lines to -------- ###
### --------- specify your build configuration --------- ###
### ---------------------------------------------------- ###

### ------------- generic input parameters ------------- ###
### _NOTE_: best to use absolute paths here
machine=${1:-"levante"}                                    # machine, either "levante" or "vanilla"
path2sdm1dkid=${2:-"/home/m/m300950/SDM-1D-KiD"}           # path to SDM-1D-KiD repository
path2build=${3:-"/home/m/m300950/SDM-1D-KiD/build"}        # path to build directory (will be created)
python=${4:-"/home/m/m300950/SDM-1D-KiD/.venv/bin/python"} # path to sdm1dkid python interpreter
yacyaxtroot=${5:-"/work/mh0731/m300950/yacyaxt/gcc"}       # path to [...]/yac and [...]/yaxt directories for Cleo
### ---------------------------------------------------- ###

### --------- extra build parameters for Cleo ---------- ###
### _NOTE_: most build flags for Cleo are already set in src/extern/cleo/CMakeLists.txt
cleo_build_flags="-DCLEO_PYTHON=${python}"
make_clean=false                                           # == "true" or otherwise false
executable="cleo_python_bindings"
path2source=${path2sdm1dkid}/src/extern/cleo
cleo_yac_module_path="${path2build}/_deps/cleo-src/libs/configuration/cmake"
cmake_nproc=16
### ---------------------------------------------------- ###

### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ------------------ check arguments ----------------- ###
if [[ "${machine}" != "levante" && "${machine}" != "vanilla" ]]; then
    echo "Error: ${machine} machine not recognised. Please choose either 'levante' or 'vanilla'"
    exit 1
fi

if [ "${path2source}" == "${path2build}" ]
then
    echo "Error: build directory cannot match source directory"
    exit 1
fi
### ---------------------------------------------------- ###

### ----------------- set Cleo flags ------------------- ###
CLEO_BUILD_FLAGS="${cleo_build_flags}"
CLEO_YAC_FLAGS="-DCLEO_YAC_MODULE_PATH=${cleo_yac_module_path} \
  -DCLEO_YAXT_ROOT=${yacyaxtroot}/yaxt \
  -DCLEO_YAC_ROOT=${yacyaxtroot}/yac"
CLEO_KOKKOS_BASIC_FLAGS="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON"
CLEO_KOKKOS_HOST_FLAGS="-DKokkos_ENABLE_OPENMP=ON"
CLEO_KOKKOS_DEVICE_FLAGS=""
### ---------------------------------------------------- ###

### ---- set relevant packages and compiler settings --- ###
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/helpers" &>/dev/null && pwd)/${machine}.sh"
load_build_packages
### ---------------------------------------------------- ###

### ----------------- print build bindings ------------- ###
echo "### --------------- Build Inputs -------------- ###"
echo "machine=${machine}"
echo "python=${python}"
echo "source directory=${path2source}"
echo "build directory=${path2build}"
echo "CMAKE_CXX_COMPILER=${cxx_compiler}"
echo "CMAKE_C_COMPILER=${c_compiler}"
echo "CMAKE_CXX_FLAGS="${cxx_flags}""

echo "CLEO_BUILD_FLAGS=${CLEO_BUILD_FLAGS}"
echo "CLEO_YAC_FLAGS=${CLEO_YAC_FLAGS}"
echo "CLEO_KOKKOS_FLAGS=${CLEO_KOKKOS_BASIC_FLAGS}\
  ${CLEO_KOKKOS_HOST_FLAGS} ${CLEO_KOKKOS_DEVICE_FLAGS}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### ------------------- build bindings ----------------- ###
cmake -DCMAKE_CXX_COMPILER=${cxx_compiler} \
    -DCMAKE_C_COMPILER=${c_compiler} \
    -DCMAKE_CXX_FLAGS="${cxx_flags}" \
    -DCMAKE_BUILD_TYPE=Release \
    -S ${path2source} -B ${path2build} \
    ${CLEO_KOKKOS_BASIC_FLAGS} ${CLEO_KOKKOS_HOST_FLAGS} ${CLEO_KOKKOS_DEVICE_FLAGS} \
    ${CLEO_BUILD_FLAGS} ${CLEO_YAC_FLAGS}
### ---------------------------------------------------- ###

### ------------------ compile bindings ---------------- ###
cmd="cd ${path2build}"
echo ${cmd}
eval ${cmd}

if [ "${make_clean}" == "true" ]
then
  cmd="make clean"
  echo ${cmd}
  eval ${cmd}
fi

cmd="make -j ${cmake_nproc} ${executable}"
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###
