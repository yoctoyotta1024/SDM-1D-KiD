#!/bin/bash

### script to launch python scripts for ensembles of runs of various 1-D KiD setups,
### either for Levante HPC or a generic "vanilla" computer.
### It assumes you have already compiled Cleo's python bindings (compile_cleo_python_bindings.sh)

### ------------------ input parameters ---------------- ###
### ------- You may want to edit these lines to -------- ###
### --------- specify your build configuration --------- ###
### ---------------------------------------------------- ###

### ------------- generic input parameters ------------- ###
### _NOTE_: best to use absolute paths here
machine=${1:-"levante"}                                    # machine, either "levante" or "vanilla"
model=${2:-"cleo"}                                         # model, either "cleo" "pysdm" or "bulk"
run_flags=${3:-""}                                         # flags to pass to run_[model].py script
path2sdm1dkid=${4:-"/home/m/m300950/SDM-1D-KiD"}           # path to SDM-1D-KiD repository
path2build=${5:-"/home/m/m300950/SDM-1D-KiD/build"}        # path to build directory (will be created)
python=${6:-"/home/m/m300950/SDM-1D-KiD/.venv/bin/python"} # path to sdm1dkid python interpreter
### ---------------------------------------------------- ###

### --------- extra run parameters for Cleo ---------- ###
cleo_ensemble_name="testing"
cleo_ensemble_size=2
### ---------------------------------------------------- ###

### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ------------------ check arguments ----------------- ###
if [[ "${machine}" != "levante" && "${machine}" != "vanilla" ]]; then
    echo "Error: ${machine} machine not recognised. Please choose either 'levante' or 'vanilla'"
    exit 1
fi

if [[ "${model}" != "cleo" && "${model}" != "pysdm" && "${model}" != "bulk" ]]; then
    echo "Error: ${model} model not recognised. Please choose either 'cleo', 'pysdm' or 'bulk'"
    exit 1
fi
### ---------------------------------------------------- ###

### --- set relevant compiler, packages and paths ------ ###
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/helpers" &>/dev/null && pwd)/${machine}.sh"
load_runtime_packages
### ---------------------------------------------------- ###

### ----------------- print build bindings ------------- ###
echo "### --------------- Build Inputs -------------- ###"
echo "machine=${machine}"
echo "python=${python}"
echo "running model=${model}"
echo "run flags=${run_flags}"
echo "SDM 1-D KiD directory=${path2sdm1dkid}"
echo "build directory=${path2build}"

echo "Cleo ensemble name=${cleo_ensemble_name}"
echo "Cleo ensemble size=${cleo_ensemble_size}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### ------------------- run Cleo 1-D KiD --------------- ###
if [[ "${model}" == "cleo" ]]; then
  cmd="${python} ${path2sdm1dkid}/scripts/run_cleo.py \
    ${path2build} \
    ${cleo_ensemble_name} \
    ${cleo_ensemble_size} \
    ${run_flags}"
  echo ${cmd}
  eval ${cmd}
fi
### ---------------------------------------------------- ###
