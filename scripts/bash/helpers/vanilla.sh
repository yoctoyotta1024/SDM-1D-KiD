#!/bin/bash

load_build_packages() {
  export cxx_compiler="$(command -v mpic++)"
  export c_compiler="$(command -v mpicc)"
  export cxx_flags="-Werror -Wall -Wextra -pedantic -Wno-unused-parameter -O3"
}

load_runtime_packages() {
  local path2sdm1dkid="$1"
  local path2build="$2"

  # (optional) prepend to python path to make import searches faster
  site_packages="$(uv run python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")"
  path2cleo_python_bindings="${path2build}/_deps/cleo-build/cleo_python_bindings"
  echo "export PYTHONPATH=${site_packages}:${path2cleo_python_bindings}:${path2sdm1dkid}:${PYTHONPATH}"
  export PYTHONPATH=${site_packages}:${path2cleo_python_bindings}:${path2sdm1dkid}:${PYTHONPATH}
}
