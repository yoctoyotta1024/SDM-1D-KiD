#!/bin/bash

load_build_packages() {
  export cxx_compiler="$(command -v mpic++)"
  export c_compiler="$(command -v mpicc)"
  export cxx_flags="-Werror -Wall -Wextra -pedantic -Wno-unused-parameter -O3"
}

load_runtime_packages() {
}
