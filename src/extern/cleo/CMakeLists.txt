cmake_minimum_required(VERSION 3.18.0)

project(CleoBindings
  LANGUAGES CXX C
  DESCRIPTION "Cleo's Python Bindings by Clara Bayley and other Developers"
  VERSION 0.0.0
)

find_package(MPI REQUIRED COMPONENTS C)

message(STATUS "CleoBindings source from CleoBindings_SOURCE_DIR: ${CleoBindings_SOURCE_DIR}")
message(STATUS "CleoBindings build in CleoBindings_BINARY_DIR: ${CleoBindings_BINARY_DIR}")

# ensure C++ compiler uses certain settings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CXX_EXTENSIONS ON)

# skip Cleo's own examples and roughpaper; we only want its libraries and python bindings
set(CLEO_NO_EXAMPLES ON)
set(CLEO_NO_ROUGHPAPER ON)
set(CLEO_NO_PYBINDINGS OFF)

# choose which of Cleo's optional libraries to build
set(CLEO_DOMAIN "cartesian")
set(CLEO_COUPLED_DYNAMICS "numpy")

# fetch and make cleo library available
include(FetchContent)
FetchContent_Declare(
  cleo
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  GIT_REPOSITORY https://github.com/yoctoyotta1024/CLEO.git
  GIT_TAG v0.71.0
)
FetchContent_MakeAvailable(cleo)
message(STATUS "Cleo version v0.71.0, installation in: ${CLEO_BINARY_DIR}")

# print default compiler flags
message(STATUS "CLEO primary CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
