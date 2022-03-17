# polyhedral-gravity-model-cpp

![Build and Test](https://github.com/schuhmaj/polyhedral-gravity-model-cpp/actions/workflows/ctest.yml/badge.svg)

Implementation of the Polyhedral Gravity Model in C++ 17.

The implementation is based of the paper [Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11.](http://dx.doi.org/10.1190/geo2010-0334.1) and its corresponding implementation in FORTRAN.


## Requirements
The project uses the following dependencies:

- GoogleTest 1.11.0 (only required for testing, automatically set-up via CMake)
- spdlog 1.9.2 (required for logging, automatically set-up via CMake)
- tetgen 1.6 (required for I/O, automatically set-up via CMake)
- yaml-cpp 0.7.0 (required for I/O, automatically set-up via CMake)

## Build
The program is build by using CMake. So first make sure that you installed
CMake and then follow these steps:

    mkdir build
    cd build
    cmake ..
    make

## Execution
After the build, the gravity model can be run by executing:

    ./polyhedralGravity [YAML-Configuration-File]

where the YAML-Configuration-File contains the required parameters.
Examples for Configuration Files can be found in this repository
in the folder `/example-config/`.

## Testing
The project uses GoogleTest for testing. In oder to execute those
tests just execute the following command in the build directory:

    ctest
