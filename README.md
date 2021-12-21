# polyhedral-gravity-model-cpp
Implementation of the Polyhedral Gravity Model in C++ 17.

The implementation is based of the paper [Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11.](http://dx.doi.org/10.1190/geo2010-0334.1) and its corresponding implementation in FORTRAN.


## Build
The program is build by using CMake. So first make sure that you installed
CMake and then follow these steps:

    mkdir build
    cd build
    cmake ..
    make

## Execution
After the build, the gravity model can be run by executing:

    ./polyhedralGravity

with TODO arguments.

## Testing
TODO Add GoogleTest