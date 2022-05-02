cmake_minimum_required(VERSION 3.18)
project(polyhedralGravity)

set(CMAKE_CXX_STANDARD 17)

# Best Combo: HOST CPP, DEVICE OMP
# Parallelization on the HOST
set(PARALLELIZATION_HOST "CPP" CACHE STRING "Host parallelization chosen by the user
 (CPP= Serial, OMP = OpenMP, TBB = Intel Threading Building Blocks")
set_property(CACHE PARALLELIZATION_HOST PROPERTY STRINGS CPP, OMP, TBB)

# Parallelization on the DEVICE
set(PARALLELIZATION_DEVICE "CPP" CACHE STRING "Device parallelization chosen by the user
 (CPP= Serial, OMP = OpenMP, TBB = Intel Threading Building Blocks, CUDA = Nvida Cuda")
set_property(CACHE PARALLELIZATION_DEVICE PROPERTY STRINGS CPP, OMP, TBB, CUDA)

# Set the Logging Level
set(LOGGING_LEVEL "2" CACHE STRING "Set the Logging level, default (INFO=2), available options:
TRACE=0, DEBUG=1, INFO=2, WARN=3, ERROR=4, CRITICAL=5, OFF=6")
set_property(CACHE PARALLELIZATION_DEVICE PROPERTY STRINGS 0, 1, 2, 3, 4, 5, 6)

add_compile_definitions(SPDLOG_ACTIVE_LEVEL=${LOGGING_LEVEL})

# Appends the the module path to contain additional CMake modules for this project
# and include everything necessary
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(spdlog)
include(tetgen)
include(yaml)
include(thrust)

# Main file, containing the gravity model's driver
set(mainFile ${PROJECT_SOURCE_DIR}/src/main.cpp)

# Including the actual Polyhedral source files
file(GLOB_RECURSE SRC
        "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.h"
        "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.cpp")

add_library(${PROJECT_NAME}_lib ${SRC})

# Adds the include Path PROJECT_SOURCE_DIR/src to the target polyhedralGravity
target_include_directories(${PROJECT_NAME}_lib PUBLIC
        "${PROJECT_SOURCE_DIR}/src"
        )

# Link libaries
target_link_libraries(${PROJECT_NAME}_lib
        spdlog
        tetgen
        yaml-cpp
        )

# Thrust set-up i.e. the parallelization library
thrust_create_target(Thrust HOST ${PARALLELIZATION_HOST} DEVICE ${PARALLELIZATION_DEVICE})
message(STATUS "Set Parallelization: HOST ${PARALLELIZATION_HOST} DEVICE ${PARALLELIZATION_DEVICE}")
target_link_libraries(${PROJECT_NAME}_lib
        Thrust
        )
# This definitions will set the main eval(..) routine to be executed on the DEVICE
if(NOT PARALLELIZATION_DEVICE STREQUAL "CPP")
    add_compile_definitions(DEVICE)
endif()

# Building the standalone Executable
add_executable(${PROJECT_NAME} ${mainFile})

# Link executable with library
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}_lib)

# Option to build tests or not
option(BUILD_POLYHEDRAL_GRAVITY_TESTS "Set to on if the tests should be built (Default: ON)" ON)
# Only build the tests if they are enabled
if (BUILD_POLYHEDRAL_GRAVITY_TESTS)
    message(STATUS "Building the Polyhedral Gravity Tests")
    # Enables CTest, must be in the top-level CMakeList.txt, otherwise it won't work
    enable_testing()

    # Subdirectory where the tests are located
    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif ()