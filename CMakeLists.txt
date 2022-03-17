cmake_minimum_required(VERSION 3.20)
project(polyhedralGravity)

set(CMAKE_CXX_STANDARD 17)

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

# TODO, Make this an option for the user
# CPP = Serial execution, We also set the DEVICE to CPP since the default CUDA should not be a dependency
thrust_create_target(Thrust HOST CPP DEVICE CPP)
target_link_libraries(${PROJECT_NAME}_lib
        Thrust
        )

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