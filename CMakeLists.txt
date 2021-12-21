cmake_minimum_required(VERSION 3.20)
project(polyhedralGravity)

set(CMAKE_CXX_STANDARD 17)

#Appends the the module path to contain additional CMake modules for this project
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(mainFile ${PROJECT_SOURCE_DIR}/src/main.cpp)

add_executable(${PROJECT_NAME} "${mainFile}")