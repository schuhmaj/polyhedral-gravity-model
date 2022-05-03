# Include pybind11
include(pybind11)

# Stores the file of the python interface in a variable
set(pythonInterface ${CMAKE_CURRENT_SOURCE_DIR}/python-interface.cpp)

pybind11_add_module(polyhedral_gravity ${pythonInterface})

target_link_libraries(polyhedral_gravity PUBLIC ${PROJECT_NAME}_lib)