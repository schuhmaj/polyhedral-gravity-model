# Adapted from https://github.com/esa/LADDS/blob/Integrate-Breathe-%26-Sphinx/cmake/modules/ladds_doxygen.cmake
# Include Doxygen
find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

if (DOXYGEN_FOUND AND SPHINX_FOUND)
    # Collect the source/ header files
    file(GLOB_RECURSE POLYHEDRAL_GRAVITY_DOC
            "${PROJECT_SOURCE_DIR}/src/*.h"
            "${PROJECT_SOURCE_DIR}/src/*.cpp"
            )
    # Settings
    set(DOXYGEN_INPUT_DIRECTORY ${POLYHEDRAL_GRAVITY_DOC})
    set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIRECTORY}/html/index.html)
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_XML YES)
    set(DOXYGEN_GENERATE_LATEX NO)
    set(DOXYGEN_HAVE_DOT YES)

    doxygen_add_docs(doc_doxygen
            ${DOXYGEN_INPUT_DIRECTORY}
            COMMENT "Generating doxygen documentation")

    # Generate documentation
    add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
            COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile.doc_doxygen
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            MAIN_DEPENDENCY Doxyfile.doc_doxygen
            COMMENT "Generating docs with doxygen")

    add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

    add_custom_target(Sphinx ALL
            COMMAND ${SPHINX_EXECUTABLE} -b html
            -Dbreathe_projects.polyhedral-gravity-model=${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml
            ${SPHINX_SOURCE} ${SPHINX_BUILD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating documentation with Sphinx")

else ()
    message(WARNING "Doxygen or Sphinx not found. No documentation will be generated. Please install doxygen, dot,
    sphinx and breathe.")
    set(BUILD_POLYHEDRAL_GRAVITY_DOCS OFF CACHE BOOL "" FORCE)
endif ()
