include(FetchContent)

message(STATUS "Setting up tetgen")

#Fetches the version 1.6 for tetgen
FetchContent_Declare(tetgen
        GIT_REPOSITORY https://github.com/libigl/tetgen.git
        )

FetchContent_GetProperties(tetgen)
if(NOT tetgen_POPULATED)
    FetchContent_Populate(tetgen)

    # This is ugly, but functional
    # We modify one source file in order to prevent console output from the printf function
    # without polluting the global include path (with would lead to more uglier issues)
    if(NOT EXISTS ${tetgen_SOURCE_DIR}/tetgen_mod.cxx)
        message(STATUS "Creating modified tetgen.cxx in order to prevent console output from library")
        file(READ ${tetgen_SOURCE_DIR}/tetgen.cxx TETGEN_CXX OFFSET 3748)

        file(APPEND ${tetgen_SOURCE_DIR}/tetgen_mod.cxx
                "#include \"tetgen.h\"\n#define printf(fmt, ...) (0)\n"
                "${TETGEN_CXX}"
                )
    endif()

    add_library(tetgen STATIC
            ${tetgen_SOURCE_DIR}/tetgen_mod.cxx
            ${tetgen_SOURCE_DIR}/predicates.cxx
            )

    target_compile_definitions(tetgen PRIVATE -DTETLIBRARY)

    target_include_directories(tetgen INTERFACE "${tetgen_SOURCE_DIR}")

endif()

# Disable warnings from the library target
target_compile_options(tetgen PRIVATE -w)