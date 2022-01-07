include(FetchContent)

#Fetches the version 1.6 for tetgen
FetchContent_Declare(tetgen
        GIT_REPOSITORY https://github.com/libigl/tetgen.git
        )

FetchContent_MakeAvailable(tetgen)

# Disable warnings from the library target
target_compile_options(tetgen PRIVATE -w)