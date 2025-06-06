# If possible, please use Fetchcontent  as described here:
# https://cmake.org/cmake/help/latest/guide/using-dependencies/index.html#fetchcontent-and-find-package-integration
# If not possible, please vendor libraries (aka copy them here to git)
add_subdirectory("geometry-central")
add_subdirectory("polyscope")
set(TBB_TEST OFF)
add_subdirectory("oneTBB")

# Suppress warnings for GCC/Clang/MSVC
function(disable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /w)
    else()
        target_compile_options(${target} PRIVATE -w)
    endif()
endfunction()

# Apply to the relevant targets
disable_warnings(geometry-central)
disable_warnings(polyscope)
disable_warnings(tbb)
