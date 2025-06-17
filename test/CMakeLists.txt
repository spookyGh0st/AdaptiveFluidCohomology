include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS "*.h" "*.cpp")
add_executable(eider-test ${SOURCE_FILES})
set_target_properties(eider-test PROPERTIES LINKER_LANGUAGE CXX)
target_compile_features(eider-test PUBLIC cxx_std_23)

target_link_libraries(eider-test PUBLIC eider polyscope implot GTest::gtest_main)
# target_compile_options(eider-test PRIVATE -ftime-trace)
include(GoogleTest)
gtest_discover_tests(eider-test)


if(MSVC)
    target_compile_options(eider PRIVATE /W4)
else()
    target_compile_options(eider PRIVATE -Wall -Wextra -Wpedantic)
endif()
