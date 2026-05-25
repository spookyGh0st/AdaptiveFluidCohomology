cmake_minimum_required(VERSION 3.22)
project(eider)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mno-omit-leaf-frame-pointer")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-omit-leaf-frame-pointer")


# After all this time, I still don't understand why you would use shared libs
# see also https://lore.kernel.org/lkml/CAHk-=whs8QZf3YnifdLv57+FhBi5_WeNTG1B-suOES=RcUSmQg@mail.gmail.com/
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

add_subdirectory(external)
add_subdirectory(eider)
add_subdirectory(test)

add_custom_target(generate_evaluation
        # The COMMAND to execute. We use absolute paths to be safe.
        COMMAND
        ${CMAKE_SOURCE_DIR}/external/gtest-parallel/gtest-parallel
        $<TARGET_FILE:eider-test> "--gtest_filter=EvaluatorTest.*"
        DEPENDS eider-test
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating Evaluations by running Evaluator Test suite..."

        # VERBATIM is good practice to ensure arguments with special
        # characters are passed correctly.
        VERBATIM
)

add_custom_target(generate_videos
        # The COMMAND to execute. We use absolute paths to be safe.
        COMMAND
        ${CMAKE_SOURCE_DIR}/external/gtest-parallel/gtest-parallel
        $<TARGET_FILE:eider-test> "--gtest_filter=VideoTest.*"
        DEPENDS eider-test
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating videos by running VideoTest suite..."
        VERBATIM
)

