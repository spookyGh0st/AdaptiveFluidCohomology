cmake_minimum_required(VERSION 3.24)
project(eider)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# After all this time, I still don't understand why you would use shared libs
# see also https://lore.kernel.org/lkml/CAHk-=whs8QZf3YnifdLv57+FhBi5_WeNTG1B-suOES=RcUSmQg@mail.gmail.com/
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

add_subdirectory(external)
add_subdirectory(eider)
add_subdirectory(test)
