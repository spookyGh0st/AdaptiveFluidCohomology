file(GLOB_RECURSE SOURCE_FILES CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(eider ${SOURCE_FILES})
set_target_properties(eider PROPERTIES LINKER_LANGUAGE CXX)
target_compile_features(eider PUBLIC cxx_std_23)

target_link_libraries(eider PUBLIC geometry-central)

set_target_properties(eider PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
target_include_directories(eider INTERFACE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(MSVC)
        target_compile_options(eider PRIVATE /W4)
else()
        target_compile_options(eider PRIVATE -Wall -Wextra -Wpedantic)
endif()
