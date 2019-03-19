cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0048 NEW)

# check supported OS
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif ()

# use same build configuration for external dependencies
# need to be called before 'project'
get_cmake_property(Build_Options CACHE_VARIABLES)
set(EXTERN_BUILD_ARGS "-Wno-dev")
foreach (cached ${Build_Options})
    get_property(helper CACHE "${cached}" PROPERTY HELPSTRING)
    if ("${helper}" MATCHES "No help, variable specified on the command line." OR "${helper}" STREQUAL "")
        list(APPEND EXTERN_BUILD_ARGS "-D${cached}=${${cached}}")
    endif ()
endforeach ()

# set project infos
project(trinity VERSION 1.0.3 LANGUAGES CXX)

# program compile-time options
option(Build_Medit    "Build medit viewer" ON)
option(Build_GTest    "Build googletest" ON)
option(Build_Main     "Build command-line tool" ON)
option(Build_Examples "Build examples sources" ON)
option(Use_Deferred   "Use deferred topology updates scheme" OFF)

# create library
add_library(trinity SHARED)

target_include_directories(trinity PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_sources(trinity PRIVATE
        src/coarsening.cpp
        src/hessian.cpp
        src/hwcount.cpp
        src/io.cpp
        src/main.cpp
        src/matching.cpp
        src/mesh.cpp
        src/metric.cpp
        src/numeric.cpp
        src/partition.cpp
        src/parser.cpp
        src/refinement.cpp
        src/rmat.cpp
        src/smoothing.cpp
        src/swapping.cpp
        src/sync.cpp
        src/tools.cpp)

set_target_properties(trinity PROPERTIES
        PUBLIC_HEADER include/trinity.h
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        CXX_EXTENSIONS ON)


# manage dependencies
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)
find_package(PAPI)
find_package(Hwloc)

if (OpenMP_CXX_FOUND)
    target_link_libraries(trinity PUBLIC OpenMP::OpenMP_CXX)
endif ()

if (PAPI_FOUND)
    target_compile_definitions(trinity PUBLIC HAVE_PAPI)
    target_include_directories(trinity PUBLIC ${PAPI_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${PAPI_LIBRARY})
endif ()

if (HWLOC_FOUND)
    target_compile_definitions(trinity PUBLIC HAVE_HWLOC)
    target_include_directories(trinity PUBLIC ${HWLOC_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${HWLOC_LIBRARY})
endif ()


# manage compilation options
target_compile_features(trinity PUBLIC cxx_std_14)
target_compile_definitions(trinity PRIVATE DEFAULT_INPUT_DIR="${CMAKE_SOURCE_DIR}/examples")
target_compile_definitions(trinity PRIVATE DEFAULT_BUILD_DIR="${CMAKE_BINARY_DIR}")

# set build-type dependent options
set(DEFAULT_BUILD_TYPE "Release")

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_VERBOSE_MAKEFILE ON)
    target_compile_options(trinity PRIVATE -Wall)
    target_compile_options(trinity PRIVATE -Wextra)
else ()
    target_compile_options(trinity PRIVATE -O3)
    target_compile_options(trinity PRIVATE -funroll-loops)
    target_compile_options(trinity PRIVATE -march=native)
    target_compile_options(trinity PRIVATE -D_REENTRANT)
    target_compile_options(trinity PRIVATE -Wunsafe-loop-optimizations)
    target_compile_options(trinity PRIVATE -Wcast-align)
    target_compile_options(trinity PRIVATE -Wpacked)
endif ()

# enable deferred updates by Rokos et al.
if (Use_Deferred)
    target_compile_definitions(trinity PRIVATE DEFERRED_UPDATES)
endif ()

if (Build_Medit)
    include(ExternalProject)
    ExternalProject_Add(medit
            GIT_REPOSITORY https://github.com/ISCDtoolbox/Medit.git
            GIT_TAG master
            PREFIX ${CMAKE_BINARY_DIR}/extern
            INSTALL_COMMAND ""
            CMAKE_ARGS ${EXTERN_BUILD_ARGS})

    ExternalProject_Get_Property(medit BINARY_DIR)
    add_dependencies(trinity medit)
    add_custom_command(TARGET trinity
            POST_BUILD
            COMMAND cmake -E copy_if_different ${BINARY_DIR}/medit bin/medit)
endif ()

if (Build_GTest)
    include(ExternalProject)
    ExternalProject_Add(googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG master
            PREFIX ${CMAKE_BINARY_DIR}/extern
            INSTALL_COMMAND ""
            CMAKE_ARGS ${EXTERN_BUILD_ARGS})

    ExternalProject_Get_Property(googletest BINARY_DIR)
    add_library(gtest STATIC IMPORTED)
    add_library(gtest_main STATIC IMPORTED)
    set_target_properties(gtest
            PROPERTIES
            IMPORTED_LOCATION ${BINARY_DIR}/lib/libgtest.a)
    set_target_properties(gtest_main
            PROPERTIES
            IMPORTED_LOCATION ${BINARY_DIR}/lib/libgtest_main.a)

    add_dependencies(gtest googletest)
    add_dependencies(gtest_main googletest)
    add_dependencies(trinity gtest)
    target_link_libraries(trinity PRIVATE gtest gtest_main)
    #add_test(NAME example_test COMMAND example)
endif ()

if (Build_Examples)
    set(EXAMPLES examples/adapt.cpp)

    foreach (path ${EXAMPLES})
        string(REPLACE ".cpp" "" trimmed ${path})
        get_filename_component(current ${trimmed} NAME)

        add_executable(${current})
        set_target_properties(${current}
                PROPERTIES
                OUTPUT_NAME "example_${current}"
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

        target_sources(${current} PUBLIC ${CMAKE_SOURCE_DIR}/${path})
        target_link_libraries(${current} PRIVATE trinity)
        target_compile_definitions(${current} PRIVATE DEFAULT_INPUT_DIR="${CMAKE_SOURCE_DIR}/examples")
        target_compile_definitions(${current} PRIVATE DEFAULT_BUILD_DIR="${CMAKE_BINARY_DIR}")
    endforeach ()
endif ()

# build the command-line tool
if (Build_Main)
    add_executable(main)
    set_target_properties(main
            PROPERTIES
            OUTPUT_NAME "trinity"
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

    target_sources(main PRIVATE src/main.cpp)
    target_link_libraries(main PRIVATE trinity)
    add_custom_command(TARGET main
            POST_BUILD
            COMMAND cmake -E make_directory ${CMAKE_BINARY_DIR}/data)
endif ()

# export target from build/install trees
export(TARGETS trinity FILE trinity-config.cmake)
export(PACKAGE trinity)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/trinity"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
        FILES_MATCHING PATTERN "*.h")

install(TARGETS trinity
        EXPORT trinity-targets
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(EXPORT trinity-targets
        FILE trinity-config.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/trinity)