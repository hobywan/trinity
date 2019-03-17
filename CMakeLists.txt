cmake_minimum_required(VERSION 3.9)

# must be done before call to 'project'
get_cmake_property(vars CACHE_VARIABLES)
foreach (var ${vars})
    get_property(currentHelpString CACHE "${var}" PROPERTY HELPSTRING)
    if ("${currentHelpString}" MATCHES "No help, variable specified on the command line." OR "${currentHelpString}" STREQUAL "")
        # message("${var} = [${${var}}]  --  ${currentHelpString}") # uncomment to see the variables being processed
        list(APPEND CL_ARGS "-D${var}=${${var}}")
    endif ()
endforeach ()

# check supported OS
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif ()

# set project infos
project(trinity VERSION 1.0.3)

# program compile-time options
option(Build_Medit "Build medit viewer" ON)
option(Deferred    "Use deferred topology updates scheme" OFF)

# define project infos
enable_language(CXX)

set(HEADERS sources/coarsening.h
            sources/core.h
            sources/hashtable.h
            sources/hashtable.impl.h
            sources/header.h
            sources/hessian.h
            sources/hwcount.h
            sources/io.h
            sources/matching.h
            sources/mesh.h
            sources/metric.h
            sources/numeric.h
            sources/optparse.h
            sources/partition.h
            sources/random_engine.h
            sources/refinement.h
            sources/rmat.h
            sources/smoothing.h
            sources/swapping.h
            sources/sync.h
            sources/timer.h
            sources/tools.h)

set(SOURCES sources/coarsening.cpp
            sources/hessian.cpp
            sources/hwcount.cpp
            sources/io.cpp
            sources/main.cpp
            sources/matching.cpp
            sources/mesh.cpp
            sources/metric.cpp
            sources/numeric.cpp
            sources/partition.cpp
            sources/refinement.cpp
            sources/rmat.cpp
            sources/smoothing.cpp
            sources/swapping.cpp
            sources/sync.cpp
            sources/tools.cpp)

# create library
add_library(trinity SHARED ${HEADERS} ${SOURCES})

target_include_directories(trinity PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/sources>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)

set_target_properties(trinity PROPERTIES
        PUBLIC_HEADER "${HEADERS}"
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
endif()

if (PAPI_FOUND)
    target_include_directories(trinity PUBLIC ${PAPI_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${PAPI_LIBRARIES})
endif ()

if (HWLOC_FOUND)
    target_include_directories(trinity PUBLIC ${HWLOC_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${HWLOC_LIBRARIES})
endif ()


# manage compilation options
target_compile_features(trinity PUBLIC cxx_std_14)

# set build-type dependent options
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
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
if (Deferred)
    target_compile_definitions(trinity PRIVATE DEFERRED_UPDATES)
endif ()

# build medit if necessary
if (Build_Medit)
    include(ExternalProject)
    ExternalProject_Add(medit
            GIT_REPOSITORY https://github.com/ISCDtoolbox/Medit.git
            GIT_TAG master
            PREFIX ${CMAKE_BINARY_DIR}/aux
            INSTALL_COMMAND ""
            CMAKE_ARGS ${CL_ARGS})

    add_dependencies(trinity medit)
    add_custom_command(
            TARGET trinity POST_BUILD
            COMMAND cmake -E copy_if_different aux/src/medit-build/medit bin/medit)
endif ()


# export target from build/install trees
export(TARGETS trinity FILE trinity-config.cmake)
export(PACKAGE trinity)

install(TARGETS trinity
        EXPORT trinity-targets
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/trinity
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(EXPORT trinity-targets
        FILE trinity-config.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/trinity)


# create the command-line tool
add_executable(main)
set_target_properties(main
        PROPERTIES
        OUTPUT_NAME "trinity"
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

target_sources(main PRIVATE sources/main.cpp)
target_link_libraries(main PUBLIC trinity)