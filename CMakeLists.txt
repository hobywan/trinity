cmake_minimum_required(VERSION 3.4)

# must be done before call to 'project'
get_cmake_property(vars CACHE_VARIABLES)
foreach (var ${vars})
    get_property(currentHelpString CACHE "${var}" PROPERTY HELPSTRING)
    if ("${currentHelpString}" MATCHES "No help, variable specified on the command line." OR "${currentHelpString}" STREQUAL "")
        # message("${var} = [${${var}}]  --  ${currentHelpString}") # uncomment to see the variables being processed
        list(APPEND CL_ARGS "-D${var}=${${var}}")
    endif ()
endforeach ()

# set project infos
project(trinity VERSION 1.0.3)

# program compile-time options
option(Build_Medit "Build medit viewer" ON)
option(Deferred "Use deferred updates mechanism" OFF)

# define project infos
enable_language(CXX)

# create binary
add_library(trinity SHARED
        sources/coarsening.cpp
        sources/coarsening.h
        sources/hashtable.h
        sources/hashtable.impl.h
        sources/header.h
        sources/hessian.cpp
        sources/hessian.h
        sources/hwcount.cpp
        sources/hwcount.h
        sources/io.cpp
        sources/io.h
        sources/main.cpp
        sources/matching.cpp
        sources/matching.h
        sources/mesh.cpp
        sources/mesh.h
        sources/metric.cpp
        sources/metric.h
        sources/numeric.cpp
        sources/numeric.h
        sources/optparse.h
        sources/partition.cpp
        sources/partition.h
        sources/random_engine.h
        sources/refinement.cpp
        sources/refinement.h
        sources/rmat.cpp
        sources/rmat.h
        sources/smoothing.cpp
        sources/smoothing.h
        sources/swapping.cpp
        sources/swapping.h
        sources/sync.cpp
        sources/sync.h
        sources/timer.h
        sources/tools.cpp
        sources/tools.h)

# headers
target_include_directories(trinity PRIVATE ${CMAKE_SOURCE_DIR}/sources)

# check supported OS
if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(FATAL_ERROR "Sorry, Windows is not supported.")
endif ()

# target options
set_target_properties(trinity
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        CXX_EXTENSIONS ON)


# Dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)
find_package(PAPI)
find_package(Hwloc)

target_link_libraries(trinity PUBLIC "${OpenMP_CXX_FLAGS}")
target_compile_options(trinity PUBLIC "${OpenMP_CXX_FLAGS}")


if (PAPI_FOUND)
    target_include_directories(trinity PUBLIC ${PAPI_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${PAPI_LIBRARIES})
endif ()

if (HWLOC_FOUND)
    target_include_directories(trinity PUBLIC ${HWLOC_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${HWLOC_LIBRARIES})
endif ()

# Compilation options
target_compile_features(trinity PUBLIC cxx_std_14)

# set default build type
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

# set compil flags according to build type
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
            SOURCE_DIR medit/src
            BINARY_DIR medit/build
            INSTALL_DIR medit/install
            INSTALL_COMMAND ""
            CMAKE_ARGS ${CL_ARGS})

    add_dependencies(trinity medit)
endif ()

# create an executable for trinity
add_executable(trinity_run)
target_sources(trinity_run PRIVATE sources/main.cpp)
target_link_libraries(trinity_run PUBLIC trinity)