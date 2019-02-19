cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 14)

project(trinity)
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 3)

# set language and paths
enable_language(CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set(EXECUTABLE_OUTPUT_PATH .)

# create binary
add_executable(trinity)

target_sources(trinity PRIVATE
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

# Dependencies
find_package(OpenMP REQUIRED)
find_package(PAPI)
find_package(Hwloc)

target_link_libraries(trinity PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(trinity PRIVATE "${OpenMP_CXX_FLAGS}")

if (PAPI_FOUND)
    target_include_directories(trinity PUBLIC ${PAPI_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${PAPI_LIBRARIES})
endif ()

if (HWLOC_FOUND)
    target_include_directories(trinity PUBLIC ${HWLOC_INCLUDE_DIR})
    target_link_libraries(trinity PUBLIC ${HWLOC_LIBRARIES})
endif ()

# Compilation options
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_BUILD_TYPE MATCHES Debug)
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

    if (CMAKE_BUILD_TYPE MATCHES Benchmark)
        target_compile_definitions(trinity PRIVATE Benchmark=1)
    endif ()
endif ()

# deferred topology updates
OPTION(DEFER "Use deferred updates mechanism" OFF)
if (DEFER)
    target_compile_definitions(trinity PRIVATE DEFERRED_UPDATES)
endif ()