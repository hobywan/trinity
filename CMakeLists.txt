cmake_minimum_required(VERSION 3.0)
set (CMAKE_CXX_STANDARD 14)

project(trinity)
set(TRINITY_VERSION_MAJOR 0)
set(TRINITY_VERSION_MINOR 3)

enable_language(CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set(EXECUTABLE_OUTPUT_PATH .)

# CXX_FLAGS
find_package(OpenMP REQUIRED)
set(WARN_FLAGS "-Wall -Wshadow -Wdeprecated -Wunsafe-loop-optimizations -Wcast-align -g -O2")
set(PERF_FLAGS "-O3 -funroll-loops -g -D_REENTRANT")
#set(PERF_FLAGS "-fopt-info-vec -Wpacked -Wpadded ${PERF_FLAGS})

# Mandatory C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=gnu++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14")
else()
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support.
  Please use a different C++ compiler.")
endif()

# Build type dependent flags
if(CMAKE_BUILD_TYPE MATCHES DEBUG)
  set(CMAKE_CXX_FLAGS_DEBUG "${OpenMP_CXX_FLAGS} ${WARN_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_VERBOSE_MAKEFILE ON)
else()
  set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PERF_FLAGS} ${CMAKE_CXX_FLAGS}")
  if(CMAKE_BUILD_TYPE MATCHES BENCHMARK)
    add_definitions(-DBENCHMARK)
  endif()
endif()

# OS dependent flags
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions( -DDARWIN )
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions( -DLINUX )
else()
  message(FATAL_ERROR "Windows not supported.")
endif()


set(TRINITY_LIBRARIES)

# Dependencies
find_package(PAPI)
if(PAPI_FOUND)
  include_directories(${PAPI_INCLUDE_DIR})
  set(TRINITY_LIBRARIES ${PAPI_LIBRARIES} ${TRINITY_LIBRARIES})
endif()

find_package(Hwloc)
if(HWLOC_FOUND)
  include_directories(${HWLOC_INCLUDE_DIR})
  set(TRINITY_LIBRARIES ${HWLOC_LIBRARIES} ${TRINITY_LIBRARIES})
endif()

# add definitions
OPTION(DEFER_FIXES "Use deferred updates mechanism" OFF)
if (DEFER_FIXES)
  add_definitions(-DDEFERRED_UPDATES)
endif()

# Build files
#file(
#  GLOB_RECURSE
#  source_files
#  sources/*.h*
#  sources/*.c*
#)

# set header and libs directories
include_directories(${CMAKE_SOURCE_DIR}/sources)
#link_directories   (${CMAKE_SOURCE_DIR}/libraries)

# create binary and link libs
set(BINARY "trinity")
add_executable(trinity
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

target_link_libraries(trinity ${TRINITY_LIBRARIES})

