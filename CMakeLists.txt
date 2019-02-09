cmake_minimum_required(VERSION 3.0)
set (CMAKE_CXX_STANDARD 14)

project(trinity)
set(TRINITY_VERSION_MAJOR 0)
set(TRINITY_VERSION_MINOR 3)

enable_language(CXX)
set(CMAKE_MODULE_PATH modules ${CMAKE_MODULE_PATH})
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

#
include_directories(${CMAKE_SOURCE_DIR}/sources)
link_directories   (${CMAKE_SOURCE_DIR}/libraries)

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
#add_definitions(-DDEFERRED_UPDATES)

# Build files
file(
  GLOB_RECURSE
  source_files
  sources/*.h*
  sources/*.c*
)

# Binary file
#set(BINARY)
#if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
#  set(BINARY "trinity_intel")
#elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#  set(BINARY "trinity_llvm")
#else()
#  set(BINARY "trinity_gnu")
#endif()

# create binary and link libs
set(BINARY "trinity")
add_executable(${BINARY} ${source_files})
target_link_libraries(${BINARY} ${TRINITY_LIBRARIES})

