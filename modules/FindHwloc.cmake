find_path(HWLOC_INCLUDE_DIR hwloc.h
  /usr/local/include
  /usr/include
  $ENV{HWLOC_DIR}/include
  )

find_library(HWLOC_LIBRARY hwloc
  /usr/local/lib
  /usr/lib
  /usr/lib64
  $ENV{HWLOC_DIR}/lib
  )

if(HWLOC_INCLUDE_DIR)
  add_definitions(-DHAVE_HWLOC)
  if(HWLOC_LIBRARY)
    set(HWLOC_LIBRARIES ${HWLOC_LIBRARY})
    set(HWLOC_FOUND "YES")
    message(STATUS "Hwloc found at ${HWLOC_INCLUDE_DIR}")
  endif(HWLOC_LIBRARY)
else()
  message(FATAL_ERROR "Hwloc not found")
endif()

mark_as_advanced( HWLOC_INCLUDE_DIR HWLOC_LIBRARY )
