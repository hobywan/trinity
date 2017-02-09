find_path(SPMP_INCLUDE_DIR CSR.hpp
  /usr/local/include
  /usr/include
  /usr/include/spmp
  $ENV{SPMP_DIR}
)

find_library(SPMP_LIBRARY spmp
  /usr/local/lib
  /usr/lib
  /usr/lib64
  $ENV{SPMP_DIR}/lib
)

if(SPMP_INCLUDE_DIR)
  add_definitions(-DHAVE_SPMP)
  if(SPMP_LIBRARY)
    set( SPMP_LIBRARIES ${SPMP_LIBRARY})
    set( SPMP_FOUND "YES" )
    message(STATUS "SpMP found at ${SPMP_INCLUDE_DIR}")
  endif(SPMP_LIBRARY)
else()
  message(FATAL_ERROR "SpMP not found")
endif()

mark_as_advanced( SPMP_INCLUDE_DIR SPMP_LIBRARY )
