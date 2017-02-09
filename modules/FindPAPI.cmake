find_path(PAPI_INCLUDE_DIR papi.h
  /usr/include
  /usr/local/include
  $ENV{PAPI_DIR}/include
)

find_library(PAPI_LIBRARY papi
  /usr/lib64
  /usr/local/lib
  /usr/lib
  $ENV{PAPI_DIR}/lib
)

if(PAPI_INCLUDE_DIR)
#  add_definitions(-DHAVE_PAPI)
  if(PAPI_LIBRARY)
    set( PAPI_LIBRARIES ${PAPI_LIBRARY})
    set( PAPI_FOUND "YES" )
    message(STATUS "PAPI found at ${PAPI_INCLUDE_DIR}")
  endif(PAPI_LIBRARY)
else()
  message(STATUS "PAPI counters not found")
endif()

mark_as_advanced( PAPI_INCLUDE_DIR PAPI_LIBRARY )
