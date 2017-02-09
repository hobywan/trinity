find_path(METIS_INCLUDE_DIR metis.h
  /usr/local/include
  /usr/include
  /usr/include/metis
  $ENV{METIS_DIR}/include
  )

find_library(METIS_LIBRARY metis
  /usr/local/lib
  /usr/lib
  /usr/lib64
  $ENV{METIS_DIR}/lib
  )

if(METIS_INCLUDE_DIR)
  add_definitions(-DHAVE_METIS)
  if(METIS_LIBRARY)
    set(METIS_LIBRARIES ${METIS_LIBRARY})
    set(METIS_FOUND "YES")
    message(STATUS "Metis found at ${METIS_INCLUDE_DIR}")
  endif(METIS_LIBRARY)
else()
  message(FATAL_ERROR "Metis not found")
endif()

mark_as_advanced( METIS_INCLUDE_DIR METIS_LIBRARY )
