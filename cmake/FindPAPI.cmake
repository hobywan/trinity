include(FindPackageHandleStandardArgs)

find_path(PAPI_INCLUDE_DIR papi.h)
find_library(PAPI_LIBRARY  papi)
find_package_handle_standard_args(PAPI DEFAULT_MSG PAPI_INCLUDE_DIR PAPI_LIBRARY)

mark_as_advanced(PAPI_INCLUDE_DIR PAPI_LIBRARY)

