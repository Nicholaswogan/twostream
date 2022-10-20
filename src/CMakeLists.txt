add_library(twostream
  twostream.f90
)
target_include_directories(twostream PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(twostream PRIVATE -fcheck=all -Wall)
  endif()
  if ("${CMAKE_BUILD_TYPE}" MATCHES "Release")
    target_compile_options(twostream PRIVATE -funroll-loops -finline-functions)
  endif()
endif()