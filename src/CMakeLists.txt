# Installs ExoCubed library
#
#   library: libexocubed_debug.o
#   library: libexocubed_release.o
#
# Define the following variables
#
#   EXOCUBED_FOUND
#   EXOCUBED_INCLUDE_DIR
#   EXOCUBED_LIBRARY_DEBUG
#   EXOCUBED_LIBRARY_RELEASE
#
# Normal usage would be:
#
#   include_directories( ${EXOCUBED_INCLUDE_DIR})
#   target_link_libraries( ${EXOCUBED_LIBRARY_DEBUG})
#

file(GLOB _src_files
  transform_ox_for_cubed_sphere.cpp
  cubed_sphere_coord.cpp
  setup_input_and_init_mesh.cpp
  gnomonic_equiangle.cpp
  correct_fluxes_cubed_sphere.cpp
  hydro/calculate_fluxes.cpp
  hydro/cubedsphere_fluxes.cpp
  )

foreach(build ${BUILD_TYPES})
  string(TOLOWER ${build} buildl)
  add_library(exocubed_${buildl}
    ${_src_files})

  set_target_properties(exocubed_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
    )

  target_include_directories(exocubed_${buildl}
    SYSTEM PRIVATE
    ${DEBUGGER_INCLUDE_DIR}
    ${ATHENAPP_INCLUDE_DIR}
    ${CLIUTILS_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${CMAKE_BINARY_DIR}
    )

  target_link_libraries(exocubed_${buildl}
    ${DEBUGGER_LIBRARY_${build}}
    ${ATHENAPP_LIBRARY_${build}}
    ${CLIUTILS_LIBRARY_${build}}
    ${MPI_CXX_LIBRARIES}
    )

  set(EXOCUBED_LIBRARY_${build} exocubed_${buildl}
    CACHE STRING 
    "exocubed library ${build}")
  mark_as_advanced(EXOCUBED_LIBRARY_${build})
endforeach()

set(EXOCUBED_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..
  CACHE PATH 
  "exocubed include directory")
set(EXOCUBED_FOUND TRUE CACHE BOOL "FIND EXOCUBED LIBRARY")
mark_as_advanced(EXOCUBED_INCLUDE_DIR)
mark_as_advanced(EXOCUBED_FOUND)
