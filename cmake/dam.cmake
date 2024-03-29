# configuration for straka hydrodynamcis

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 3)

# canoe configure
set(EOS "shallow_yz")
set(NETCDF ON)
set(MPI ON)
set(HYDROSTATIC ON)
set(AFFINE ON)
