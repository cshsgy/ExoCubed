# configuration for ktable test

macro(SET_IF_EMPTY _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# athena variables
set_if_empty(NUMBER_GHOST_CELLS 2)
set_if_empty(NVAPOR 0)

# canoe variables
set_if_empty(NCLOUD 0)
set_if_empty(NTRACER 2)

# canoe task set(TASKLIST InversionTasks)

# canoe configure
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/opt/")
set(USER_MODS ${CMAKE_SOURCE_DIR}/examples/giants_add_absorber)
set(HYDROSTATIC ON)
set(RFM ON)
set(NETCDF ON)
set(DISORT ON)
