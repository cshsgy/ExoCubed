## set up model configuration ##

if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "DebugRelease")
endif()

# equation of state
set(EquationOfState adiabatic
  CACHE STRING "Choose the equation of state for primitive-conserved conversion")
set_property(CACHE EquationOfState
  PROPERTY STRINGS
  adiabatic
  shallow_water
  )

# coordinate system
set(CoordinateSystem cartesian
  CACHE STRING "Choose the coordinate system of the problem")
set_property(CACHE CoordinateSystem
  PROPERTY STRINGS
  cartesian
  affine
  )

if (${EquationOfState} STREQUAL "shallow_water")
  set(RiemannSolver roe_shallow_water
    CACHE STRING "Choose the Riemann Solver")
  option(Barotropic "Barotropic equation of state" ON)
else()
  set(RiemannSolver hllc
    CACHE STRING "Choose the Riemann Solver")
  option(Barotropic "Barotropic equation of state" OFF)
endif()

# riemann solver
set_property(CACHE RiemannSolver
  PROPERTY STRINGS
  hllc
  roe_shallow_water
  )

# hydrostatic flag
option(Hydrostatic "Turn on hydrostatic assumption" OFF)

# NetCDF output flag
option(UseNetCDF "Enable NetCDF output" OFF)

if (${UseNetCDF})
  find_package(NetCDF REQUIRED)
  set(NetCDFOption NETCDF_OUTPUT)
else()
  set(NetCDFOption NOT_NETCDF_OUTPUT)
endif()

# PNetCDF output flag
option(UsePNetCDF "Enable NetCDF output" OFF)

if (${UsePNetCDF})
  find_package(PNetCDF REQUIRED)
  set(PNetCDFOption PNETCDF_OUTPUT)
else()
  set(PNetCDFOption NOT_PNETCDF_OUTPUT)
endif()

# MPI flag
option(UseMPI "Enable MPI" ON)

# MPI flag
if (${UseMPI})
  find_package(MPI REQUIRED)
  set(MPIOption MPI_PARALLEL)
endif()

# CubedSphere flag
option(UseCubedSphere "Enable CubedSphere" OFF)

# CubedSphere flag
if (${UseCubedSphere})
  set(CubedSphereOption CUBED_SPHERE)
  set(CoordinateSystem gnomonic_equiangle)
endif()

# Affine Coordinates flag
option(UseAffine "Enable Affine Coordinate" OFF)

# Affine Coordinates flag
if (${UseAffine})
  set(UseAffine AFFINE)
  set(CoordinateSystem affine_coordinates)
endif()

# ghost zone size
set(GhostZoneSize 2
  CACHE STRING "Set ghose zone size")

# configure athenapp
#message(STATUS "Include ${CMAKE_SOURCE_DIR}/athenapp/cmake/setup_configure.cmake")
#include(${CMAKE_SOURCE_DIR}/athenapp/cmake/setup_configure.cmake)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  if (NOT "DEBUG" IN_LIST BUILD_TYPES)
    list(APPEND BUILD_TYPES "DEBUG")
  endif()
endif()

if (CMAKE_BUILD_TYPE MATCHES "Release")
  if (NOT "RELEASE" IN_LIST BUILD_TYPES)
    list(APPEND BUILD_TYPES "RELEASE")
  endif()
endif()
