## set up model configuration ##

if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "DebugRelease")
endif()

# equation of state
set(EquationOfState adiabatic_hydro
  CACHE STRING "Choose the equation of state for primitive-conserved conversion")

# riemann solver
set(RiemannSolver hllc
  CACHE STRING "Choose the Riemann Solver")

# MPI flag
option(UseMPI "Enable MPI" ON)

# MPI flag
if (${UseMPI})
  find_package(MPI REQUIRED)
  set(MPIOption MPI_PARALLEL)
endif()

# CubedSphere flag
option(UseCubedSphere "Enable CubedSphere" ON)

# CubedSphere flag
if (${UseCubedSphere})
  set(CubedSphereOption CUBED_SPHERE)
  set(CoordinateSystem gnomonic_equiangle)
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
