# Checks for compiler features (such as C++14 support) and compiler
# specific bugs that
#   - usually set up further configuration (such as preprocessor
#     definitions)
#   - disable a specific flag for a specific compiler version.
#
# belong the corresponding file:
#
#   ./cmake/checks/check_01_cpu_features.cmake
#   ./cmake/checks/check_01_cxx_features.cmake
#   ./cmake/checks/check_02_compiler_features.cmake
#   ./cmake/checks/check_02_system_features.cmake
#   ./cmake/checks/check_03_compiler_bugs.cmake
#

# General setup for GCC and compilers sufficiently close to GCC:
#
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  set(CMAKE_CXX_FLAGS_RELEASE 
    "-O2 -funroll-loops -funroll-all-loops -fstrict-aliasing"
    )

  set(CMAKE_CXX_FLAGS_DEBUG
    "-g3"
    )
  set(CMAKE_C_FLAGS_RELEASE 
    "-O2 -funroll-loops -funroll-all-loops -fstrict-aliasing"
    )

  set(CMAKE_C_FLAGS_DEBUG
    "-g3"
    )

  #set(CMAKE_Fortran_FLAGS_RELEASE
  #  "-O3"
  #  )
  set(KNOWN_COMPILER TRUE)
endif()

#
# Setup for ICC compiler (version >= 10):
#
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(_flags ${CMAKE_SOURCE_DIR}/cmake/compiler_flags_intel.cmake)
  message(STATUS "Include ${_flags}")
  include(${_flags})
  set(KNOWN_COMPILER TRUE)
endif()

#
# Setup for MSVC compiler (version >= 2012):
#
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(_flags ${CMAKE_SOURCE_DIR}/cmake/compiler_flags_msvc.cmake)
  message(STATUS "Include ${_flags}")
  include(${_flags})
  set(KNOWN_COMPILER TRUE)
endif()

if (NOT KNOWN_COMPILER)
  message(FATAL_ERROR "\n"
    "Unknown compiler!\n"
    "If you're serious about it, set SETUP_DEFAULT_COMPILER_FLAGS=OFF "
    "and set the relevant compiler options by hand.\n\n"
    )
endif()
