cmake_minimum_required(VERSION 3.18)

project(
  ExoCubed
  LANGUAGES CXX
  )
set(CMAKE_CXX_STANDARD 14)
message(STATUS "")
message(STATUS "== Setting up ExoCubed library ==")

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

## 1. set up compiler flags ##
message(STATUS "1. Set up project compiler flags ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake)

## 2. set up project configuration ##
message(STATUS "2. Set up project configuration flags ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake)

## 3. set up project paths ##
message(STATUS "3. Set up project paths ...")
set(EXECUTABLE_RELDIR "bin")
set(INCLUDE_RELDIR "include")
set(SOURCE_RELDIR "source")

## 4. set up project system libraries ##
message(STATUS "4. Set up project system libraries ...")

## 5. set up project external libraries ##
message(STATUS "5. Set up project external libraries ...")
set(HAS_EXTERNAL_LIBRARY FALSE)

#find_package(Eigen3 3.3 REQUIRED HINTS ${HOME}/opt)
find_package(Eigen3 3.3 REQUIRED)

# compile debugger
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(debugger)
  message(STATUS "== Done.")
  set(HAS_EXTERNAL_LIBRARY TRUE)
endif()

# compile athenapp
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(athenapp)
  message(STATUS "== Done.")
  set(HAS_EXTERNAL_LIBRARY TRUE)
endif()

# compile cliutils
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(cliutils)
  message(STATUS "== Done.")
  set(HAS_EXTERNAL_LIBRARY TRUE)
endif()

## 6. set up project library ##
if (HAS_EXTERNAL_LIBRARY)
  message(STATUS "")
endif()
message(STATUS "6. Set up project complication and installation ...")

# set up configure file
configure_file(${CMAKE_SOURCE_DIR}/configure.hpp.in configure.hpp @ONLY)
add_subdirectory(src)
add_subdirectory(tests)
