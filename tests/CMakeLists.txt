# set up the following tests
#
# test_ghost_cells.cpp

file(GLOB _main ${CMAKE_SOURCE_DIR}/src/main.cpp)
file(GLOB _pgens *.cpp)

include_directories(
  SYSTEM 
  ${DEBUGGER_INCLUDE_DIR}
  ${ATHENAPP_INCLUDE_DIR}
  ${CLIUTILS_INCLUDE_DIR}
  ${EXOCUBED_INCLUDE_DIR}
  ${CMAKE_BINARY_DIR}
  )

foreach(_pgen ${_pgens})
  foreach(build ${BUILD_TYPES})
    string(TOLOWER ${build} buildl)
    get_filename_component(_name ${_main} NAME_WE)
    get_filename_component(_prob ${_pgen} NAME_WE)

    add_library(${_prob}.${buildl} OBJECT ${_pgen})

    add_executable(${_name}.${_prob}.${buildl} ${_main})

    add_dependencies(${_name}.${_prob}.${buildl} ${_prob}.${buildl})

    set_target_properties(${_name}.${_prob}.${buildl}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${EXECUTABLE_RELDIR}"
      COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
      )

    target_link_libraries(${_name}.${_prob}.${buildl}
      ${DEBUGGER_LIBRARY_${build}}
      ${ATHENAPP_LIBRARY_${build}}
      ${CLIUTILS_LIBRARY_${build}}
      ${EXOCUBED_LIBRARY_${build}}
      )
  endforeach()
endforeach()

# model input files
file(GLOB _inputs *.inp)
foreach(_input ${_inputs})
  get_filename_component(_name ${_input} NAME)

  configure_file(${_input} ${CMAKE_BINARY_DIR}/${EXECUTABLE_RELDIR}/${_name} COPYONLY)
endforeach()
