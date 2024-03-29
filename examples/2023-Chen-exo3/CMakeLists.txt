string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

set(test_cases)

if(${HYDROSTATIC})
  if(${EQUATION_OF_STATE} STREQUAL "shallow_yz")
    list(APPEND test_cases dam_break)
    list(APPEND test_cases test_blast)
    list(APPEND test_cases W92)
    list(APPEND test_cases test_injection)
    list(APPEND test_cases steady_zonal)
  endif()
else()
  list(APPEND test_cases test_adiabat)
  list(APPEND test_cases hs94)
  list(APPEND test_cases hot_jupiter)
  list(APPEND test_cases polar_dry)
endif()

foreach(namel ${test_cases})
  add_executable(${namel}.${buildl} ${namel}.cpp
    ${CMAKE_BINARY_DIR}/main_${TASK}.cpp)

  set_target_properties(
      ${namel}.${buildl}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
      COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${buildu}}
      )

  target_link_libraries(${namel}.${buildl}
      ${CANOE_LIBRARY_${buildu}}
      ${NETCDF_LIBRARIES}
      ${PNETCDF_LIBRARY}
      ${MPI_CXX_LIBRARIES}
      ${GLOG_LIBRARY}
      )

  target_include_directories(${namel}.${buildl}
      PRIVATE
      ${CMAKE_BINARY_DIR}
      ${CANOE_INCLUDE_DIR}
      ${MPI_CXX_INCLUDE_PATH}
      ${EIGEN3_INCLUDE_DIR}
      ${NETCDF_INCLUDES}
      ${PNETCDF_INCLUDE_DIR}
      )
endforeach()


# 1. Copy input file to run directory

file(GLOB inputs *.inp *.dat)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
