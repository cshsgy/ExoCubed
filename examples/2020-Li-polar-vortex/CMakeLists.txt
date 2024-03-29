string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# 1. Compile straka.cpp

#set(namel straka)

foreach(namel shallow_xy polar_vortex)
  find_package(Boost REQUIRED)

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
      )

  target_include_directories(${namel}.${buildl}
      PRIVATE
      ${CMAKE_BINARY_DIR}
      ${CANOE_INCLUDE_DIR}
      ${MPI_CXX_INCLUDE_PATH}
      ${NETCDF_INCLUDES}
      ${PNETCDF_INCLUDE_DIR}
      )
endforeach()


# 1. Copy input file to run directory

file(GLOB inputs *.inp *.dat)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
