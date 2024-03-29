string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

# 1. Compile plume.cpp

set(namel plume)

add_executable(${namel}.${buildl} ${namel}_pgen.cpp ${namel}_mods.cpp
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
    ${GLOG_LIBRARIES}
    )

target_include_directories(${namel}.${buildl}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CANOE_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_PATH}
    ${NETCDF_INCLUDES}
    ${PNETCDF_INCLUDE_DIR}
    )

# 1. Copy input file to run directory

file(GLOB inputs *.inp *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
