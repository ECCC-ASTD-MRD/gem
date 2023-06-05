message(STATUS "(EC) Generating gemdyn libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

file(GLOB F_FILES */*.F90)
file(GLOB C_FILES  */*.c)
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/gemdyn/src/main/main*" )
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/gemdyn/src/no_mpi/main*" )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )
set(PROJECT_INCLUDE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(gemdyn STATIC ${F_FILES} ${C_FILES})
add_dependencies(gemdyn ${PROJECT_NAME}_build_info)

target_include_directories(gemdyn PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRECTORY}>)
target_include_directories(gemdyn PRIVATE base)

target_link_libraries(gemdyn PUBLIC
  rmn::rmn-ompi
  vgrid::vgrid
  tdpack::tdpack
  rpncomm::rpncomm
  modelutils
  ${phy_LIBRARIES}
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  MPI::MPI_Fortran
  MPI::MPI_C
  )

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  if (EXTRA_CHECKS AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_compile_options(gemdyn PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:Fortran>>:-fcheck=bounds>" )
    message(STATUS "(EC) Extra debug flags for gemdyn: -fcheck=bounds")
  endif()
endif()

install(TARGETS gemdyn
    LIBRARY       DESTINATION lib
    ARCHIVE       DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
