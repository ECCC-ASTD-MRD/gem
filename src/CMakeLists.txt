message(STATUS "(EC) Generating rpnphy libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

file(GLOB F_FILES main/*.F90 api/*.F90 base/*.F90 series/*.F90 surface/*.F90 utils/*.F90)
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(PROJECT_INCLUDE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(rpnphy STATIC ${F_FILES})

target_compile_definitions(rpnphy PUBLIC ECCCGEM)

target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRECTORY}>)
target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

target_include_directories(rpnphy PRIVATE
  utils
  base
  )

target_link_libraries(rpnphy PUBLIC
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  vgrid::vgrid
  tdpack::tdpack
  modelutils
  rpnphy_cpl_stubs
  rpnphy_chm_stubs 
  )

set_target_properties(rpnphy PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(rpnphy ${PROJECT_NAME}_build_info)

file(GLOB F_CHMSTUBS_FILES chm_stubs/*.F90)
add_library(rpnphy_chm_stubs STATIC ${F_CHMSTUBS_FILES})
add_dependencies(rpnphy_chm_stubs ${PROJECT_NAME}_build_info)
target_link_libraries(rpnphy_chm_stubs PUBLIC modelutils)

add_library(rpnphy_cpl_stubs STATIC cpl_stubs/cpl_itf_stubs.F90)
add_dependencies(rpnphy_cpl_stubs ${PROJECT_NAME}_build_info)

install(TARGETS rpnphy
    LIBRARY       DESTINATION lib
    ARCHIVE       DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
