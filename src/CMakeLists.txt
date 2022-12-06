message(STATUS "(EC) Generating rpnphy libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

file(GLOB F_FILES main/*.F90 api/*.F90 base/*.F90 series/*.F90 surface/*.F90 utils/*.F90)
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*" )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )

add_library(rpnphy STATIC ${F_FILES})
target_include_directories(rpnphy PRIVATE
  utils
  base
  ${CMAKE_CURRENT_SOURCE_DIR}/../include 
  ${modelutils_INCLUDE_DIRS}
  ${modelutils_MODULES_DIR}
  ${vgrid_INCLUDE_DIRS}
  ${tdpack_INCLUDE_DIRS}
  ${rmn_INCLUDE_DIRS}
  )

set_target_properties(rpnphy PROPERTIES VERSION ${PROJECT_VERSION} )
add_dependencies(rpnphy rpnphy_cpl_stubs ${PROJECT_NAME}_build_info)
add_dependencies(rpnphy rpnphy_chm_stubs ${PROJECT_NAME}_build_info)

file(GLOB F_CHMSTUBS_FILES chm_stubs/*.F90)
add_library(rpnphy_chm_stubs STATIC ${F_CHMSTUBS_FILES})
add_library(rpnphy_cpl_stubs STATIC cpl_stubs/cpl_itf_stubs.F90)
add_dependencies(rpnphy_chm_stubs ${PROJECT_NAME}_build_info)
add_dependencies(rpnphy_cpl_stubs ${PROJECT_NAME}_build_info)
target_include_directories(rpnphy_chm_stubs PRIVATE
  ${modelutils_INCLUDE_DIRS}
  ${modelutils_MODULES_DIR}
  )

