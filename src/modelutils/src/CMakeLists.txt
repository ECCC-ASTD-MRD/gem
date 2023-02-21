message(STATUS "(EC) Generating modelutils librairies Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

file(GLOB C_FILES */*.c)
file(GLOB F_FILES main/flipit.F90 main/yy2global.F90 main/yydecode.F90 main/yyencode.F90 base/*F90 driver/*.F90 utils/*.F90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )

add_library(modelutils STATIC ${F_FILES} ${C_FILES})
add_dependencies(modelutils ${PROJECT_NAME}_build_info)
target_include_directories(modelutils PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${rmn_INCLUDE_DIRS}
  ${vgrid_INCLUDE_DIRS}
  ${rpncomm_INCLUDE_DIRS}
  ${rpncomm_GEN_INCLUDE_DIR}
  ${tdpack_INCLUDE_DIRS}
  )

target_link_libraries(modelutils PUBLIC
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  )

set_target_properties(modelutils PROPERTIES VERSION ${PROJECT_VERSION} )

add_library(modelutils_ov_ifort STATIC ov_ifort/ifort_fpe_handler.F90)
add_dependencies(modelutils_ov_ifort modelutils)
add_dependencies(modelutils_ov_ifort ${PROJECT_NAME}_build_info)

add_library(modelutils_stubs STATIC stubs/drv_dyn_itf_stubs.F90)
target_include_directories(modelutils_stubs PRIVATE ${vgrid_INCLUDE_DIRS} )
add_dependencies(modelutils_stubs modelutils)
add_dependencies(modelutils_stubs ${PROJECT_NAME}_build_info)

add_library(modelutils_tmg_stubs STATIC tmg_stubs/tmg_stubs.F90)
add_dependencies(modelutils_tmg_stubs ${PROJECT_NAME}_build_info)

install(TARGETS modelutils
    PUBLIC_HEADER DESTINATION include
    LIBRARY       DESTINATION lib
    ARCHIVE       DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
