message(STATUS "(EC) Generating maingem* Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

if ( (DEFINED ENV{EC_ARCH}) AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 2021) )
  if (NOT DEFINED ENV{CRAYPE_VERSION})
    set(EXTRA_LIBRARIES "-lopen-pal")
  endif()
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/maingemdm.in ${CMAKE_CURRENT_SOURCE_DIR}/maingemdm.F90)

add_executable(maingemdm maingemdm.F90)
target_link_libraries(maingemdm
   gemdyn
   rpnphy
   rpnphy_cpl_stubs
   rpnphy_chm_stubs
   modelutils
   modelutils_tmg_stubs
   modelutils_ov_ifort
   vgrid
   tdpack
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
if (NOT "${gemdyn_SRC}" STREQUAL "" OR NOT "${rpncomm_SRC}" STREQUAL "" OR NOT "${modelutils_SRC}" STREQUAL ""
      OR NOT "${vgrid_SRC}" STREQUAL "" OR NOT "${rpncomm_SRC}" STREQUAL "" OR NOT "${rmn_SRC}" STREQUAL ""  OR NOT "${massv_SRC}" STREQUAL "")
    add_dependencies(maingemdm ${gemdyn_SRC} ${rpnphy_SRC} ${modelutils_SRC} ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC} ${massv_SRC})
endif()
install(TARGETS maingemdm RUNTIME DESTINATION bin)
