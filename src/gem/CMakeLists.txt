message(STATUS "(EC) Generating gemdm Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

if (DEFINED ENV{EC_ARCH})
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
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(maingemdm gemdyn rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)

install(TARGETS maingemdm RUNTIME 
  DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
