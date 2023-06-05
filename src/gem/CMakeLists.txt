message(STATUS "(EC) Generating maingem* Makefile")

if ((WITH_MACH) OR (CACHE WITH_MACH) OR (ENV{WITH_MACH}))
  message(STATUS "(EC) chm_LIBRARIES= ${chm_LIBRARIES}")
  message(STATUS "(EC) List of GEMMOD object files: ${GEMMOD_OBJ_FILES}")
endif()

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILER_FLAGS}")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILER_FLAGS}")
endif()

add_executable(maingemdm maingemdm.F90)

# Debug: generate linker map file
target_link_options(maingemdm PRIVATE "-Wl,-Map=maingemdm.map")

target_link_libraries(maingemdm
   ${GEMMOD_OBJ_FILES}
   gemdyn
   ${phy_LIBRARIES}
   ${cpl_LIBRARIES}
   ${chm_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   modelutils_ov_ifort
   vgrid::vgrid
   tdpack::tdpack
   rpncomm::rpncomm
   rmn::rmn-ompi
   OpenMP::OpenMP_Fortran
   MPI::MPI_Fortran
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

install(TARGETS maingemdm RUNTIME DESTINATION bin)
