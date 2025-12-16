message(STATUS "(EC) Generating maingem* Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILER_FLAGS}")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILER_FLAGS}")
endif()

add_executable(${gem_exe} maingemdm.F90)

target_compile_definitions(${gem_exe} PUBLIC
  dyn_VERSION="${dyn_VERSION}"
  modelutils_VERSION="${modelutils_VERSION}"
  phy_VERSION="${phy_VERSION}"
)

if (("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM") AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2025))
  # link with libgcc_s to get a traceback
  set(GCC_S gcc_s)
endif()

target_link_libraries(${gem_exe}
   gemdyn
   ${GCC_S}
   ${phy_LIBRARIES}
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

install(TARGETS ${gem_exe} RUNTIME DESTINATION bin)
