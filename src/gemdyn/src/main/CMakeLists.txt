message(STATUS "(EC) Generating gemdyn binaries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILER_FLAGS}")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILER_FLAGS}")
endif()

if (("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM") AND (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2025))
  # link with libgcc_s to get a traceback
  set(GCC_S gcc_s)
endif()

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   vgrid::vgrid
   tdpack::tdpack
   rpncomm::rpncomm
   rmn::rmn-ompi
   OpenMP::OpenMP_Fortran
   MPI::MPI_Fortran
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${GCC_S}
)

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   vgrid::vgrid
   tdpack::tdpack
   rpncomm::rpncomm
   rmn::rmn-ompi
   OpenMP::OpenMP_Fortran
   MPI::MPI_Fortran
   ${GCC_S}
)

add_executable(prgemnml prgemnml.F90 mainprgemnml.F90)
target_link_libraries(prgemnml
   gemdyn
   modelutils
   modelutils_tmg_stubs
   rmn::rmn
   ${GCC_S}
)

add_executable(serverOUT mainserverout.F90)
target_link_libraries(serverOUT
  gemdyn
  modelutils
  modelutils_tmg_stubs
  rmn::rmn-ompi
  vgrid::vgrid
  tdpack::tdpack
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  MPI::MPI_Fortran
  MPI::MPI_C
  ${GCC_S}
  )

add_executable(serverIN mainserverIN.F90)
target_link_libraries(serverIN
  gemdyn
  modelutils
  modelutils_tmg_stubs
  rmn::rmn-ompi
  vgrid::vgrid
  tdpack::tdpack
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  MPI::MPI_Fortran
  MPI::MPI_C
  ${GCC_S}
  )

add_executable(check_usr_int_dic main_check_usr_int_dic.F90)
target_link_libraries(check_usr_int_dic
  gemdyn
  ${GCC_S}
  )

install(TARGETS checkdmpart gemgrid prgemnml serverOUT serverIN RUNTIME DESTINATION bin)
