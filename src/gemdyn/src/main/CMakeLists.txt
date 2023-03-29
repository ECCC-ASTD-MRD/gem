message(STATUS "(EC) Generating gemdyn binaries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gemdyn
   ${phy_LIBRARIES}
   ${cpl_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   vgrid
   tdpack
   rpncomm
   rmn
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)
add_dependencies(checkdmpart gemdyn ${phy_SRC} modelutils ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC})

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gemdyn
   ${phy_LIBRARIES}
   ${cpl_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   vgrid
   tdpack
   rpncomm
   rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)
add_dependencies(gemgrid gemdyn ${phy_SRC} modelutils ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC})

add_executable(prgemnml mainprgemnml.F90)
target_link_libraries(prgemnml
   gemdyn
   ${phy_LIBRARIES}
   ${cpl_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   vgrid
   tdpack
   rpncomm
   rmn
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)
add_dependencies(prgemnml gemdyn ${phy_SRC} modelutils ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC})

install(TARGETS checkdmpart gemgrid prgemnml RUNTIME DESTINATION bin)
