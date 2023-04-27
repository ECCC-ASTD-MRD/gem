message(STATUS "(EC) Generating maingem* Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

add_executable(maingemdm maingemdm.F90)
target_include_directories(maingemdm PUBLIC
   ${phy_INCLUDE_DIRS}
   ${dyn_INCLUDE_DIRS}
   ${modelutils_INCLUDE_DIRS})

target_link_libraries(maingemdm
   gemdyn
   ${phy_LIBRARIES}
   ${cpl_LIBRARIES}
   ${chm_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   modelutils_ov_ifort
   vgrid
   tdpack
   rpncomm
   rmn
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

add_dependencies(maingemdm gemdyn ${phy_SRC} modelutils ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC})

install(TARGETS maingemdm RUNTIME DESTINATION bin)
