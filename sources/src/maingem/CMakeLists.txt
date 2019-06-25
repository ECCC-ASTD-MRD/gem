include(FindMPI OPTIONAL)
#include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

add_executable(maingemdm maingemdm.F90)
target_link_libraries(maingemdm gem phy modelutils vgrid rpncomm rmn massv
               ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} -lfftw3_omp -lfftw3 )

#target_link_libraries(maingemdm gem phy modelutils
#               -ldescrip -lrpn_comm -lrmnMP -lmassv
#               ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} -lfftw3_omp -lfftw3 )

add_dependencies(maingemdm gem phy modelutils vgrid rpncomm rmn massv)


add_executable(maingemgrid maingemgrid.F90)
target_link_libraries(maingemgrid gem phy modelutils vgrid rpncomm rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(maingemgrid gem phy modelutils vgrid rpncomm rmn massv)

install(TARGETS maingemdm maingemgrid RUNTIME 
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
