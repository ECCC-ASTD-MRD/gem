include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()


message(STATUS "COMPILER=${COMPILER}")
if(COMPILER STREQUAL "gfortran")
   message(STATUS "COMPILER is gfortran")

#   find_package(MPI REQUIRED)
#   find_package(BLAS REQUIRED)
   find_package(LAPACK REQUIRED)
endif()

foreach(_varName BLAS_FOUND BLAS_LINKER_FLAGS BLAS_LIBRARIES BLAS95_LIBRARIES BLAS95_FOUND  BLA_STATIC BLA_VENDOR BLA_F95)
   message(STATUS "${_varName}=${${_varName}}")
endforeach()

message("LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")

# We get the following warning with the Upstream::<lib> syntax:

# Policy CMP0028 is not set: Double colon in target name means ALIAS or
# IMPORTED target.  Run "cmake --help-policy CMP0028" for policy details.
# Use the cmake_policy command to set the policy and suppress this warning.
# 
# Target "maingemdm" links to target "Upstream::LAPACK" but the target was
# not found.  Perhaps a find_package() call is missing for an IMPORTED
# target, or an ALIAS target is missing?
# This warning is for project developers.  Use -Wno-dev to suppress it.

add_executable(maingemdm maingemdm.F90)
target_link_libraries(maingemdm
   gem
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(maingemdm gem phy modelutils vgrid rpncomm rmn massv)


add_executable(maingemgrid maingemgrid.F90)
target_link_libraries(maingemgrid
   gem
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(maingemgrid gem phy modelutils vgrid rpncomm rmn massv)

install(TARGETS maingemdm maingemgrid RUNTIME 
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
