message(STATUS "(EC) Generating mach library Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )
set(PROJECT_INCLUDE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

include(sourcelist.cmake)

# Prepare a list of the object files to be passed on to compile maingemdm
set(GEMMOD_OBJ_FILES "")
foreach(file ${GEMMODSRC})
  list(APPEND GEMMOD_OBJ_FILES "${file}.o")
endforeach ()
# Prepend path of the binary dir to the list of object files
set(GEMMOD_BINARY_DIR "${CMAKE_BINARY_DIR}/src/mach/mach/CMakeFiles/mach.dir")
list(TRANSFORM GEMMOD_OBJ_FILES PREPEND "${GEMMOD_BINARY_DIR}/")
# Send the list upwards
set(GEMMOD_OBJ_FILES ${GEMMOD_OBJ_FILES} CACHE INTERNAL "GEMMOD object files")

# Debug: uncomment to print the list of object files
message(STATUS "(EC) List of GEMMOD object files = ${GEMMOD_OBJ_FILES}")

add_library (mach STATIC ${CHEMSRC} ${GEMMODSRC})

target_include_directories(mach PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRECTORY}>)
target_include_directories(mach PRIVATE base
  ${CMAKE_SOURCE_DIR}/src/rpnphy/src/base
  ${CMAKE_SOURCE_DIR}/src/rpnphy/src/utils
  )

target_link_libraries(mach PUBLIC
  rmn::rmn-ompi
  vgrid::vgrid
  tdpack::tdpack
  rpncomm::rpncomm
  modelutils
  ${phy_LIBRARIES}
  OpenMP::OpenMP_Fortran
  MPI::MPI_Fortran
  )

target_compile_definitions(mach PUBLIC "$<$<CONFIG:Debug>:DEBUG_OUTPUT>")
