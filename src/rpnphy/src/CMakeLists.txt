message(STATUS "(EC) Generating rpnphy libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

file(GLOB F_FILES main/*.F90 api/*.F90 base/*.F90 series/*.F90 surface/*.F90 utils/*.F90 ${MACH_FILES})
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(PROJECT_INCLUDE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(rpnphy STATIC ${F_FILES})

target_compile_definitions(rpnphy PUBLIC ECCCGEM)
if ((WITH_MACH) OR (CACHE WITH_MACH) OR ($ENV{WITH_MACH}))
   target_compile_definitions(rpnphy PUBLIC HAVE_MACH)
   message(INFO " -- (EC RPNPHY) target_compile_definitions HAVE_MACH")
   target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${MACH_INCLUDE_DIRECTORY}>)
endif()

target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRECTORY}>)
target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

target_include_directories(rpnphy PRIVATE
  utils
  base
  )

target_link_libraries(rpnphy PUBLIC
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  vgrid::vgrid
  tdpack::tdpack
  modelutils
  )

set_target_properties(rpnphy PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(rpnphy ${PROJECT_NAME}_build_info)


install(TARGETS rpnphy
    LIBRARY       DESTINATION lib
    ARCHIVE       DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
