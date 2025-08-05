message(STATUS "(EC) Generating rpnphy libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

file(GLOB F_FILES main/*.F90 api/*.F90 base/*.F90 series/*.F90 surface/*.F90 utils/*.F90)
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*")
if(WITH_CPL)
  file(GLOB CPLOCN_SOURCES ${CPLOCN_GLOB})
  list(APPEND F_FILES ${CPLOCN_SOURCES})
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(PROJECT_INCLUDE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(rpnphy STATIC ${F_FILES})

target_compile_definitions(rpnphy PUBLIC ECCCGEM)
if(WITH_CPL)
  target_compile_definitions(rpnphy PUBLIC HAVE_CPL)
  message(INFO " -- (EC RPNPHY) target_compile_definitions HAVE_CPL")
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
  rpnphy_chm_stubs 
  )

if(WITH_CPL)
  string(REPLACE ":" ";" LIBRARY_PATH_AS_LIST $ENV{LIBRARY_PATH})
  list(APPEND CMAKE_LIBRARY_PATH ${LIBRARY_PATH_AS_LIST})
  find_library(LIBCPLGEMNEMO cplgemnemo REQUIRED)
  target_link_libraries(rpnphy PUBLIC ${LIBCPLGEMNEMO})
endif()

set_target_properties(rpnphy PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(rpnphy ${PROJECT_NAME}_build_info)

file(GLOB F_CHMSTUBS_FILES chm_stubs/*.F90)
add_library(rpnphy_chm_stubs STATIC ${F_CHMSTUBS_FILES})
add_dependencies(rpnphy_chm_stubs ${PROJECT_NAME}_build_info)
target_link_libraries(rpnphy_chm_stubs PUBLIC modelutils)

# Don't install lib and include if not installing GEM
# SUPER_PROJECT_NAME should be have been set in GEM main CMakeLists.txt
if ("${SUPER_PROJECT_NAME}" STREQUAL "gem" OR (NOT SUPER_PROJECT_NAME))
  install(TARGETS rpnphy
    LIBRARY       DESTINATION lib
    ARCHIVE       DESTINATION lib)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
endif()
