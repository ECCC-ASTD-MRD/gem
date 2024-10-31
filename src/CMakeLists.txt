if(DEFINED ENV{WITH_SYSTEM_RPN})
  set(WITH_SYSTEM_RPN $ENV{WITH_SYSTEM_RPN} CACHE BOOL "Control whether to use installed RPN-SI libs")
else()
  set(WITH_SYSTEM_RPN FALSE CACHE BOOL "Control whether to use installed RPN-SI libs")
endif()
message(STATUS "(EC) Using RPN-SI installed libraries = ${WITH_SYSTEM_RPN}")

# Use installed RPN libraries (WITH_SYSTEM_RPN=TRUE)
if(WITH_SYSTEM_RPN)
  find_package(rmn ${rmn_REQ_VERSION} COMPONENTS static REQUIRED)
  find_package(tdpack ${tdpack_REQ_VERSION} COMPONENTS static REQUIRED)
  find_package(rpncomm ${rpncomm_REQ_VERSION} REQUIRED)
  find_package(vgrid ${vgrid_REQ_VERSION} COMPONENTS static REQUIRED)
  
# Use sources for RPN libraries (WITH_SYSTEM_RPN=FALSE)
else()
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/rpn-si/rmn/CMakeLists.txt)
    add_subdirectory(rpn-si/rmn rmn)
    set(rmn_FOUND TRUE CACHE BOOL "Control whether rmn library was found")
   else()
     message(FATAL_ERROR "(EC) No rmn directory found!")
  endif()
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/rpn-si/tdpack/CMakeLists.txt)
    add_subdirectory(rpn-si/tdpack tdpack)
    set(tdpack_FOUND TRUE CACHE BOOL "Control whether tdpack was found")
  else()
    message(FATAL_ERROR "(EC) No tdpack directory found!")
  endif()
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/rpn-si/vgrid/CMakeLists.txt)
    add_subdirectory(rpn-si/vgrid vgrid)
    set(vgrid_FOUND TRUE CACHE BOOL "Control whether vgrid was found")
  else()
    message(FATAL_ERROR "(EC) No vgrid directory found!")
  endif()
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/rpn-si/rpncomm/CMakeLists.txt)
    add_subdirectory(rpn-si/rpncomm rpncomm)
    set(rpncomm_FOUND TRUE CACHE BOOL "Control whether rpncomm library was found")
  else()
    message(FATAL_ERROR "(EC) No rpncomm directory found!")
  endif()
endif()

add_subdirectory(modelutils modelutils)

set(WITH_MACH OFF CACHE BOOL "Control whether to build GEM-MACH")
if (WITH_MACH)
  # Add gem-mach source files, if source directory exists
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/mach/CMakeLists.txt)
    add_subdirectory(mach mach)
    set(MACH_FILES ${CMAKE_SOURCE_DIR}/src/mach/src/base/*.F90)
    set(MACH_INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src/mach/include)
    set(MACH_VERSION ${mach_VERSION} CACHE STRING "mach version")
    message(STATUS "(EC) Compiling GEM-MACH")
  else()
    message(ERROR "(EC) You set WITH_MACH=TRUE, but no mach directory was found!")
  endif()
endif()

set(WITH_NEMO OFF CACHE BOOL "Control whether to build with NEMO coupling (\"GEM-NEMO\")")
if (WITH_NEMO)
  # Add NEMO coupling source files, if source directory exists
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/cplocn/)
    set(CPLOCN_GLOB ${CMAKE_SOURCE_DIR}/src/cplocn/*.F90)
    message(STATUS "(EC) Compiling GEM-NEMO")
  else()
    message(ERROR "(EC) You set WITH_NEMO=TRUE, but no cplocn directory was found!")
  endif()
endif()

if (EXISTS ${CMAKE_SOURCE_DIR}/src/rpnphy/CMakeLists.txt)
  set(phy_LIBRARIES rpnphy CACHE STRING "rpnphy library")
  add_subdirectory(rpnphy rpnphy)
  set(phy_VERSION ${rpnphy_VERSION} CACHE STRING "rpnphy version")
endif()

if (EXISTS ${CMAKE_SOURCE_DIR}/src/gemdyn/CMakeLists.txt)
  add_subdirectory(gemdyn gemdyn)
  set(dyn_VERSION ${gemdyn_VERSION} CACHE STRING "gemdyn version")
endif()

add_subdirectory(gem gem)

# Add those tools only for external users, or if WITH_SYSTEM_RPN is FALSE
if((NOT DEFINED ENV{EC_ARCH}) OR (NOT WITH_SYSTEM_RPN))
  add_subdirectory(rpn-si/cclargs cclargs)
  add_subdirectory(rpn-si/fst-tools fst-tools)
  add_subdirectory(rpn-si/filetype filetype)
endif()

# To be installed only with package
install(DIRECTORY ../scripts/ DESTINATION scripts USE_SOURCE_PERMISSIONS PATTERN "rpy" EXCLUDE PATTERN "support" EXCLUDE)
install(DIRECTORY ../share/configurations/ DESTINATION share/configurations USE_SOURCE_PERMISSIONS)

if (EXISTS ${CMAKE_SOURCE_DIR}/share/gem-maestro)
  install(DIRECTORY ../share/gem-maestro/ DESTINATION share/gem-maestro USE_SOURCE_PERMISSIONS PATTERN "version.in" EXCLUDE)
endif()

# Generate version specific nml
install(CODE "
   message(STATUS \"Generating nml\")
   execute_process(
     COMMAND ${CMAKE_SOURCE_DIR}/scripts/gem_nml_mkref ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref
     COMMAND_ERROR_IS_FATAL ANY
)
   execute_process(
     COMMAND ${CMAKE_SOURCE_DIR}/scripts/rpy/rpy.nml_get -v -f ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref
     OUTPUT_FILE ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref.kv
     COMMAND_ERROR_IS_FATAL ANY
     )
   execute_process(
     COMMAND cut ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref.kv -d= -f1
     OUTPUT_FILE ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref.k
     COMMAND_ERROR_IS_FATAL ANY)
     ")

install(FILES "${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref" DESTINATION share/nml_ref)
install(FILES "${CMAKE_SOURCE_DIR}/src/rpnphy/include/physics_input_table_GENERIC" DESTINATION share/rpnphy RENAME physics_input_table)
install(PROGRAMS "${CMAKE_SOURCE_DIR}/src/rpnphy/bin/rpnphy_nml_mkref" DESTINATION scripts)
install(PROGRAMS "${CMAKE_SOURCE_DIR}/src/gemdyn/bin/gemdyn_nml_mkref" DESTINATION scripts)
