set(WITH_SYSTEM_RPN FALSE CACHE BOOL "Control whether to use installed RPN-SI libs")
if ($ENV{WITH_SYSTEM_RPN})
  message(STATUS "(EC) Using RPN-SI installed libraries = $ENV{WITH_SYSTEM_RPN}")
else()
  message(STATUS "(EC) Using RPN-SI installed libraries = ${WITH_SYSTEM_RPN}")
endif()

# Use installed RPN libraries (WITH_SYSTEM_RPN=TRUE)
if ((WITH_SYSTEM_RPN) OR (CACHE WITH_SYSTEM_RPN) OR ($ENV{WITH_SYSTEM_RPN}))
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

if (EXISTS ${CMAKE_SOURCE_DIR}/src/rpnphy/CMakeLists.txt)
  set(phy_LIBRARIES rpnphy CACHE STRING "rpnphy library")
  set(cpl_LIBRARIES rpnphy_cpl_stubs CACHE STRING "rpnphy coupling library")
  set(chm_LIBRARIES rpnphy_chm_stubs CACHE STRING "rpnphy chemistry stub library")
  add_subdirectory(rpnphy rpnphy)
endif()

set(WITH_MACH FALSE CACHE BOOL "Control whether to compile with mach")
if ((WITH_MACH) OR (CACHE WITH_MACH) OR ($ENV{WITH_MACH}))
  # Add gem-mach source files, if source directory exists
  if (EXISTS ${CMAKE_SOURCE_DIR}/src/mach/CMakeLists.txt)
    set(chm_LIBRARIES mach CACHE STRING "mach chemistry library" FORCE)
    add_subdirectory(mach mach)
  else()
    message(ERROR " -- (EC) You set WITH_MACH=TRUE, but no mach directory was found!")
  endif()
endif()

add_subdirectory(gemdyn gemdyn)
add_subdirectory(gem gem)

# Add those tools only for external users, or if WITH_SYSTEM_RPN is FALSE
if ((NOT DEFINED ENV{EC_ARCH} ) OR (NOT ((WITH_SYSTEM_RPN) OR (CACHE WITH_SYSTEM_RPN) OR ($ENV{WITH_SYSTEM_RPN}))))
  add_subdirectory(rpn-si/cclargs cclargs)
  add_subdirectory(rpn-si/fst-tools fst-tools)
  add_subdirectory(rpn-si/filetype filetype)
endif()

# To be installed only with package
install(DIRECTORY ../scripts/ DESTINATION scripts USE_SOURCE_PERMISSIONS PATTERN "rpy" EXCLUDE PATTERN "support" EXCLUDE)
install(DIRECTORY ../share/gem-maestro/ DESTINATION share/gem-maestro USE_SOURCE_PERMISSIONS PATTERN "version.in" EXCLUDE)
install(DIRECTORY ../share/configurations/ DESTINATION share/configurations USE_SOURCE_PERMISSIONS)

# Generate version specific nml
install(CODE "
   message(STATUS \"Generating nml\")
   execute_process(
     COMMAND ${CMAKE_SOURCE_DIR}/scripts/gem_nml_mkref ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref)
   execute_process(
     COMMAND ${CMAKE_SOURCE_DIR}/scripts/rpy/rpy.nml_get -v -f ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref
     OUTPUT_FILE ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref.kv
     )
   execute_process(
     COMMAND cut ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref.kv -d= -f1
     OUTPUT_FILE ${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref.k)
     ")

install(FILES "${CMAKE_SOURCE_DIR}/share/nml_ref/gem_settings.${GEM_VERSION}.ref" DESTINATION share/nml_ref)
install(FILES "${CMAKE_SOURCE_DIR}/src/rpnphy/include/physics_input_table_GENERIC" DESTINATION share/rpnphy RENAME physics_input_table)
