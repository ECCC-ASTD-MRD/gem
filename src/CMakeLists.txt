# Build the following components
add_subdirectory(rpn-si/rmn rmn)
add_subdirectory(rpn-si/mpi-extras mpi-extras)
add_subdirectory(rpn-si/rpncomm rpncomm)
add_subdirectory(rpn-si/massv massv)
add_subdirectory(rpn-si/vgrid vgrid)
add_subdirectory(modelutils modelutils)
add_subdirectory(rpnphy rpnphy)
add_subdirectory(gemdyn gemdyn)
add_subdirectory(gem gem)
add_subdirectory(rpn-si/cclargs cclargs)
add_subdirectory(rpn-si/filetype filetype)
add_subdirectory(rpn-si/fst-tools fst-tools)

# Declare the dependencies of all the components
add_dependencies(modelutils vgrid)
add_dependencies(rpnphy modelutils vgrid)
add_dependencies(gemdyn rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)
add_dependencies(rpnphy_chm_stubs modelutils vgrid rpncomm rmn)
