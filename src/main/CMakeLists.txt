message(STATUS "(EC) Generating rpnphy binaries Makefile")

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   rpnphy
   rmn::rmn
)
add_dependencies(feseri rpnphy ${modelutils_SRC})
 
add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
   rmn::rmn
)
add_dependencies(prphynml rpnphy ${modelutils_SRC})

install(TARGETS feseri prphynml RUNTIME DESTINATION bin)
