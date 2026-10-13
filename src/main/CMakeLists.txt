message(STATUS "(EC) Generating rpnphy binaries Makefile")

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   rpnphy
   rmn::rmn-ompi
)
 
add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
   rmn::rmn-ompi
)

add_executable(test_phyutils maintest_phyutils.F90)
target_link_libraries(test_phyutils
   rpnphy
   rmn::rmn-ompi
)

install(TARGETS feseri prphynml test_phyutils RUNTIME DESTINATION bin)
