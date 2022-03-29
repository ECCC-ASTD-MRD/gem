# Build the following components

add_subdirectory(fstinfo fstinfo)
# add some tools only outside ECCC
if (NOT DEFINED ENV{EC_ARCH})
  add_subdirectory(editfst editfst)
  add_subdirectory(fstcomp fstcomp)
  add_subdirectory(fststat fststat)
  add_subdirectory(pgsm pgsm)
  add_subdirectory(voir voir)
endif()
