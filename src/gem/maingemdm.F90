program gem
  use app
  implicit none

#include <gem_build_info.h>
#include <gemdyn_version.inc>
#include <rpnphy_version.inc>
#include <modelutils_version.inc>

  app_ptr=app_init(0,PROJECT_NAME_STRING,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
  call app_libregister(APP_LIBVGRID,HAVE_VGRID)
  call app_libregister(APP_LIBGEMDYN,GEMDYN_VERSION_S)
  call app_libregister(APP_LIBRPNPHY,RPNPHY_VERSION_S)
  call app_libregister(APP_LIBMDLUTIL,MODELUTILS_VERSION_S)
  call app_start()

  call gemdm

  app_status=app_end(-1)
end program gem
