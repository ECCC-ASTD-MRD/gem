program maingemgrid
   use app
   use iso_fortran_env
   implicit none
   integer  err
   integer :: colors(3), COMMs(3), wnum, wme, cnum, cme
   character(len=256) :: component_S

#include <gemdyn_build_info.h>

   component_S= 'GEMDM' ; colors= (/1,3,4/)
   call MiMd_init (component_S,colors,3,COMMs,wnum, wme, cnum, cme)

   app_ptr=app_init(0,"gemgrid"//c_null_char,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBDYN,VERSION//c_null_char)
 
   call app_start()

   call init_component(COMMs,3)

   call gemgrid

   app_status=app_end(-1)
   call rpn_comm_FINALIZE(err)

end program maingemgrid
