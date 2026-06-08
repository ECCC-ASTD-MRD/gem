program maincheckdmpart
   use app
   use app_mpmd
   use iso_fortran_env
   implicit none
   integer  err
   integer :: colors(3), COMMs(3), wnum, wme, cnum, cme
   character(len=256) :: component_S

#include <gemdyn_build_info.h>

   app_ptr=app_init(0,"checkdmpart"//c_null_char,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBDYN,VERSION//c_null_char)

   component_S= 'GEMDM' ; colors= (/1,3,4/)
   call MiMd_init (component_S,colors,3,COMMs,cnum, cme)

   call init_component(COMMs,3)

   call checkdmpart

   app_status=app_end(-1)
   call App_MPMD_Finalize()
   call MPI_FINALIZE(err)

end program maincheckdmpart
