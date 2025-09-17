program mainprgemnml
   use app
   use iso_fortran_env
   implicit none

#include <gemdyn_build_info.h>

   app_ptr=app_init(0,"prgemnml"//c_null_char,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBDYN,VERSION//c_null_char)

   call app_start()

   call prgemnml

   app_status=app_end(-1)

end program mainprgemnml
