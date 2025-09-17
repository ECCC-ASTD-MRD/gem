program maincheckdmpart
   use app
   use iso_fortran_env
   implicit none
   integer  err

#include <gemdyn_build_info.h>

   call MPI_INIT(err)

   app_ptr=app_init(0,"checkdmpart"//c_null_char,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBDYN,VERSION//c_null_char)

   call app_start()

   call init_component()

   call checkdmpart

   app_status=app_end(-1)
   call rpn_comm_FINALIZE(err)

end program maincheckdmpart
