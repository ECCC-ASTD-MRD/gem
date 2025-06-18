program mainprphynml
   use app

#include <rpnphy_build_info.h>

   app_ptr=app_init(0,"prphynml"//c_null_char,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBPHY,VERSION//c_null_char)

   call app_start()

   call prphynml

   app_status=app_end(-1)

end program mainprphynml
