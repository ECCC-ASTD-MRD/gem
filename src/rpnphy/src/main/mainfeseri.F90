program mainfeseri
   use app
   use feseri_mod, only: feseri

#include <rpnphy_build_info.h>

   app_ptr=app_init(0,"feseri"//c_null_char,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBPHY,VERSION//c_null_char)

   call app_start()

   call feseri

   app_status=app_end(-1)
end program mainfeseri
