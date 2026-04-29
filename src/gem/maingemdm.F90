program gem
      
   use app
   use iso_fortran_env
#if defined(HAVE_NEMO) && !defined(HAVE_GOSSIP)
   use iris_mod
#endif
   implicit none

#ifdef HAVE_MACH
#include <gemmach_build_info.h>
#else
#include <gem_build_info.h>
#endif

   integer(kind=int32) ierror
   integer :: colors(3), COMMs(3), wnum, wme, cnum, cme, un_out
   character(len=256) :: component_S
#if defined(HAVE_NEMO) && !defined(HAVE_GOSSIP)
   integer :: model_comm
   integer :: iris_component
   logical :: with_iris
#endif

   app_ptr=app_init(0,PROJECT_NAME_STRING,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_logstream("stdout")
   call app_libregister(APP_LIBVGRID,HAVE_VGRID//c_null_char)
   call app_libregister(APP_LIBTDPACK,HAVE_TDPACK//c_null_char)
   call app_libregister(APP_LIBDYN,dyn_VERSION//c_null_char)
   call app_libregister(APP_LIBPHY,phy_VERSION//c_null_char)
   call app_libregister(APP_LIBMDLUTIL,modelutils_VERSION//c_null_char)
#ifdef WITH_MACH
   call app_libregister(APP_LIBMACH,mach_VERSION//c_null_char)
#endif
#ifdef WITH_NEMO
   call app_libregister(APP_LIBIRIS,iris_VERSION//c_null_char)
#endif

   component_S= 'GEMDM' ; colors= (/1,3,4/)
   call MiMd_init (component_S,colors,3,COMMs,wnum, wme, cnum, cme)

   call app_start()

#if defined(HAVE_NEMO) && !defined(HAVE_GOSSIP)
   iris_component = App_MPMD_GetComponentId("iris")
   with_iris = iris_component >= 0
   if(with_iris) then
       ! call app_libregister(APP_LIBIRIS, iris_VERSION//c_null_char)
       model_comm = iris%init(ismodel=1)
       if(model_comm == 0) then
          write(error_unit,*) "Error in iris%init"
          call MPI_Abort()
       endif
   endif
#endif

   call gemdm (COMMs,cme,3)

#if defined(HAVE_NEMO) && !defined(HAVE_GOSSIP)
   if(with_iris) then
       call iris%model_finalize()
   endif
#endif

   app_status=app_end(-1)

   call MPI_FINALIZE(ierror)  

end program gem
