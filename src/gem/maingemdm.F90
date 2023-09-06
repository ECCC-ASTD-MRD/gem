program gem
   use app
   use iso_fortran_env
   implicit none

#include <gem_build_info.h>

   integer(kind=int32) ierror

   app_ptr=app_init(0,PROJECT_NAME_STRING,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBVGRID,HAVE_VGRID)
   call app_libregister(APP_LIBTDPACK,HAVE_TDPACK)
   call app_libregister(APP_LIBDYN,dyn_VERSION)
   call app_libregister(APP_LIBPHY,phy_VERSION)
   call app_libregister(APP_LIBMDLUTIL,modelutils_VERSION)

   call MPI_INIT(ierror)
   call app_start()
 
   ! Initialize: Domain, MPI, processor topology and ptopo.cdk
   call init_component()

   ! Establish: model configuration, domain decomposition and model geometry
   call set_world_view()
  
   ! Initialize the ensemble prevision system
   call itf_ens_init()
  
   ! Initialize the physics parameterization package
   call itf_phy_init()
  
   ! Initialize tracers
   call tracers()
  
   ! Setup main memory
   call main_gmm_storage()
   call set_dyn_opr()
  
   ! Run GEM
   call gem_ctrl()
  
   ! Terminate
   call stop_world_view()
!   call MPI_FINALIZE(ierror)  

   app_status=app_end(-1)
end program gem
