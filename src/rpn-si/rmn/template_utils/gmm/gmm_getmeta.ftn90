#ifdef INTERFACEONLY
      interface
      integer function gmm_getmeta(varname, meta)
#include "gmm_definitions.inc"
      character(len=*), intent(in) :: varname
      type(gmm_metadata), intent(out) :: meta
      end function gmm_getmeta
      end interface
#endif

#ifndef INTERFACEONLY
   integer function gmm_getmeta(varname, meta)
   use gmm_internals
   implicit none
   character(len=*), intent(in) :: varname
   type(gmm_metadata), intent(out) :: meta

   include 'gmm_directory_interface.inc'
   integer   :: i
   integer*8 :: key

   key = 0
   call find_directory_entry(varname,key)
!   print *, varname, key
   if (key == GMM_KEY_NOT_FOUND) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
        print *, '(GMM_GETMETA) Variable ', varname, ' not found'
      endif
      gmm_getmeta = GMM_ERROR
      return
   endif

   meta%a = directory(cur_page)%entry(cur_entry)%a
   meta%l = directory(cur_page)%entry(cur_entry)%l
   gmm_getmeta = 0
!
   return
   end  function gmm_getmeta
#endif

