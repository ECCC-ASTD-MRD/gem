#ifdef INTERFACEONLY
      interface
      integer function gmm_rename(old_varname, new_varname)
      character(len=*), intent(in) :: old_varname, new_varname
      end function gmm_rename
      end interface
#endif

#ifndef INTERFACEONLY
   integer function gmm_rename(old_varname, new_varname)
   use gmm_internals
   implicit none

   character(len=*), intent(in) :: old_varname, new_varname
#include "gmm_directory_interface.inc"
   integer*8 :: key

   key = GMM_KEY_NOT_FOUND
   call find_directory_entry(new_varname,key)
   if (key >= 0) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
         print *, '(GMM_RENAME) Variable ', trim(new_varname), ' is already defined'
      endif
      gmm_rename = GMM_ERROR
      return
   endif

   key = GMM_KEY_NOT_FOUND
   call find_directory_entry(old_varname,key)
   if (key == GMM_KEY_NOT_FOUND) then
      if (gmm_verbose_level <= GMM_MSG_WARN) then
         print *, '(GMM_RENAME) Variable ', trim(old_varname), ' not defined'
      endif
      gmm_rename = GMM_ERROR
      return
   endif

   directory(cur_page)%entry(cur_entry)%name = new_varname
   if (gmm_verbose_level == GMM_MSG_DEBUG) then
      print *, '(GMM_RENAME) Variable ', trim(old_varname), ' renamed to ', trim(new_varname)
   endif
   gmm_rename = 0
   return
   end  function gmm_rename
#endif
