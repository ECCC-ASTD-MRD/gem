
!---------------------------------- LICENCE BEGIN ------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END --------------------------------

!/@*
subroutine gmmx_name_ci(F_nameout_S,F_namein_S,F_update_L)
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_tolower
   use rmn_gmm
   implicit none
   !@objective Search GMM labels in a case insensitive way
   !@arguments
   logical,intent(in) :: F_update_L
   character(len=*),intent(in) :: F_namein_S
   character(len=*),intent(out):: F_nameout_S
   !@author S.Chamberland, 2012-06
   !*@/

   integer,parameter :: NMAX = 4096
   integer,save :: nkeys = 0
   character(len=GMM_MAXNAMELENGTH),save :: keylist_S(NMAX)
   character(len=GMM_MAXNAMELENGTH) :: str1_S,str2_S
   integer :: istat,ikey
   ! ---------------------------------------------------------------------
   if (nkeys == 0 .or. F_update_L) then !TODO: should we update every time?
      keylist_S = ' '
      nkeys = gmm_keys(keylist_S,'')
   endif
   F_nameout_S = F_namein_S
   str1_S = F_namein_S
   istat = clib_tolower(str1_S)
   do ikey=1,nkeys
      if (keylist_S(ikey) == ' ') exit
      str2_S = keylist_S(ikey)
      istat = clib_tolower(str2_S)
      if (str1_S == str2_S) then
         F_nameout_S = keylist_S(ikey)
         exit
      endif
   enddo
   ! ---------------------------------------------------------------------
   return
end subroutine gmmx_name_ci


!/@*
subroutine gmmx_name_parts(F_name_S,F_prefix_S,F_basename_S,F_time_S,F_ext_S)
   implicit none
!!!#include <arch_specific.hf>
   !@objective Split a GMM name into convential parts "prefix/basename:time,ext"
   !@arguments
   character(len=*),intent(in) :: F_name_S
   character(len=*),intent(out):: F_prefix_S,F_basename_S,F_time_S,F_ext_S
   !@author S.Chamberland, 2012-06
   !*@/
   integer,parameter :: NMAX = 16
   character(len=1),parameter :: SEP_PREFIX = '/'
   character(len=1),parameter :: SEP_TIME = ':'
   character(len=1),parameter :: SEP_EXT = ','
   character(len=32) :: parts_S(NMAX),tmp_S
   integer :: n
   !---------------------------------------------------------------
   parts_S = ' '
   call str_split2list(parts_S,F_name_S,SEP_PREFIX,NMAX)
   F_prefix_S = ' '
   do n=1,NMAX-1
      if (parts_S(n+1) == ' ') exit
      F_prefix_S = trim(F_prefix_S)//trim(parts_S(n))//SEP_PREFIX
   enddo
   tmp_S = parts_S(n)
   call str_split(F_basename_S,F_time_S,tmp_S,SEP_TIME)
   tmp_S = F_time_S
   call str_split(F_time_S,F_ext_S,tmp_S,SEP_EXT)
   if (F_time_S /= ' ') F_time_S = SEP_TIME//F_time_S
   if (F_ext_S /= ' ') F_ext_S = SEP_EXT//F_ext_S
   !---------------------------------------------------------------
   return
end subroutine gmmx_name_parts
