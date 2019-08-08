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
function phygetvarlist2(F_varlist_S,F_varname_S,F_outname_S,F_inname_S, &
     F_bus_S,F_init,F_maxvars,F_shortmatch_L) result(F_nvars)
   use clib_itf_mod, only: clib_tolower, clib_toupper
   use phy_typedef
   implicit none
   !@objective Return list of var matching given conditions
   !@arguments
   integer,intent(in) :: F_init,F_maxvars
   character(len=*),intent(out) :: F_varlist_S(F_maxvars)
   character(len=*),intent(in) :: F_varname_S,F_outname_S,F_inname_S,F_bus_S
   logical,intent(in) :: F_shortmatch_L
   !@return
   integer :: F_nvars
   !@author Stephane Chamberland, 2012-04
   !@description
   !  Find matching gesdict entry
   !  strings are "short matched" (empty string match all)
   !  F_init could be: 0=no init; 1=init; -1=match all
   !  return nb of match items (up to size(F_varlist_S)
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   include "buses.cdk"
   integer :: istat,vnlen,onlen,inlen,iadd
   character(len=PHY_MAXNAMELENGTH) :: varname_S,outname_S,inname_S,bus_s
   character(len=PHY_MAXNAMELENGTH), external :: bvstrip
   !-------------------------------------------------------------------
   F_nvars = 0
   F_varlist_S = ' '
   varname_S = bvstrip(F_varname_S)
   outname_S = bvstrip(F_outname_S)
   inname_S  = bvstrip(F_inname_S)
   iadd = 1
   if (F_shortmatch_L) iadd = 0
   vnlen = min(max(1,len_trim(varname_S)+iadd),PHY_MAXNAMELENGTH)
   onlen = min(max(1,len_trim(outname_S)+iadd),PHY_MAXNAMELENGTH)
   inlen = min(max(1,len_trim(inname_S)+iadd),PHY_MAXNAMELENGTH)
   bus_S = F_bus_S
   istat = clib_tolower(varname_S)
   istat = clib_tolower(outname_S)
   istat = clib_tolower(inname_S)
   istat = clib_toupper(bus_S)

   if (any(bus_S==(/' ','E'/))) &
        call priv_getlist(entnm,enttop,varname_S,outname_S,inname_S,entpar)
   if (any(bus_S==(/' ','D'/))) &
        call priv_getlist(dynnm,dyntop,varname_S,outname_S,inname_S,dynpar)
   if (any(bus_S==(/' ','P'/))) &
        call priv_getlist(pernm,pertop,varname_S,outname_S,inname_S,perpar)
   if (any(bus_S==(/' ','V'/))) &
        call priv_getlist(volnm,voltop,varname_S,outname_S,inname_S,volpar)
   !-------------------------------------------------------------------
   return

contains

   subroutine priv_getlist(my_busnm_S,my_bustop,my_varname_S,my_outname_S,my_inname_S,my_params)
      implicit none
      integer,intent(in) :: my_bustop,my_params(:,:)
      character(len=*),intent(in) :: my_busnm_S(:,:)
      character(len=*),intent(in) :: my_varname_S,my_outname_S,my_inname_S

      integer :: i,istat
      character(len=PHY_MAXNAMELENGTH) :: str1_S,str2_S,str3_S
      logical :: vn_L,on_L,in_L,init_L
      !-------------------------------------------------------------------
      do i=1,min(my_bustop,size(my_busnm_S,1))
         str1_S = bvstrip(my_busnm_S(i,BUSNM_VN))
         str2_S = bvstrip(my_busnm_S(i,BUSNM_ON))
         str3_S = bvstrip(my_busnm_S(i,BUSNM_IN))

         istat = clib_tolower(str1_S)
         istat = clib_tolower(str2_S)
         istat = clib_tolower(str3_S)
         vn_L = (my_varname_S == str1_S(1:vnlen) .or. my_varname_S == ' ')
         on_L = (my_outname_S == str2_S(1:onlen) .or. my_outname_S == ' ')
         in_L = (my_inname_S  == str3_S(1:inlen) .or. my_inname_S == ' ')
         init_L = (F_init == -1 .or. F_init == my_params(i,BUSPAR_INIT))
         if (vn_L .and. on_L .and. in_L .and. init_L .and. F_nvars<F_maxvars) then
            F_nvars = F_nvars + 1
            F_varlist_S(F_nvars) = str1_S
         endif
         if (F_nvars >= F_maxvars) exit
      end do
      !-------------------------------------------------------------------
      return
   end subroutine priv_getlist

end function phygetvarlist2
