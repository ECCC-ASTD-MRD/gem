!---------------------------------- LICENCE BEGIN -------------------------------
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
!---------------------------------- LICENCE END ---------------------------------

!/@
recursive function phygetindx(F_varname_S,F_outname_S,F_inname_S,F_bus_S,F_index,F_param,F_npar) result(F_istat)
   use clib_itf_mod, only: clib_tolower, clib_toupper
   implicit none
   !@objective Return params of var 
   !@arguments
   integer,intent(in) :: F_npar
   character(len=*),intent(inout) :: F_varname_S,F_outname_S,F_inname_S,F_bus_S
   integer,intent(out) :: F_index,F_param(F_npar)
   !@return
   integer :: F_istat
   !@author Michel Desgagne - April 2011
   !@revision
   !  2011-06 Stephane Chamberland 
   !@description
   !  Find matching F_varname_S (or F_outname_S if F_varname_S=='')
   !  return 
   !      RMN_ERR on "no match"
   !      nb params on match (then F_varname_S,F_outname_S,F_bus_S,F_index,F_param are filled with values from phys gesdict provided values)
   !@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   include "buses.cdk"

   integer :: nparams,idx,istat
   character(len=32) :: varname_S,outname_S,inname_S,bus_S
   character(len=32),save :: ext_S = ' '
   !-------------------------------------------------------------------
   F_istat = RMN_ERR
   nparams = min(size(F_param),BUSPAR_MAXPAR)
   F_param = RMN_ERR
   F_index = RMN_ERR

   varname_S = F_varname_S
   if (varname_S /= ' ') varname_S = trim(varname_S)//ext_S
   outname_S = F_outname_S
   inname_S  = F_inname_S
   bus_s = F_bus_S
   istat = clib_tolower(varname_S)
   istat = clib_tolower(outname_S)
   istat = clib_tolower(inname_S)
   istat = clib_toupper(bus_S)

   if (any(bus_S==(/' ','E'/))) then
      idx = priv_getidx(entnm,enttop,varname_S,outname_S,inname_S)
      if (idx > 0) then
         F_bus_S = 'E'
         F_index = entpar(idx,BUSPAR_I0)
         F_param(1:nparams)  = entpar(idx,1:nparams)
      endif
   endif

   if (.not.RMN_IS_OK(F_index) .and. any(bus_S==(/' ','D'/))) then
      idx = priv_getidx(dynnm,dyntop,varname_S,outname_S,inname_S)
      if (idx > 0) then
         F_bus_S = 'D'
         F_index = dynpar(idx,BUSPAR_I0)
         F_param(1:nparams)  = dynpar(idx,1:nparams)
      endif
   endif

   if (.not.RMN_IS_OK(F_index) .and. any(bus_S==(/' ','P'/))) then
      idx = priv_getidx(pernm,pertop,varname_S,outname_S,inname_S)
      if (idx > 0) then
         F_bus_S = 'P'
         F_index = perpar(idx,BUSPAR_I0)
         F_param(1:nparams)  = perpar(idx,1:nparams)
      endif
   endif

   if (.not.RMN_IS_OK(F_index) .and. any(bus_S==(/' ','V'/))) then
     idx = priv_getidx(volnm,voltop,varname_S,outname_S,inname_S)
      if (idx > 0) then
         F_bus_S = 'V'
         F_index = volpar(idx,BUSPAR_I0)
         F_param(1:nparams)  = volpar(idx,1:nparams)
      endif
   endif

   if (RMN_IS_OK(F_index)) then
      if (F_varname_S == ' ') F_varname_S = varname_S !TODO: should we remove ',w'
      if (F_outname_S == ' ') F_outname_S = outname_S
      if (F_inname_S == ' ') F_inname_S  = inname_S
      F_istat = nparams
   else
      if (ext_S == ' ' .and. F_varname_S /= ' ')  then
         ext_S = ',w'
         F_istat = phygetindx(F_varname_S,F_outname_S,F_inname_S,F_bus_S,F_index,F_param,F_npar)
      endif
   endif
   ext_S = ' '
   !-------------------------------------------------------------------
   return

contains

   function priv_getidx(my_busnm_S,my_bustop,my_varname_S,my_outname_S,my_inname_S) result(my_idx)
      implicit none
      integer,intent(in) :: my_bustop
      character(len=*),intent(in) :: my_busnm_S(:,:)
      character(len=*),intent(inout) :: my_varname_S,my_outname_S,my_inname_S
      integer :: my_idx

      integer :: i,istat
      character(len=32) :: str1_S,str2_S,str3_S
      logical :: vn_L,on_L,in_L
      !-------------------------------------------------------------------
      my_idx = RMN_ERR
      do i=1,min(my_bustop,size(my_busnm_S,1))
         str1_S = my_busnm_S(i,BUSNM_VN)
         str2_S = my_busnm_S(i,BUSNM_ON)
         str3_S = my_busnm_S(i,BUSNM_IN)
         istat = clib_tolower(str1_S)
         istat = clib_tolower(str2_S)
         istat = clib_tolower(str3_S)
         vn_L = (my_varname_S == str1_S .or. my_varname_S == ' ')
         on_L = (my_outname_S == str2_S .or. my_outname_S == ' ')
         in_L = (my_inname_S  == str3_S .or. my_inname_S == ' ')
         if (vn_L .and. on_L .and. in_L) then
            my_varname_S = str1_S
            my_outname_S = str2_S
            my_inname_S  = str3_S
            my_idx = i
            return
         endif
      end do
      !-------------------------------------------------------------------
      return
   end function priv_getidx

end function phygetindx

