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
subroutine drv_ptopo_ndoms(F_ndomains,F_dom_deb,F_istat)
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_getenv
   use wb_itf_mod, only: wb_put
   use str_mod
   implicit none
   !@arguments
   integer,intent(out) :: F_ndomains,F_dom_deb,F_istat
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
 
   character(len=16) :: ndomains_S,tmp1_S,tmp2_S
   integer :: dom_fin
   !-------------------------------------------------------------------
   call msg(MSG_DEBUG,'[BEGIN] drv_ptopo_ndoms')
   F_istat = clib_getenv('UM_EXEC_NDOMAINS',ndomains_S)
   if (.not.RMN_IS_OK(F_istat)) then
      call msg(MSG_ERROR,'(drv_ptopo_ndoms) Env variable UM_EXEC_NDOMAINS is undefined')
      return
   endif

   F_dom_deb = -1
   dom_fin = -1
   call str_split(tmp1_S,tmp2_S,ndomains_S,':')
   F_istat = str_toint(F_dom_deb,tmp1_S)
   if (RMN_IS_OK(F_istat) .and. tmp1_S==ndomains_S) then
      dom_fin = F_dom_deb
   else
      F_istat = min(str_toint(dom_fin,tmp2_S),F_istat)
   endif
   if (any((/F_dom_deb,dom_fin/) < 0) .or. .not.RMN_IS_OK(F_istat)) then
      call msg(MSG_ERROR,'(drv_ptopo_ndoms) Env variable UM_EXEC_NDOMAINS incorrectly defined: '//trim(ndomains_S))
      F_istat = RMN_ERR
      return
   endif

   F_ndomains = dom_fin - F_dom_deb + 1
   if (F_ndomains < 1) then
      write(tmp1_S,'(i4)') F_ndomains
      call msg(MSG_ERROR,'(drv_ptopo_ndoms)  ndomains_S= '//trim(tmp1_S)//' Check Env variable UM_EXEC_NDOMAINS: '//trim(ndomains_S))
      F_istat = RMN_ERR
      return
   endif

   F_istat = wb_put('ptopo_cfgs/ndomains',F_ndomains)
   F_istat = 0
   call msg(MSG_DEBUG,'[END] drv_ptopo_ndoms')
   !-------------------------------------------------------------------
   return
end subroutine drv_ptopo_ndoms

