!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2020 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : chm_init.ftn90
! Creation       : A. Akingunola (Summer 2020)
! Description    : Initialization of the sundry chemistry parameters at the beginning
!                  of each execution of the model
!
!
!==============================================================================
!!if_on
function chm_init(F_path_S) result(F_istat)
!!if_off
   use chm_utils_mod,        only: chm_lun_out
   use chm_nml_mod,          only: chm_master, chm_pkg_gas_s
   use chm_headers_mod,      only: chm_set_dbg_point
   use mach_adom2_rates_mod, only: mach_adom2_chemi

   implicit none

#include <rmn/msg.h>

   include 'mach_version.inc'

!!if_on
   character(len=*), intent(in) :: F_path_S !# data/tables dir
   integer(kind=4) ::  F_istat
!!if_off
!
! Local variables
   integer(kind=4) :: istat
   character(len=1024) :: msg_S
!
   F_istat = -1

   !# Print Banner
   call msg(MSG_INFO, '   ************************************************************')
   write(msg_S, "(3x,'Package: ',a,5x,'version: ',a)") &
           trim(MACH_NAME_S), trim(MACH_VERSION_S)
   call msg(MSG_INFO, msg_S)
   write(msg_S, "(3x,'Release of: ',a,5x,'COMPILER: ',a)") &
           trim(MACH_DSTP_S), trim(MACH_EC_ARCH_S)
   call msg(MSG_INFO, msg_S)
   call msg(MSG_INFO, '   ************************************************************')

!
!  Detect Master switch. If false, NORMAL EXIT WITH MESSAGE
!
   if (.not. chm_master) then
      if (chm_lun_out > 0) write(chm_lun_out, *) &
             'CHM_INIT -> DETECTED CHEMICAL MASTER KILL'
      F_istat = 0
      return
   end if
!
!  Set debug grid point
   istat = chm_set_dbg_point()
   if (istat < 0) then
      write(chm_lun_out, *) '### Error in chm_set_dbg_point ###'
      return
   end if
!  get rpn surface bus indexes
   call chm_getphybus_struct()

   ! Initialize some values needed for ADOM2 Y&B solver
   if (chm_pkg_gas_s(1:5) == 'ADOM2') call mach_adom2_chemi()
!
   F_istat = 1

   return
end function chm_init
