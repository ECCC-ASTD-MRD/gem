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
! Fichier/File   : chm_snapshot.ftn90
! Creation       : A. Akingunola (Spring 2021)
! Description    : Save and restore snapshots of chemistry parameters at
!                  specified breakpoints
!
!
!==============================================================================
!!if_on
function chm_snapshot(F_mode) result(F_istat)
!!if_off
   use chm_utils_mod,  only: chm_lun_out, chm_rstn_l
   use chm_nml_mod,    only: chm_master, chm_cffeps_online_l
   use mach_cffeps_mod ! Inherits rmn_gmm

   implicit none

!!if_on
   character(len=*), intent(in) :: F_mode
   integer(kind=4) ::  F_istat
!!if_off
!
! Local variables
   integer(kind=4) :: istat, flag_r_n
   logical(kind=4), save :: do_once = .true.
   real(kind=4), pointer, dimension(:,:) :: lfire_emis_rstr, lfire_info_rstr
!
   F_istat = -1
!
!  Detect Master switch. If false, NORMAL EXIT WITH MESSAGE
!
   if (.not. chm_master) then
      if (chm_lun_out > 0) write(chm_lun_out, *) &
             'CHM_SNAPSHOT -> DETECTED CHEMICAL MASTER KILL'
      F_istat = 0
      return
   end if
!
!  Enable reading/writing snapshots (for instance, in case of digital filtering), and reading restarts
!  for CFFEPS fire emissions and hotspots state info.
   if (chm_cffeps_online_l .and. lc_hotspots > 0) then

      nullify(lfire_emis_rstr, lfire_info_rstr)

      if (trim(F_mode) == 'W') then
         if (do_once) then
            flag_r_n = GMM_FLAG_RSTR + GMM_FLAG_INAN
            istat = gmm_create('CFFEPS_INFO_rstr', lfire_info_rstr, meta_fire_info, flag_r_n)
            istat = gmm_create('CFFEPS_EMIS_rstr', lfire_emis_rstr, meta_fire_emis, flag_r_n)
            do_once = .false.
         end if
         istat = gmm_get('CFFEPS_INFO_rstr', lfire_info_rstr, meta_fire_info)
         istat = gmm_get('CFFEPS_EMIS_rstr', lfire_emis_rstr, meta_fire_emis)

         if (associated(lfire_info_rstr) .and. associated(lfire_emis_rstr)) then
            lfire_info_rstr = lfire_info
            lfire_emis_rstr = lfire_emissions
         else
            write(*, *)'(chm_snapshot) Cannot save CFFEPS FIRE INFO & EMISSIONS'
         end if
!
      else if (trim(F_mode) == 'R') then
         chm_rstn_l = .true.
         istat = gmm_get('CFFEPS_INFO_rstr', lfire_info_rstr, meta_fire_info)
         istat = gmm_get('CFFEPS_EMIS_rstr', lfire_emis_rstr, meta_fire_emis)
         if (associated(lfire_info_rstr) .and. associated(lfire_emis_rstr)) then
            lfire_info = lfire_info_rstr
            lfire_emissions = lfire_emis_rstr
         else
            write(*, *)'(chm_snapshot) Cannot restore CFFEPS FIRE INFO & EMISSIONS'
         end if
!
      end if
   end if
!
   F_istat = 1

   return
 end function chm_snapshot
