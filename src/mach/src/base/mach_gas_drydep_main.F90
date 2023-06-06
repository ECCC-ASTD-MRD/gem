!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
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
! Fichier/File   : mach_gas_drydep_main.ftn90
! Creation       : A. Kallaur S. Menard, H. Landry, P.A. Beaulieu  - Jan 2008
! Description    : Prepare to compute dry deposition resistances and velocities
!                  for selected gas species.
!
! Arguments:
!            IN

!              busper      --> Permanent bus
!              metvar2d    --> Array of 2-D MET fields from the Physics buses
!              lfu         --> Land-use fractions
!              iseasn      --> Assigned season descriptors
!
!            IN/OUT
!              busvol --> Volatile bus
!
!==============================================================================
!
!!if_on
subroutine mach_gas_drydep_main(busper, busvol, metvar2d, lfu, iseasn)
   use chm_ptopo_grid_mod,   only: chm_ni
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_utils_mod,        only: global_debug, ik, CHM_MSG_DEBUG
   use chm_nml_mod,          only: chm_timings_L, chm_gas_drydep_s, chm_ammonia_bidi_s
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_LU15, sp_LAI, sp_NH3
   use mach_gas_headers_mod, only: mach_gas_drydep_solver, mach_gas_drydep_stat, &
                                   mach_gas_drydep_ra, mach_gas_drydep_ra2,      &
                                   mach_gas_drydep_solver2, mach_gas_bidi
   use mach_drydep_mod,      only: nb_gas_depo, gas_depo
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: iseasn  (chm_ni)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
!!if_off
!
!  Local variables
!
   integer(kind=4) :: ii, lk, il, sp_index, busid
   real(kind=4)    :: lai_2d(chm_ni)
   real(kind=4)    :: vd(nb_gas_depo, chm_ni), diff_resist(nb_gas_depo, chm_ni)
   real(kind=4)    :: aero_resist(chm_ni, lucprm)
   real(kind=4)    :: surf_resist(lucprm, nb_gas_depo, chm_ni)
   logical(kind=4) :: local_dbg

   real(kind=4), allocatable :: vdg(:,:)
!
!  External subroutines
!
   external msg_toall, timing_start_omp, timing_stop_omp
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_gas_drydep [BEGIN]')
   if (chm_timings_L) call timing_start_omp(310, 'mach_gas_drydep', 480)

!  Calculate aerodynamic resistance for each landuse category

   select case (chm_gas_drydep_s)
      case ('ROBICHAUD', 'ROBICHAUD3')
         call mach_gas_drydep_ra(aero_resist, iseasn, lfu, metvar2d)
      case ('ROBICHAUD2')
         call mach_gas_drydep_ra2(aero_resist, iseasn, lfu, metvar2d)
   end select

   do ii = 1, chm_ni
      lai_2d(ii) = busper(sm(sp_LAI) % per_offset + ii - 1)
   end do

!  Compute dry deposition velocity for all chemical species of interest

   if (trim(chm_ammonia_bidi_s) /= 'OFF' .or. sm(sp_NH3) % vdg_offset > 0) then

      allocate(vdg(lucprm, chm_ni))

      select case (chm_gas_drydep_s)
      case ('ROBICHAUD','ROBICHAUD3')
         call mach_gas_drydep_solver(vd, aero_resist, diff_resist, surf_resist, &
                                     iseasn, lfu, lai_2d, metvar2d, vdg=vdg)
      case ('ROBICHAUD2')
         call mach_gas_drydep_solver2(vd, aero_resist, diff_resist, surf_resist, &
                                      iseasn, lfu, metvar2d, vdg=vdg)
      end select

   else

      select case (chm_gas_drydep_s)
      case ('ROBICHAUD','ROBICHAUD3')
         call mach_gas_drydep_solver(vd, aero_resist, diff_resist, surf_resist, &
                                     iseasn, lfu, lai_2d, metvar2d)
      case ('ROBICHAUD2')
         call mach_gas_drydep_solver2(vd, aero_resist, diff_resist, surf_resist, &
                                      iseasn, lfu, metvar2d)
      end select

   end if

   local_dbg = .false.
   if (local_dbg) then
      call mach_gas_drydep_stat(vd, aero_resist, diff_resist, surf_resist, &
                                lfu, metvar2d)
   end if

   if (trim(chm_ammonia_bidi_s) /= 'OFF') then
      call mach_gas_bidi(vdg, metvar2d, busper, busvol)
   end if

   do lk = 1, lucprm
      do ii = 1, chm_ni
         il = (lk - 1) * chm_ni + (ii - 1)
         busvol(sm(sp_LU15) % ra_offset + il) = aero_resist(ii, lk)
      end do
   end do
   do sp_index = 1, nb_gas_depo
      busid = gas_depo(sp_index) % sp_id
      do ii = 1, chm_ni
         busvol(sm(busid) % vd_offset + ii - 1) = vd(sp_index, ii)
      end do
      if (sm(busid) % rb_offset > 0) then
         do ii = 1, chm_ni
            busvol(sm(busid) % rb_offset + ii - 1) = diff_resist(sp_index, ii)
         end do
      end if
      if (sm(busid) % rc_offset > 0) then
         do lk = 1, lucprm
             do ii = 1, chm_ni
                il = (lk - 1) * chm_ni + (ii - 1)
                busvol(sm(busid) % rc_offset + il) = surf_resist(lk, sp_index, ii)
             end do
          end do
       end if
   end do
   if (sm(sp_NH3) % vdg_offset > 0) then
      do lk = 1, lucprm
         do ii = 1, chm_ni
            il = (lk - 1) * chm_ni + (ii - 1)
            busvol(sm(sp_NH3) % vdg_offset + il) = vdg(lk, ii)
         end do
      end do
   end if

   if (allocated(vdg)) deallocate(vdg)

   call msg_toall(CHM_MSG_DEBUG, 'mach_gas_drydep [END]')
   if (chm_timings_L) call timing_stop_omp(310)
   !-----------------------------------------------------------------

   return
end subroutine mach_gas_drydep_main
