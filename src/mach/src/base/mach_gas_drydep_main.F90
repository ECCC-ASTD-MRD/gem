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
!              chem_tr     --> Chemical tracers concentrations (ug/kg)
!              metvar2d    --> Array of 2-D MET fields from the physics buses
!              metvar3d    --> Array of 3-D MET fields from the physics buses
!              lfu         --> Land-use fractions
!              iseasn      --> Assigned season descriptors
!
!            IN/OUT
!              pvars --> Phybus pointer
!
!==============================================================================
!
!!if_on
subroutine mach_gas_drydep_main(pvars, chem_tr, metvar2d, metvar3d, &
                                lfu, iseasn)
   use phymem,               only: phyvar
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
   use mach_drydep_mod,      only: lucprm
!!if_off
   use chm_utils_mod,        only: ik, chm_lun_out, chm_msg_debug
   use chm_nml_mod,          only: chm_timings_L, chm_gas_drydep_s, chm_nh3_bidi_s
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_LU15, sp_LAI, sp_NH3
   use mach_gas_headers_mod, only: mach_gas_drydep_solver, mach_gas_drydep_stat, &
                                   mach_gas_drydep_ra, mach_gas_drydep_ra2,      &
                                   mach_gas_drydep_solver2, mach_gas_bidi
   use mach_drydep_mod,      only: nb_gas_depo, gas_depo
   use timing_omp,           only: timing_start_omp, timing_stop_omp
   implicit none
!!if_on
   type(phyvar),    pointer, contiguous :: pvars(:)
   real(kind=4),    intent   (in) :: chem_tr (chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: lfu     (chm_ni, lucprm)
   integer(kind=4), intent   (in) :: iseasn  (chm_ni)
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
#include <rmn/msg.h>
   !-----------------------------------------------------------------
   call msg_toall(chm_msg_debug, 'mach_gas_drydep [BEGIN]')
   if (chm_timings_L) call timing_start_omp(310, 'mach_gas_drydep', 480)
   local_dbg = (.false. .and. (chm_lun_out>0))

!  Calculate aerodynamic resistance for each landuse category

   select case (chm_gas_drydep_s)
      case ('ROBICHAUD', 'ROBICHAUD3')
         call mach_gas_drydep_ra(aero_resist, iseasn, lfu, metvar2d)
      case ('ROBICHAUD2')
         call mach_gas_drydep_ra2(aero_resist, iseasn, lfu, metvar2d)
   end select

   do ii = 1, chm_ni
      lai_2d(ii) = pvars(sm(sp_LAI)%per_pvarid)%data(ii)
   end do

!  Compute dry deposition velocity for all chemical species of interest

   if (trim(chm_nh3_bidi_s) /= 'OFF' .or. sm(sp_NH3)%vdg_pvarid > 0) then

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

   if (local_dbg) then
      call mach_gas_drydep_stat(vd, aero_resist, diff_resist, surf_resist, &
                                lfu, metvar2d)
   end if

   if (trim(chm_nh3_bidi_s) /= 'OFF') then
      call mach_gas_bidi(pvars, chem_tr, metvar2d, metvar3d, vdg)
   end if

   do lk = 1, lucprm
      do ii = 1, chm_ni
         il = (lk - 1) * chm_ni + ii
         pvars(sm(sp_LU15)%ra_pvarid)%data(il) = aero_resist(ii, lk)
      end do
   end do
   do sp_index = 1, nb_gas_depo
      busid = gas_depo(sp_index)%sp_id
      do ii = 1, chm_ni
         pvars(sm(busid)%vd_pvarid)%data(ii) = vd(sp_index, ii)
      end do
      if (sm(busid)%rb_pvarid > 0) then
         do ii = 1, chm_ni
            pvars(sm(busid)%rb_pvarid)%data(ii) = diff_resist(sp_index, ii)
         end do
      end if
      if (sm(busid)%rc_pvarid > 0) then
         do lk = 1, lucprm
             do ii = 1, chm_ni
                il = (lk - 1) * chm_ni + ii
                pvars(sm(busid)%rc_pvarid)%data(il) = surf_resist(lk, sp_index, ii)
             end do
          end do
       end if
   end do
   if (sm(sp_NH3)%vdg_pvarid > 0) then
      do lk = 1, lucprm
         do ii = 1, chm_ni
            il = (lk - 1) * chm_ni + ii
            pvars(sm(sp_NH3)%vdg_pvarid)%data(il) = vdg(lk, ii)
         end do
      end do
   end if

   if (allocated(vdg)) deallocate(vdg)

   call msg_toall(chm_msg_debug, 'mach_gas_drydep [END]')
   if (chm_timings_L) call timing_stop_omp(310)
   !-----------------------------------------------------------------

   return
end subroutine mach_gas_drydep_main
