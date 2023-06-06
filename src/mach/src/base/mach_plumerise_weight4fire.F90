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
!                            http: /  / www.ec.gc.ca                             !
!============================================================================!
!
! Projet / Project : GEM-MACH
! Fichier / File   : mach_plumerise_weight4fire.ftn90
! Creation       : Stephane Gaudreault, Alexander Kallaur, mars 2008 adaptation for GEM-MACH
!                  Wanmin Gong, for AURAMS 2004
!                  Original version from Janusz Pudykiewicz for CHRONOS 1995
!
! Description    : Assembly of major-point-source emission data for calculation
!                  with Eulerian transport model simulating atmospheric
!                  oxidants. Plume rise formula based on Briggs (1984).
!
! Arguments:
!            IN
!                cur_source          -> ID of the current source
!                z_magl              -> height of the model levels (m above ground)
!                boundlayer_height   -> eight of the boundary layer (m above ground)
!                pbl_indx            -> Index above boundary layer height
!                safe_inv_mo_length  -> Inverse of Monin obukhov_length
!                F_bio               -> Dominant Biomass vegetation fractions
!
!            IN / OUT
!                weight              -> Weight of each cell in the column over the point source
!
! Extra info     :
!            2019 Feb.  add a case wf_case=4 to use CFFEPS input
!
!============================================================================
!
!!if_on
subroutine mach_plumerise_weight4fire(cur_source, z_magl, pbl_indx, pbl_hgt, &
                                      safe_inv_mo_length, weight, F_bio, rho)
   use chm_ptopo_grid_mod,      only: chm_nk
!!if_off
   use chm_nml_mod,             only: wf_case
   use chm_utils_mod,           only: global_debug, chm_lun_out
   use chm_mjrpts_sortinfo_mod, only: lstack_info, i_gilc, i_vel, i_hgt

   implicit none
!!if_on
   integer(kind=4), intent (in) :: cur_source
   real(kind=4),    intent (in) :: z_magl(chm_nk+1)
   real(kind=4),    intent(out) :: weight(chm_nk)
   real(kind=4),    intent (in) :: pbl_hgt
   real(kind=4),    intent (in) :: F_bio (4)
   integer(kind=4), intent (in) :: pbl_indx
   real(kind=4),    intent (in) :: safe_inv_mo_length
   real(kind=4),    intent (in) :: rho   (chm_nk)
!!if_off
!
!  Local variables
!
   integer(kind=4) :: k                       ! loop index
   integer(kind=4) :: biome_index
   integer(kind=4) :: pbl_nlay
   integer(kind=4) :: zplm_indx, zplm_nlay
   real(kind=4)    :: stack_zplm, stack_rsmk
   integer(kind=4) :: unstable
   real(kind=4)    :: column_sum, column_sum_1
   real(kind=4)    :: Gaus_f
   real(kind=4)    :: increment
   real(kind=4)    :: stdev                   ! standard deviation used in the gaussian case
   real(kind=4)    :: stack_magl              ! Stack height in meter above ground
   real(kind=4)    :: stack_magl_f
                         ! based on Freitas et al 2007: ACP, 7, 3385-3398
   real(kind=4), dimension(4), parameter :: f_smolder = (/.55, .55, .25, .03/)
                         ! based on Val Martin et al 2010: ACP, 10, 1491-1510
   real(kind=4), dimension(4), parameter :: pl_med_hgt_m = (/1040., 781., 935., 691./)
   real(kind=4), dimension(4), parameter :: pl_med_hgt_std = (/646., 544., 604., 408./)
   real(kind=4), dimension(4), parameter :: pl_depth_m = (/1100., 828., 967., 625./)
   real(kind=4), dimension(4), parameter :: pl_depth_std = (/703., 653., 695., 468./)
   real(kind=4), dimension(chm_nk)       :: weight_s, weight_f
   real(kind=4)    :: weight_f_sum

   logical(kind=4) :: local_dbg

   local_dbg = (.false. .or. global_debug)

   weight = 0.
   weight_s = 0.
   weight_f = 0.
   unstable = 0
!
   if (wf_case < 3) then  ! Gaussian distribution within PBL
!
!  wf_case = 1 or 2   based on landuse
!
      biome_index = 4
      do k = 1,3
         if (F_bio(k) >= F_bio(biome_index)) biome_index = k
      end do
!
      stack_magl = pl_med_hgt_m(biome_index)

! determine final plume height and spread based on stability conditions

      if (stack_magl < pbl_hgt) then
         if (safe_inv_mo_length < 0.0 .and. (-stack_magl * safe_inv_mo_length) > 4.0) then
! Unstable condition
            stack_magl_f = stack_magl + pl_med_hgt_std(biome_index)
            stdev = 0.5 * (pl_depth_m(biome_index) + pl_depth_std(biome_index))
            unstable = 1
         else if (safe_inv_mo_length >= 0.0 .and. (stack_magl * safe_inv_mo_length) > 0.5) then
! Stable condition
            stack_magl_f = stack_magl - pl_med_hgt_std(biome_index)
            stdev = 0.5 * (pl_depth_m(biome_index) - pl_depth_std(biome_index))
            stdev = min(stdev, (pbl_hgt - stack_magl_f))
         else
! Neutral condition
            stack_magl_f = stack_magl
            stdev = 0.5 * pl_depth_m(biome_index)
         end if
      else
! Outside the boundary layer_height
         stack_magl_f = stack_magl
         stdev = 0.5 * (pl_depth_m(biome_index) - pl_depth_std(biome_index))
      end if
!
! Distributing flaming portion of the fire emission using Gaussian distribution
!
      column_sum = 0.0
      column_sum_1 = 0.0
      do k = 1, chm_nk-1
         increment = exp( -(z_magl(k) - stack_magl)**2 / (2.0 * stdev * stdev)) * (z_magl(k) - z_magl(k+1))
         column_sum = column_sum + increment
         if (k < pbl_indx) &
            column_sum_1 = column_sum_1 + increment
      end do
      column_sum = column_sum + exp( -(z_magl(chm_nk) - stack_magl)**2 / (2.0 * stdev * stdev)) * z_magl(chm_nk) 

      do k = 1, chm_nk
         weight_f(k) = exp(-(z_magl(k) - stack_magl)**2 / (2.0 * stdev * stdev)) / column_sum
      end do
!
   end if 

! evenly distributing smoldering portion of the fire emission within PBL
   if (wf_case == 3) then
      pbl_nlay = chm_nk - pbl_indx + 1  ! find # of layers below pbl
      do k = pbl_indx, chm_nk
         weight_s( k ) = 1.0 / pbl_nlay
      end do
   end if

!! CFFEPS case for distribution below ZPLM
   if (wf_case == 4) then
      stack_zplm = lstack_info(i_hgt, cur_source)
      stack_rsmk = lstack_info(i_vel, cur_source)
! Evaluate position zplm from the ground up
      do k = chm_nk, 1, -1
         zplm_indx = k
         if (z_magl(k) >= stack_zplm) then
            exit
         end if
      end do
      zplm_nlay = chm_nk - zplm_indx + 1  ! find # of layers below ZPLM
! calculate smoke density below zplm
      weight_f = 0.0
      weight_f(zplm_indx:chm_nk) = rho(zplm_indx:chm_nk) * stack_rsmk
      weight_f_sum = sum(weight_f)
      do k = zplm_indx, chm_nk
         if (weight_f_sum > 0.0) then
            weight_s(k) = weight_f(k) / weight_f_sum
         else
            weight_s(k) = 1.0 / zplm_nlay
         end if
      end do
! debug statements
      if (local_dbg) then
         write(chm_lun_out, *) 'wf_case: ',wf_case
         if (weight_f_sum <= 0.0) &
            write(chm_lun_out, *) 'WARNING:sum(weight)=0 ', stack_rsmk, stack_zplm
         write(chm_lun_out, *) 'z_magl: stk:', z_magl
         write(chm_lun_out, 111) 'weight_wf: stk:', weight_s
         write(chm_lun_out, 111) 'rho: stk:', rho
      end if
111 format(a, 100(f8.3))
   end if

! normalize weighing factor by layer thickness for subsequent S/R
   do k = 1,chm_nk-1
      weight_s(k) = weight_s(k) / (z_magl(k) - z_magl(k+1))
   end do
   weight_s(chm_nk) = weight_s(chm_nk) / z_magl(chm_nk)

!
!  Combining flaming and smoldering portions
!
   if (wf_case == 2 .and. unstable == 1)  then
      do k = 1, chm_nk
         weight(k) = weight_s(k) * f_smolder(biome_index) + &
                     (weight_s(k) * (1.0 - Gaus_f) + weight_f(k) * Gaus_f) * &
                     (1.0 - f_smolder(biome_index))
      end do
   else if (wf_case == 3 .or. wf_case == 4)  then
      do k = 1, chm_nk
         weight(k) = weight_s(k)
      end do
   else  ! either wf_case == 1 or wf_case == 2 and unstable /= 1
      do k = 1, chm_nk
         weight(k) = weight_s(k) * f_smolder(biome_index) + &
                     weight_f(k) * (1.0 - f_smolder(biome_index))
      end do
   end if

   return
end subroutine mach_plumerise_weight4fire
