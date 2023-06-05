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
! Fichier/File   : mach_calc_diag.ftn90
! Creation       : I. Ivanova  May 2015
!                  A. Kallaur
! Description    : Calculates averages and accumulators of chemical tendencies
!                  and diagnostics (similarly to 'calcdiag' in physics)
! Revisions
! 001    I. Ivanova    (May 2015) - add gas-phase & PM total column amounts
! 002    A. Akingunola (Mar.2018) - Merge chm_diag_colum_L diagnostics from then
!                                   various gas phase chemistry main routines
!
!Arguments
!
!          - Input/Output -
! busper   Permanent bus for chemistry
! busvol   Volatile bus for chemistry
!
!          - Input -
! step            Flag for first chem. step in current
! slab_index      Slice number
! chem_tr         Species concentrations (ug/kg)
! metvar2d        2D meteorology variables
! metvar3d        3D meteorology variables
!
!!if_on
subroutine mach_calc_diag(busvol, busper, chem_tr, metvar2d, metvar3d, step)
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
!!if_off
   use chm_utils_mod,        only: ik, chm_timestep
   use chm_datime_mod,       only: secondsin1hour
   use chm_consphychm_mod,   only: consth, mwt_air
   use chm_nml_mod,          only: chm_diag_accum_L, chm_diag_colum_L, &
                                   chm_moyhr, chm_acchr, chm_strato_s
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_SO2,  sp_SO4,  sp_NO,  sp_NO2, sp_O3,  &
                                   sp_N2O5, sp_HNO3, sp_PAN, sp_CO, sp_HCHO, &
                                   sp_ISOP, sp_col_O3,                       &
                                   sp_col_hu, sp_col_no, sp_col_no2,         &
                                   sp_tcol_CO, sp_tcol_HNO3, sp_tcol_ISOP,   &
                                   sp_tcol_N2O5, sp_tcol_NO, sp_tcol_NO2,    &
                                   sp_tcol_O3, sp_tcol_PAN, sp_tcol_SO2,     &
                                   sp_tcol_SO4, sp_tcol_HCHO, sp_tend_o3
   use chm_metvar_mod,       only: MV2D_PPLUS, MV3D_SIGT, MV3D_SIGM, MV3D_HUPLUS
   use linoz_param,          only: p_linoz_tropo, hu_linoz_tropo
   implicit none
!
!  Declaration of subroutine arguments
!
!
! Declare Subroutine Arguments
!
!!if_on
   integer(kind=4), intent(in)    :: step
   real   (kind=4), dimension(:), pointer, contiguous :: busper
   real   (kind=4), dimension(:), pointer, contiguous :: busvol
   real   (kind=4), intent   (in) :: chem_tr(chm_ni, chm_nk+1, nb_dyn_tracers)
   real   (kind=4), intent   (in) :: metvar2d(chm_ni, SIZE_MV2D)
   real   (kind=4), intent   (in) :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
!!if_off
!
! Declare Local variables
!
   integer(kind=4) :: i, k, isp, this_ik, sp_indx
   integer(kind=4), parameter :: nspecies = 11
   integer(kind=4) :: sp_conc(nspecies), sp_col(3), sp_tcol(nspecies)
   real(kind=4)    :: moyhr, moyhr_steps, hr_step, moysti, acchr, acchr_steps
!
   real(kind=4)    :: p_top(chm_ni, chm_nk), p_bot(chm_ni, chm_nk)
   real(kind=4)    :: p2d(chm_ni, chm_nk), hu_vmr(chm_ni, chm_nk)
   real(kind=4)    :: var_over(chm_ni, chm_nk), var_vmr(chm_ni, chm_nk)
   real(kind=4)    :: sp_over(chm_ni, chm_nk, nspecies), hu_over(chm_ni, chm_nk)
!
! Declare external subroutine
!
   external ghg_xcol
!
   sp_conc = (/sp_O3,   sp_NO2,  sp_NO,   sp_PAN,  sp_HNO3, sp_N2O5, &
               sp_SO2,  sp_SO4,  sp_CO,   sp_HCHO, sp_ISOP/)
   sp_col  = (/sp_col_O3,    sp_col_NO2,   sp_col_NO/)
   sp_tcol = (/sp_tcol_O3,   sp_tcol_NO2,  sp_tcol_NO,   sp_tcol_PAN,  &
               sp_tcol_HNO3, sp_tcol_N2O5, sp_tcol_SO2,  sp_tcol_SO4,  &
               sp_tcol_CO,   sp_tcol_HCHO, sp_tcol_ISOP/)
!
   do i = 1, chm_ni
      do k = 1, chm_nk
         hu_vmr(i, k) = consth * metvar3d(i, k, MV3D_HUPLUS)  ! units ppmv
         p2d(i, k)    = metvar2d(i, MV2D_PPLUS) * metvar3d(i, k, MV3D_SIGT)
         ! p_top & p_bot: Pressure(Pa)  on momentum levels
         p_top(i, k) = metvar3d(i, k, MV3D_SIGM) * metvar2d(i, MV2D_PPLUS)
         if (k == chm_nk) then
            p_bot(i, k) = metvar2d(i, MV2D_PPLUS)
         else
            p_bot(i, k) = metvar3d(i, k+1, MV3D_SIGM) * metvar2d(i, MV2D_PPLUS)
         endif
      end do
   end do
!
! Overhead column water vapor
   call ghg_xcol(hu_vmr, p_bot, p_top, hu_over, chm_ni, chm_nk)
!
! Overhead column of O3 NO2 NO PAN HNO3 NO3 N2O5 SO2 SO4 CO HCHO ISOP
   do isp = 1, nspecies
      sp_indx = sp_conc(isp)
      do k = 1, chm_nk
         do i = 1, chm_ni
            var_vmr(i, k) = chem_tr(i, k, sp_indx) * mwt_air / &
                                (sm(sp_indx) % mol_wt) * 1.E-9  !mole/mole vmr
         end do
      end do

      call ghg_xcol(var_vmr, p_bot, p_top, var_over, chm_ni, chm_nk)
      sp_over(:, :, isp) = var_over
   end do
!
   if (chm_diag_colum_L) then
!  Output 3D total columns of HU, NO2, NO, and O3 1E+15 (Peta) molecules cm-2;
      do k = 1, chm_nk
         do i = 1, chm_ni
            this_ik = ik(i, k, chm_ni)
!  YHU: Column water vapor
            busvol(sm(sp_col_hu) % out_offset + this_ik) = hu_over(i, k)
!  Column O3, NO2, and NO 1E+15 (Peta) molecules cm-2 (YO3, YNO2, YNO);
            do isp = 1, 3
               sp_indx = sp_col(isp)
               busvol(sm(sp_indx) % out_offset + this_ik) = sp_over(i, k, isp)
               busvol(sm(sp_indx) % out_offset + this_ik) = sp_over(i, k, isp)
               busvol(sm(sp_indx) % out_offset + this_ik) = sp_over(i, k, isp)
            end do
         end do
      end do
!
!  Output 2D Total (Tropospheric) Column amounts 1E+15 (Peta) molecules cm-2;
! ZO3 ZNO2 ZNOX ZPAN ZHN3 ZN25 ZSO2 ZSO4 ZCO ZHCH ZISO
      do isp = 1, nspecies
         sp_indx = sp_tcol(isp)
         do i = 1, chm_ni
            busvol(sm(sp_indx) % out_offset + i - 1) = sp_over(i, chm_nk, isp)
         end do
      end do
!
   end if
!
   if (.not. chm_diag_accum_L) return
!
!****************************************************************
!     ACCUMULATIONS                                             *
!     -------------                                             *
!****************************************************************
!
!  compute the averaging interval (inverse), with all available data averaged when
!  moving through driver step=0
   hr_step = chm_timestep * float(step) / sngl(secondsin1hour) ! hrs

   if (chm_moyhr > 0) then
      moyhr = float(chm_moyhr)
      moyhr_steps = moyhr / chm_timestep * sngl(secondsin1hour)
      moysti = 1.0 / moyhr_steps
   end if
!
!****************************************************************
!     AVERAGES                                                  *
!     --------                                                  *
!****************************************************************
!
!     set averages to zero every chm_moyhr hours
   if (chm_moyhr > 0) then
!
      RESET_AVERAGES: if (step == 0 .or. mod(step-1, int(moyhr_steps)).eq.0) then
!
         do isp = 1, nspecies
            do i = 1, chm_ni
               busper(sm(sp_tcol(isp)) % per_offset + i - 1) = 0.0
            end do
         end do
!
         do k = 1, chm_nk
            do i = 1, chm_ni
               this_ik = ik(i, k, chm_ni)
!
               busper(sm(sp_col_O3)  % per_offset + this_ik) = 0.0
               busper(sm(sp_col_NO2) % per_offset + this_ik) = 0.0
               busper(sm(sp_col_NO)  % per_offset + this_ik) = 0.0
               busper(sm(sp_col_HU)  % per_offset + this_ik) = 0.0
               busper(sm(sp_tend_O3) % per_offset + this_ik) = 0.0
!
               do isp = 1, nspecies
                  busper(sm(sp_conc(isp)) % per_offset + this_ik) = 0.0
               end do
            end do
         end do
!
      end if RESET_AVERAGES
!
   end if
!
   if (chm_moyhr > 0 .and. step /= 0) then
!
!   2D (AZO3 AZN2 AZNO AZPN AZH3 AZ25 AZS2 AZS4 AZCO AZHC AZIS)
!
      do isp = 1, nspecies
         sp_indx = sp_tcol(isp)
         do i = 1, chm_ni
            busper(sm(sp_indx) % per_offset + i - 1) = sp_over(i, chm_nk, isp) + &
                                  busper(sm(sp_indx) % per_offset + i - 1)
         end do
      end do
!
      if (mod(step, int(moyhr_steps)) == 0) then
         do isp = 1, nspecies
            sp_indx = sp_tcol(isp)
            do i = 1, chm_ni
               busper(sm(sp_indx) % per_offset + i - 1) = &
                         busper(sm(sp_indx) % per_offset + i - 1) * moysti
            end do
         end do
      end if
!
!  3D
!
      do k = 1, chm_nk
         do i = 1, chm_ni
            this_ik = ik(i, k, chm_ni)
! AYHU (Peta molec. cm-2)
            busper(sm(sp_col_HU) % per_offset + this_ik) = hu_over(i, k) + &
                      busper(sm(sp_col_HU) % per_offset + this_ik)
!
! AYO3 (D.U.); AYN2 & AYNO (Peta molec. cm-2)
            do isp = 1, 3
               busper(sm(sp_col(isp)) % per_offset + this_ik) = &
                 sp_over(i, k, isp) + busper(sm(sp_col(isp)) % per_offset + this_ik)
            end do
!
! AGOZ  (ug /kg /sec)
            busper(sm(sp_tend_O3) % per_offset + this_ik) = &
                      busper(sm(sp_tend_O3) % per_offset + this_ik) + &
                      busvol(sm(sp_tend_O3) % out_offset + this_ik)
!
! AMO3 AMN2 AMNO AMPN AMH3 AM25 AMS2 AMS4 AMCO AMHC AMIS  (ug /kg)
            do isp = 1, nspecies
               sp_indx = sp_conc(isp)
               busper(sm(sp_indx) % per_offset + this_ik) = &
                            busper(sm(sp_indx) % per_offset + this_ik) + &
                            chem_tr(i, k, sp_indx)
            end do

         end do  !i = 1, chm_ni
      end do    !k = 1, chm_nk
!
      if (mod(step, int(moyhr_steps)) == 0) then
         do k = 1, chm_nk
            do i = 1, chm_ni
               this_ik = ik(i, k, chm_ni)
!
               busper(sm(sp_tend_O3) % per_offset + this_ik) = &
                         busper(sm(sp_tend_O3) % per_offset + this_ik) * moysti
               do isp = 1, nspecies
                  sp_indx = sp_conc(isp)
                  busper(sm(sp_indx) % per_offset + this_ik) = &
                         busper(sm(sp_indx) % per_offset + this_ik) * moysti
               end do
!
            end do  !i = 1, chm_ni
         end do    !k = 1, chm_nk
      end if
!
   end if
!
!*************************************************************
!  ACCUMULATORS                                              *
!*************************************************************
!
!  set accumulators to zero at the beginning and after every chm_acchr hours, and by
!  default (chm_acchr=0) as the model step goes through 0.

   if (chm_acchr > 0) then
      acchr = float(chm_acchr)
      acchr_steps = acchr / chm_timestep * sngl(secondsin1hour)
   end if

   RESET_ACCUMULATORS: if (step == 0. .or. (chm_acchr > 0 .and. &
                           mod(step-1, int(acchr_steps)) == 0) .or. &
                           (chm_acchr == 0 .and. step-1 == 0)) then

!      do i = 1, chm_ni

! 2D

!         busper(sm(sp_tcol_O3) % per_offset + i - 1)  = 0.0
!      end do

!      do k = 1, chm_nk
!         do i = 1, chm_ni
! 3D
!            this_ik = ik(i, k, chm_ni)
! AGO3
!            busper(sm(sp_tend_O3 ) % per_offset + this_ik) = 0.0
!
!         end do
!      end do
!
   end if RESET_ACCUMULATORS
!
!   if (step /= 0) then
!
!VDIR NODEP
!      do i = 1, chm_ni
!
!                        Accumulation of 2D chemical species

! AZO3 (E+15 Peta molec cm-2 per acchr)
!         busper(sm(sp_tcol_O3 ) % per_offset + i - 1 ) = &
!                   busper(sm(sp_tcol_O3 ) % per_offset + i - 1)  + &
!                   busvol(sm(sp_tcol_O3 ) %  out_offset + i - 1) * chm_timestep

!      end do
!
!                     Accumulation of 3D chemical species
!
!      do k = 1, chm_nk
!         do i = 1, chm_ni

!            this_ik = ik(i, k, chm_ni)

! AGO3 (ug /kg /sec per acchr)
!            busper(sm(sp_tend_O3)  % per_offset + this_ik) = &
!                   busper(sm(sp_tend_O3  ) % per_offset + this_ik) + &
!                   busvol(sm(sp_tend_O3  ) % out_offset + this_ik) * chm_timestep
!         end do
!      end do
!
!   end if
!
!
   return
end subroutine mach_calc_diag
