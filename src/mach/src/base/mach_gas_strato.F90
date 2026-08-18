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
! Fichier/File   : mach_gas_strato.ftn90
! Creation       : J. de Grandpre ; Feb 2013
! Description    : Call to Stratospheric Ozone Chemistry (LINOZ)
!
! Extra info     : Update to use GEM-Linoz subroutines (D. Akingunola, Jul. 2020)
!                  see also $rpnphy/src/base/linoz.F90
! Arguments:
!            IN
!               busper       -> Physics permanent bus
!
!               zhuplus      -> Specific humidity vmr
!               ztplus       -> Temperature (K)
!               sigt         -> Local sigma values on thermodynamic levels
!               zpplus       -> Surface Pressure (Pa)
!               o3           -> Ozone concentration  (mmr in ug/kg)
!               step         -> Flag for first chem. step in current
!
!           IN/OUT
!               lino3_new     -> Updated Stratospheric Linoz ozone at all levels (ug/kg)
!
!================================================================================================
!!if_on
subroutine mach_gas_strato(busper, hu_ppm, ztplus, sigt, o3, zpplus, lino3_new, &
                           gni, gnk, nmod)
!!if_off
   use chm_ptopo_grid_mod,   only: chm_ni
   use chm_utils_mod       , only: ik, chm_timestep, CHM_MSG_DEBUG
   use chm_consphychm_mod  , only: mwt_air, consth
   use chm_nml_mod,          only: chm_timings_L
   use chm_species_info_mod, only: sm
   use chm_phyvar_mod,       only: o3ce, o3s
   use chm_species_idx_mod , only: sp_O3,  lin_c2, lin_c4, lin_c5, lin_c6, &
                                   lin_c7
   use linoz_param,          only: p_linoz_tropo, p_linoz_meso, p_linoz_c4, &
                                   QepsO3
   implicit none
!
!  Declaration of subroutine arguments
!
!!if_on
   integer(kind=4), intent (in) :: gni, gnk
   integer(kind=4), intent (in) :: nmod     (gni)
   real   (kind=4), dimension(:), pointer, contiguous :: busper
   real   (kind=4), intent (in) :: zpplus   (gni)
   real   (kind=4), intent (in) :: o3       (gni, gnk)  !ug/kg
   real   (kind=4), intent (in) :: hu_ppm   (gni, gnk)  !ppmv
   real   (kind=4), intent (in) :: ztplus   (gni, gnk)  !K
   real   (kind=4), intent (in) :: sigt     (gni, gnk)  !
   real   (kind=4), intent(out) :: lino3_new(gni, gnk)  !ug /kg
!!if_off
!
!  Declaration of local variables.
!
   integer(kind=4) :: i, k, this_ik
   real   (kind=4), dimension(gni, gnk+1) :: shtj, ztplus1
   real   (kind=4), dimension(gni, gnk)   :: pbot, ptop
   real   (kind=4), dimension(gni, gnk)   :: zttce, zlin4, zlin5, zlin6, &
                                             zlin7, cnull
   real   (kind=4), dimension(gni, gnk)   :: o3_vmr, o3c_vmr, zo3col,    &
                                             zo3ccol, o3_tend
   real   (kind=4), dimension(gni, gnk)   :: zhuplus
   real   (kind=4) :: zo3ce, zo3s, ugkg2vmr
!
!  Declaration of external subroutines
!
   external msg_toall, timing_start_omp, vssqrt, timing_stop_omp, &
            linoz_xcol, linoz_tend
!
   !-----------------------------------------------------------------
   call msg_toall(CHM_MSG_DEBUG, 'mach_gas_strato [BEGIN]')
   if (chm_timings_L) call timing_start_omp(338, 'mach_gas_strato', 330)

   ugkg2vmr = mwt_air / sm(sp_O3) % mol_wt
!
   ztplus1(:, 1:gnk) = ztplus(:, 1:gnk)
   do i = 1, gni
      shtj(i,1) = sigt(i,1) * sqrt(sigt(i,1) / sigt(i,2))
      shtj(i,gnk) = 1.0
   end do

   do k = 2, gnk
      do i = 1, gni
         shtj(i,k) = sqrt(sigt(i,k-1) * sigt(i,k))
      end do
   end do

   do i = 1, gni
      shtj(i, gnk+1) = shtj(i, gnk)
      ztplus1(i, gnk+1) = ztplus(i, gnk)
   end do

   do k = 1, gnk
      do i = 1, gni
         this_ik = ik(nmod(i), k, chm_ni)
         zttce(i, k) = busper(sm(lin_c2) % per_offset + this_ik)
         zlin4(i, k) = busper(sm(lin_c4) % per_offset + this_ik)
         zlin5(i, k) = busper(sm(lin_c5) % per_offset + this_ik)
         zlin6(i, k) = busper(sm(lin_c6) % per_offset + this_ik)
         zlin7(i, k) = busper(sm(lin_c7) % per_offset + this_ik)

         ptop(i, k) = shtj(i, k)   * zpplus(i)  ! pressure (Pa) at the upper interface
         pbot(i, k) = shtj(i, k+1) * zpplus(i)  ! pressure (Pa) at the bottom interface

            ! Set lower limit on species as in BIRA (units mole /mole)
         o3_vmr(i, k) = o3(i, k) * 1.0E-9 * ugkg2vmr      !mole /mole vmr <-- micro g/kg air
         o3_vmr(i, k) = max(o3_vmr(i, k), QepsO3)

            ! FK-HALOE ozone climatology from cccmarad blended with ERA5 below 1 hPa
         zo3s  = busper(o3s  + this_ik)
         o3c_vmr(i, k) = zo3s * ugkg2vmr                  !mole /mole vmr  <-- kg /kg

         zo3ce = busper(o3ce + this_ik)
         if (zo3ce >= 0.0 .and. ptop(i, k) > p_linoz_meso) then
            o3c_vmr(i, k) = zo3ce * ugkg2vmr              !mole /mole vmr <-- kg /kg air
         end if

         ! Ignore P-L term on RHS above p_linoz_c4=10hPa
         if (ptop(i, k) < p_linoz_c4 ) zlin4(i, k) = 0.
      end do
   end do

   ! Total Column LINOZ ozone (D.U.)
   call linoz_xcol(o3_vmr, pbot, ptop, zo3col, gni, gnk)

      ! Total column ERA5-FK-HALOE ozone climatology (D.U.)
   call linoz_xcol(o3c_vmr, pbot, ptop, zo3ccol, gni, gnk)
!
      !    Linoz tendencies
   cnull = 0.0
   zhuplus = hu_ppm / consth
   call linoz_tend(                                         &
           o3_vmr                                         , & !input, mole/mole vmr
           zo3col                                         , & !input, D.U.
           ztplus1, zpplus, shtj, zhuplus                 , & !input
           o3c_vmr,zttce,zo3ccol,zlin4,zlin5,zlin6,zlin7  , & !input mole/mole vmr ERA-3 ozone climato in troposphere
           lino3_new                                      , & !output, mole /mole
           cnull ,cnull, cnull ,cnull                     , & !output, mole /mole /sec
           o3_tend                                        , & !output, mole /mole /sec
           chm_timestep, gni, gnk, gnk+1)                     !input

   lino3_new = lino3_new / ugkg2vmr * 1.0E9

   call msg_toall(CHM_MSG_DEBUG, 'mach_gas_strato [END]')
   if (chm_timings_L) call timing_stop_omp(338)
   !-----------------------------------------------------------------
   return
end subroutine mach_gas_strato
