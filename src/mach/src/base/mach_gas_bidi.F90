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
! Projet / Project : GEM-MACH
! Fichier / File   : mach_gas_bidi.ftn90
! Creation         : M. Sitwell - Jan 2022
! Description      : Calculates emissions from bidirectional flux for ammonia.
!
! Extra Info       : Current implementation allows for constant emissions potential
!                    values that are (1) a function of land-use category or (2) on
!                    the 2D horizontal grid input from FST files.
!
!                    This subroutine assumes that 15 land-use categories used in the
!                    ROBICHAUD' dry deposition scheme is used for the bidirectional
!                    flux as well.
!
! Arguments:  IN
!
!               chem_tr ->                  Chemical tracers concentrations (ug/kg)
!               metvar2d(:, MV2D_TSURF) ->  Surface temperature (K)
!               metvar2d(:, MV2D_DXDY)  ->  Grid cell area (m^2)
!               metvar3d(:, :, MV3D_RHO) -> Air density (kg/m^3)
!               vdg ->                      Deposition velocity for ammonia through ground (m/s)
!
!             IN/OUT
!
!               busper ->                   Permanent bus
!
!               busvol ->                   Volatile bus
!
!============================================================================!
!
!!if_on
subroutine mach_gas_bidi(busper, busvol, chem_tr, metvar2d, metvar3d, vdg)
   use mach_drydep_mod,      only: lucprm
   use chm_ptopo_grid_mod,   only: chm_ni, chm_nk
   use chm_species_info_mod, only: nb_dyn_tracers
   use chm_metvar_mod,       only: SIZE_MV2D, SIZE_MV3D
!!if_off
   use chm_metvar_mod,       only: MV2D_TSURF, MV2D_DXDY, MV3D_RHO, MV2D_WSOIL
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_NH3
   use chm_nml_mod,          only: chm_nh3_bidi_s, chm_nh3_gep2d_l, chm_nh3_gep, &
                                   chm_nh4_soil_loss, chm_soil_ph2d_l, chm_soil_ph
   use chm_utils_mod,        only: chm_timestep

   implicit none
!!if_on
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
   real(kind=4),    intent   (in) :: chem_tr  (chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4),    intent   (in) :: metvar2d (chm_ni, SIZE_MV2D)
   real(kind=4),    intent   (in) :: metvar3d (chm_ni, chm_nk, SIZE_MV3D)
   real(kind=4),    intent   (in) :: vdg      (lucprm, chm_ni)
!!if_off

   ! Local variables

   integer(kind=4) :: i, nlus, bus_index
   real(kind=4)    :: tsurf, xg, xa, mwt_nh3, dxdy, wsoil, rho, gamma, gamma_s
   real(kind=4)    :: emis, coeff, tau_soil, pH, lambda_a, tau_a
   logical(kind=4) :: gepdyn

   ! coefficients for compensation point calculation for ammonia (see Nemitz et al. 2000, doi.org/10.1016/S0168-1923(00)00206-9)
   real(kind=4), parameter :: A_CP_NH3 = 1.615E8  ! K*mol/m^3
   real(kind=4), parameter :: B_CP_NH3 = 10380.   ! K

   real(kind=4), parameter :: dsoil = 0.02  ! soil layer depth (m)

   gepdyn = trim(chm_nh3_bidi_s) == 'DYNAMIC'

   mwt_nh3 = sm(sp_NH3) % mol_wt  ! g/mol

   tau_soil = chm_nh4_soil_loss * 3600.  ! loss time constant in soil (s)

   ! loop over domain

   do i = 1, chm_ni

      tsurf = metvar2d(i, MV2D_TSURF)  ! K
      dxdy = metvar2d(i, MV2D_DXDY)    ! m^2

      coeff = (A_CP_NH3 / tsurf) * exp(-B_CP_NH3 / tsurf)  ! mol/m^3

      if (gepdyn) then

         rho = metvar3d(i, chm_nk, MV3D_RHO)  ! kg/m^3, lowest model level
         wsoil = metvar2d(i, MV2D_WSOIL)      ! m^3 water / m^3 soil
         xa = chem_tr(i, chm_nk, sp_NH3)      ! ug/kg, lowest model level

         wsoil = max(wsoil, tiny(wsoil))      ! impose minimum to prevent divide by zero in lambda_a computation
         xa = 1.E-6 * xa * rho / mwt_nh3      ! mol/m^3

         busvol(sm(sp_NH3) % epa_offset + i - 1) = xa / coeff  ! diagnostic atmospheric deposition potential

      end if

      emis = 0.0

      ! loop over land-use categories

      do nlus = 1, lucprm

         if (chm_nh3_gep2d_l .or. gepdyn) then
            bus_index = (nlus - 1) * chm_ni + i - 1
         end if

         ! get emissions potential value

         if (chm_nh3_gep2d_l) then
            gamma_s = busper(sm(sp_NH3) % gep_offset + bus_index)
         else
            gamma_s = chm_nh3_gep(nlus)
         end if

         if (gepdyn) then
            gamma = busper(sm(sp_NH3) % epd_offset + bus_index)
         else
            gamma = gamma_s
         end if

         ! calculate ground concentration and emissions

         xg = coeff * gamma  ! mol/m^3

         emis = emis + mwt_nh3 * vdg(nlus, i) * xg * dxdy  ! g/s

         ! update dynamic emissions potential if used

         if (gepdyn) then

            if (chm_soil_ph2d_l) then
               pH = busper(sm(sp_NH3) % sph_offset + bus_index)
            else
               pH = chm_soil_ph(nlus)
            end if

            lambda_a = vdg(nlus, i) * 10**(pH-3.) / (dsoil * wsoil)  ! 1/s * m^3/mol

            busper(sm(sp_NH3) % epd_offset + bus_index) = gamma + chm_timestep * ( lambda_a * (xa - xg) + (gamma_s - gamma) / tau_soil )

            ! compute bidirectional flux time scale for diagnostics
            if (lambda_a > tiny(lambda_a)) then
               tau_a = 1. / (lambda_a * coeff * 3600.)  ! hours
            else
               tau_a = huge(lambda_a)
            end if
            busvol(sm(sp_NH3) % bdt_offset + bus_index) = tau_a

         end if

      end do

      busvol(sm(sp_NH3) % bd_offset + i - 1) = emis

   end do

   return

end subroutine mach_gas_bidi
