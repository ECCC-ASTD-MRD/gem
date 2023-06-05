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
!               vdg ->                      Deposition velocity for ammonia through ground (m/s)
!
!               metvar2d(:, MV2D_TSURF) ->  Surface temperature (K)
!
!               metvar2d(:, MV2D_DXDY) ->   Grid cell area (m^2)
!
!               busper ->                   Permanent bus
!
!             IN/OUT
!
!               busvol ->                   Volatile bus
!
!============================================================================!
!
!!if_on
subroutine mach_gas_bidi(vdg, metvar2d, busper, busvol)
   use chm_metvar_mod,       only: SIZE_MV2D
   use mach_drydep_mod,      only: lucprm
   use chm_ptopo_grid_mod,   only: chm_ni
!!if_off
   use chm_metvar_mod,       only: MV2D_TSURF, MV2D_DXDY
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_NH3
   use chm_nml_mod,          only: chm_ammonia_bidi_s, chm_ammonia_gep

   implicit none
!!if_on
   real(kind=4),    intent   (in) :: vdg      (lucprm, chm_ni)
   real(kind=4),    intent   (in) :: metvar2d (chm_ni, SIZE_MV2D)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    dimension(:), pointer, contiguous :: busvol
!!if_off

   ! Local variables

   integer(kind=4) :: i, nlus, il, gopt
   real(kind=4)    :: tsurf, xg, mwt_nh3, dxdy, gamma, emis, coeff

   ! coefficients for compensation point calculation for ammonia (see Nemitz et al. 2000, doi.org/10.1016/S0168-1923(00)00206-9)
   real(kind=4), parameter :: A_CP_NH3 = 1.615E8  ! K*mol/m^3
   real(kind=4), parameter :: B_CP_NH3 = 10380.   ! K

   integer(kind=4), parameter :: GOPT_GEP    = 1
   integer(kind=4), parameter :: GOPT_GEP2D  = 2

   select case(trim(chm_ammonia_bidi_s))
   case('GEP')
      gopt = GOPT_GEP
   case('GEP2D')
      gopt = GOPT_GEP2D
   end select

   mwt_nh3 = sm(sp_NH3) % mol_wt

   ! loop over domain

   do i = 1, chm_ni

      tsurf = metvar2d(i, MV2D_TSURF)  ! K
      dxdy = metvar2d(i, MV2D_DXDY)    ! m^2

      coeff = mwt_nh3 * (A_CP_NH3/tsurf) * exp(-B_CP_NH3/tsurf)  ! g/m^3

      emis = 0.0

      ! loop over land-use categories

      do nlus = 1, lucprm

         select case(gopt)
         case(GOPT_GEP)
            gamma = chm_ammonia_gep(nlus)
         case(GOPT_GEP2D)
            il = (nlus - 1) * chm_ni + i - 1
            gamma = busper(sm(sp_NH3) % gep_offset + il)
         end select

         xg = gamma * coeff  ! g/m^3

         emis = emis + vdg(nlus, i) * xg * dxdy  ! g/s

      end do

      busvol(sm(sp_NH3) % bd_offset + i - 1) = emis

   end do

   return

end subroutine mach_gas_bidi
