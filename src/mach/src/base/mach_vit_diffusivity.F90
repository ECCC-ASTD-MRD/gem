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
! Fichier/File   : mach_vit_diffusivity.ftn90
! Creation       : P. A Makar, A. Akingunola, Feb 2019
! Description    : Calculate vehicle-induced turbulence diffusivity.
!
! Extra info     :
!  References underlying the TKE formulae:
!  Miller et al., A Study of the Spatial Variation of Vehicle-Induced
!     Turbulence on Highways Using Measurements from a Mobile Platform,
!     Boundary Layer Meteorology, https://doi.org/10.1007/s10546-018-0416-9, 2018.
!  Gordon et al., Measurements of Enhanced Turbulent Mixing near Highways,
!     Journal of Applied Meteorology, 51, 1618-1632, 2012.
!--------------------------------------------------------------
!
! Arguments:
!           IN
!              dxdy   -> Averaged model grid area
!              zmom   -> Momentum level model heights
!              busper -> Permanent bus
!              busvol -> Volatile bus
!
!           OUT
!              kvkt   -> Vehicle Induced Turbulence diffusivity
!=============================================================================
!
!!if_on
subroutine mach_vit_diffusivity(busper, dxdy, zmom, kvkt, nmod, ni)
   use chm_ptopo_grid_mod,   only: chm_nk
!!if_off
   use chm_species_info_mod, only: sm
   use chm_species_idx_mod,  only: sp_CVKT, sp_MVKT, sp_TVKT
   implicit none
!
!!if_on
   integer(kind=4), intent (in) :: ni
   integer(kind=4), intent (in) :: nmod(ni)
   real(kind=4),    dimension(:), pointer, contiguous :: busper
   real(kind=4),    intent (in) :: dxdy(ni)
   real(kind=4),    intent (in) :: zmom(ni, chm_nk)
   real(kind=4),    intent(out) :: kvkt(ni, chm_nk)
!!if_off
!
! local variables:
!  Length scales for turbulence for cars (lc), mid-sized vehicles (lm) and trucks (lt),
!  based on the data of Gordon et al. (2012) and Miller et al. (2018):
!
   real(kind=4), parameter     :: lc = 13.859804
   real(kind=4), parameter     :: lm = 6.254106
   real(kind=4), parameter     :: lt = 11.288571
!
!  interval for estimation of vkt in the vertical, m:
   integer(kind=4), parameter  :: iznt = 50
!
   integer(kind=4)             :: ii, im, kk, k
   real(kind=4), dimension(ni) :: Fc, Fm, Ft, a
   real(kind=4)                :: odx, ccar, cmid, ctru, z, dz, sumF
!
! (1)
!  Extract the vehicle-km-travelled (km/s) values for car, mid-sized vehicle,
!  and truck for the current hour. DX is the dimension of the grid cell in km.
!
!  Fc, Fm, Ft are the number of cars, mid-sized, and trucks per second,
!  assuming the VKT values may be divided by the grid cell size to get the
!  number of vehicles going past a single point in the grid cell.
   do ii = 1, ni
      im = nmod(ii) - 1
      odx = 1.0 / (sqrt(dxdy(ii)) * 1.0e-3)
      Fc(ii) = busper(sm(sp_CVKT) % per_offset + im) * odx
      Fm(ii) = busper(sm(sp_MVKT) % per_offset + im) * odx
      Ft(ii) = busper(sm(sp_TVKT) % per_offset + im) * odx
   end do
!
!-------
!  (2) Divide each of the first three model levels up from the surface into 50 intervals,
!  calculate the VIT KT values at those levels, and average them to get the KT enhancement in that layer.
   kvkt = 0.0
   a = 0.0
! Lowest layer:
   do ii = 1, ni
      dz = zmom(ii, chm_nk) / real(iznt)
      sumF = Fc(ii) + Fm(ii) + Ft(ii)
      if (sumF > 1.E-22) then
         a(ii) = 0.4* (lc * Fc(ii) + lm * Fm(ii) + lt * Ft(ii)) / sumF
         do kk = 1, iznt
            z = real(kk - 1) * dz
! Following lines make the exponential term have a floor of 1E-22, to prevent
! underflow errors.  Term within the sqrt is the VIT TKE.
            ccar = max(-50.656, -2.39735E-02 * (z - 1.5 ) * (z - 1.5 ))
            cmid = max(-50.656, -1.17738E-01 * (z - 1.9 ) * (z - 1.9 ))
            ctru = max(-50.656, -3.61383E-02 * (z - 4.11) * (z - 4.11))
            kvkt(ii, chm_nk) = kvkt(ii, chm_nk) + a(ii) *          &
                               sqrt(2.42895E+00 * Fc(ii) * exp(ccar) + &
                                    1.55775E+01 * Fm(ii) * exp(cmid) + &
                                    2.04285E+01 * Ft(ii) * exp(ctru))
         end do  !kk
         kvkt(ii, chm_nk) = kvkt(ii, chm_nk) / real(iznt)
      end if
   end do ! ii

! Next two layers up:
   do k = chm_nk-1, chm_nk-2, -1
      do ii = 1, ni
         if (a(ii) > 0.0) then
            dz = (zmom(ii, k) - zmom(ii, k+1)) / real(iznt)
            do kk = 1, iznt
               z = real(kk - 1) * dz + zmom(ii, k)
               ccar = max(-50.656, -2.39735E-02 * (z - 1.5 ) * (z - 1.5 ))
               cmid = max(-50.656, -1.17738E-01 * (z - 1.9 ) * (z - 1.9 ))
               ctru = max(-50.656, -3.61383E-02 * (z - 4.11) * (z - 4.11))
               kvkt(ii, k) = kvkt(ii, k) + a(ii) *                   &
                             sqrt(2.42895E+00 * Fc(ii) * exp(ccar) + &
                                  1.55775E+01 * Fm(ii) * exp(cmid) + &
                                  2.04285E+01 * Ft(ii) * exp(ctru))
            end do  !kk
            kvkt(ii, k) = kvkt(ii, k) / real(iznt)
         end if
      end do ! ii
   end do  ! k
!
   return
end subroutine mach_vit_diffusivity
