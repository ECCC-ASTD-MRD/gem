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
! Fichier/File   : mach_incld_intrqf.ftn90
! Creation       : S. Menard, GEM-MACH, June 2008
! Description    : Compute intra-hydrometeor flux terms - this version only
!                  computes riming terms.  It is called only when f13 > 0.
!
! Extra info     : QA  ADOM  VERSION: ADOMIIB(CLEAN)  LEVEL: 06/09/89  INTRQF  ENSR(AES)
!
! Arguments  IN
!               ideriv -->
!               temp   --> Atmospheric Temperature (Kelvin)
!               f13    -->
!               aq     --> Aqueous sp. conc(molar) in cloud water & Ice/Snow
!               q      --> Cloudwater, ice/snow
!            OUT
!               PPAQ   --> For cloudwater
!               DDAQ   --> For ice/snow
!               CCDAQ  --> Change
!=============================================================================
!
!!if_on
subroutine mach_incld_intrqf (aq, ppaq, ddaq, ccdaq, q, f13, ideriv, temp)
   use mach_cam_utils_mod, only:maxnsaq
!!if_off
   implicit none
! !!if_ons
   integer(kind=4), intent (in) :: ideriv
   real(kind=4),    intent (in) :: temp
   real(kind=4),    intent (in) :: f13
   real(kind=4),    intent (in) :: q(2)
   real(kind=4),    intent (in) :: aq(maxnsaq, 2)
   real(kind=4),    intent(out) :: ppaq(maxnsaq, 2)
   real(kind=4),    intent(out) :: ddaq(maxnsaq, 2)
   real(kind=4),    intent(out) :: ccdaq(maxnsaq, 2)
!!if_off
!
! Local variables
!
   integer(kind=4) :: naqend, i
   real(kind=4)    :: delt, b13h, dem
!
!  aqueous species
!  HSO3 H2O2 ROOH SO4= NO3- NH4+ CAT1 HCO3- H+   OH- FEMN O3  H2OA
!  1    2    3    4    5    6    7     8    9   10   11  12   (constant = 1.0)
!
!  fraction of aq cloudwater species frozen in ice/snow
!  HSO3, H202, ROOH, SO4=, NO3-, NH4+, CAT1, HCO3-, H+, OH-, FEMN, O3
!
   real(kind=4), dimension(12) :: b13frz = (/.25, .25, .25, 1., .25, .25, 1., .25, 1., 1., 1., .25/)
   real(kind=4)                :: aa=5.8e-3, bb=1.2e-2, tzero=273.16

!  zero out
   naqend = maxnsaq
   if (ideriv < 0) naqend = 3
   do i = 1, naqend
      ppaq(i, 1) = 0.
      ddaq(i, 2) = 0.
   end do

!  compute fraction of so2 entrapped by rime ice using lamb and
!  blumenstein's equation (atmos. environ., vol. 21, no. 8, pp.
!  1765-1772, 1987).
   delt = max(tzero - temp, 0.)
   b13frz(1) = aa * delt + bb

!  compute fraction of protons to transfer in riming (f13)
      b13h = b13frz(9)
      dem = aq(4, 1) * 2. + aq(1, 1) + aq(5, 1) - aq(6, 1) + aq(8, 1)
      dem = max(dem, aq(9, 1))
      dem = max(dem, 1.e-12)
      b13frz(9) = (b13frz(4) * aq(4, 1) * 2. + b13frz(1) * aq(1, 1)  &
                 + b13frz(5) * aq(5, 1) - b13frz(6) * aq(6, 1)       &
                 + b13frz(8) * aq(8, 1)) / dem

!  cloudwater
   do i = 1, naqend
      ddaq(i, 1) = b13frz(i) * f13 / q(1)

!  ice/snow
      ppaq(i, 2) = b13frz(i) * aq(i, 1) * f13 / q(2)

!  compute cdot
      ccdaq(i, 1) = ppaq(i, 1) - ddaq(i, 1) * aq(i, 1)
      ccdaq(i, 2) = ppaq(i, 2) - ddaq(i, 2) * aq(i, 2)
   end do
   b13frz(9) = b13h

   return
end subroutine mach_incld_intrqf
