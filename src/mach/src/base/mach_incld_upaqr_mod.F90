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
! Fichier/File   : mach_incld_upaqr_mod.ftn90
! Creation       : S. Menard, S. Gravel, GEM-MACH, June 2008
! Description    : Calculations of equilibrium & kinetic reactions rates for q > 0
!
! Extra info     : ADOM  VERSION: ADOMIIB(CLEAN)  LEVEL: 06/09/89  STEADY  ENSR(AES)
!
! Arguments  IN
!             q(1)     --> cloudwater conc (g-w/m3 air)
!             q(2)     --> Ice/Snow conc (g-w/m3 air)
!             ti       --> Inverse temperature (1/K)
!             rtw      --> TEMP*RGAS/WDEN
!              AQ      --> Aqueous sp. conc(molar) in cloud water & Ice/Snow
!             pgscav   --> Gas/part scavenging (or diffusion)
!
!            OUT
!              R       --> Rate constants for aqueous-phase
!
!=============================================================================
!
!
!
!                 SPECIES LIST
!       AQUEOUS                    GAS/PART.
!       -
!       1.  HSO3                   1.  SO2G
!       2.  H202                   2.  HPXG
!       3.  ROOH                   3.  RPXG
!       4.  SO4=                   4.  H2SO4 = SO4P1
!       5.  NO3-                   5.  NH4HSO4 = SO4P2
!       6.  NH4+                   6.  (NH4)2SO4 = SO4P3
!       7.  CAT1                   7.  HNO3G
!       8.  HCO3-                  8.  NH3G
!       9.  H+                     9.  NH4NO3 = NNO3P
!      10.  OH-                   10.  DUST
!      11.  FEMN                  11.  O3G
!      12.  O3                    12.  CO2G (CONSTANT)
!           H2OA=1.0 (CONSTANT)
!
!                  REACTION LIST
!                  -
!
!  1.     1.SO4P1          -->    B1*SO4= +   B3*H+
!  2.     1.SO2G           -->    B1*HSO3 +   B1*H+
!  3.     1.HSO3 +  1.H+   -->    B2*SO2G
!  4.     1.O3G            -->    B1*O3
!  5.     1.O3             -->    B2*O3G
!  6.     1.HPXG           -->    B1*H2O2
!  7.     1.H2O2           -->    B2*HPXG
!  8.     1.HNO3           -->    B1*NO3- +   B1*H+
!  9.     1.NO3- +  1.H+   -->    B2*HNO3
! 10.     1.RPXG           -->    B1*ROOH
! 11.     1.ROOH           -->    B2*RPXG
! 12.     1.NH3G           -->    B1*NH4+ +   B1*OH-
! 13.     1.NH4+ +  1.OH-  -->    B2*NH3G
! 14.     1.DUST           -->    B5*FEMN +   B4*HCO3 + B4*CAT1+
! 15.     1.CO2G           -->    B1*HCO3 +   B1*H+
! 16.     1.HCO3 +  1.H+   -->    B2*CO2G
! 17.     1.H+   +  1.OH-  -->   1.00H2OA
! 18.     1.H2OA           -->   1.00H+   +  1.00OH-
! 19.     1.HSO3 +  1.O3   -->   1.00SO4= +  1.00H+
! 20.     1.HSO3 +  1.H2O2 -->   1.00SO4= +  1.00H+
! 21.     1.HSO3 +  1.ROOH -->   1.00SO4= +  1.00H+
! 22.     1.HSO3    (FEMN) -->   1.00SO4= +  1.00H+
! 23.     1.SO4P2          -->    B1*SO4= +   B1*H+ + B1*NH4+
! 24.     1.SO4P3          -->    B1*SO4= + B3*NH4+
! 25.     1.NNO3P          -->    B1*NO3- + B1*NH4+
!
 module mach_incld_upaqr_mod
   use mach_cam_utils_mod, only: maxnsaq, maxnsg
   use mach_incld_mod,     only: doh2o2, dorooh, domnfe
   private
   public :: mach_incld_upaqr
!
   interface mach_incld_upaqr
      module procedure mach_incld_upaqr1
      module procedure mach_incld_upaqr2
      module procedure mach_incld_upaqr3
   end interface mach_incld_upaqr

 contains

!  compute rate constants that depend on [c]'s (each time step)
   subroutine mach_incld_upaqr3(aq, ti, r)
      implicit none
      real(kind=4),    intent   (in) :: ti
      real(kind=4),    intent(inout) :: r (25, 2)
      real(kind=4),    intent   (in) :: aq(maxnsaq, 2)
!
!  Local variables
      real(kind=4)       :: xhp, xfe, xmn, xph, tfact, h17

!
!  evaluate special rate constants
      h17 = 10. ** (-3.83 - 3030. * ti)
      r(18, 1) = 0.1 * aq(9, 1)
      r(17, 1) = r(18, 1) / h17
!
!  hso3 + o3
      r(19, 1) = (9.560496e9 + 2641716. / aq(9, 1)) * exp(-3019.6 * ti)
!  hso3 + h2o2
      if (doh2o2) then
         r(20, 1) = (1.81507711e13 * aq(9, 1)) * exp(-3673.8 * ti)
      end if
!  hso3 + rooh
      if (dorooh) then
!  organic peroxide oxidation rate = 50% of lind et al. 1984
!  paa rate.  2/21/85
         r(21, 1) = exp(-3994.*(ti - (1. / 298.))) * (1.82e7 * aq(9, 1) + 300.3)
      end if
!  hso3 --> in presence of fe/mn
      if (domnfe) then
!  use ibusuki and takeuchis's expression for fe/mn catalyzed
!  oxidation rate (atmos. environ., vol. 21, no. 7, pp.1555-
!  1560, 1987).
         xhp = aq(9, 1)
         xfe = .95 * aq(11, 1)
         xmn = .05 * aq(11, 1)
         xph = -alog10(xhp)
!
!  note: 0.0033693 = 1/296.8k
         tfact = exp(-8431.16 * (ti - 0.0033693))
         if (xph < 4.2) then
            r(22, 1) = 3.72e7 * xmn * xfe * tfact * xhp ** (-0.74)
         else
            r(22, 1) = 2.51e13 * xmn * xfe * tfact * xhp ** 0.67
         end if
      end if
   
      return
   end subroutine mach_incld_upaqr3

!  A vectorized version of compute rate constants that depend on [c]'s (each time step)
   subroutine mach_incld_upaqr2(aq, ti, r, q, nptsnz)
      implicit none
      integer(kind=4), intent   (in) :: nptsnz
      real(kind=4),    intent   (in) :: aq(nptsnz, maxnsaq, 2)
      real(kind=4),    intent   (in) :: ti(nptsnz)
      real(kind=4),    intent   (in) :: q (nptsnz, 2)
      real(kind=4),    intent(inout) :: r (nptsnz, 25, 2)
!
!  Local variables
      integer(kind=4) :: ik
      real(kind=4)    :: xhp, xfe, xmn, xph, tfact, h17

!
      do ik = 1, nptsnz
         if (q(ik, 1) <= 0.0) cycle
!  evaluate special rate constants
         h17 = 10. ** (-3.83 - 3030. * ti(ik))
         r(ik, 18, 1) = 0.1 * aq(ik, 9, 1)
         r(ik, 17, 1) = r(ik, 18, 1) / h17
!
!  hso3 + o3
         r(ik, 19, 1) = (9.560496e9 + 2641716. / aq(ik, 9, 1)) * &
                        exp(-3019.6 * ti(ik))
!  hso3 + h2o2
         if (doh2o2) then
            r(ik, 20, 1) = (1.81507711e13 * aq(ik, 9, 1)) * exp(-3673.8 * ti(ik))
         end if
!  hso3 + rooh
         if (dorooh) then
!  organic peroxide oxidation rate = 50% of lind et al. 1984
!  paa rate.  2/21/85
            r(ik, 21, 1) = exp(-3994.*(ti(ik) - (1. / 298.))) * &
                           (1.82e7 * aq(ik, 9, 1) + 300.3)
         end if
!  hso3 --> in presence of fe/mn
         if (domnfe) then
!  use ibusuki and takeuchis's expression for fe/mn catalyzed
!  oxidation rate (atmos. environ., vol. 21, no. 7, pp.1555-
!  1560, 1987).
            xhp = aq(ik, 9, 1)
            xfe = .95 * aq(ik, 11, 1)
            xmn = .05 * aq(ik, 11, 1)
            xph = -alog10(xhp)
!
!  note: 0.0033693 = 1/296.8k
            tfact = exp(-8431.16 * (ti(ik) - 0.0033693))
            if (xph < 4.2) then
               r(ik, 22, 1) = 3.72e7 * xmn * xfe * tfact * xhp ** (-0.74)
            else
               r(ik, 22, 1) = 2.51e13 * xmn * xfe * tfact * xhp ** 0.67
            end if
         end if
      end do
   
      return
   end subroutine mach_incld_upaqr2

! Initial evaluation of rate constants
   subroutine mach_incld_upaqr1(aq, ti, r, q, nptsnz, pgscav, rtw)
      implicit none
      integer(kind=4), intent (in) :: nptsnz
      real(kind=4),    intent (in) :: aq    (nptsnz, maxnsaq, 2)
      real(kind=4),    intent (in) :: ti    (nptsnz)
      real(kind=4),    intent (in) :: q     (nptsnz, 2)
      real(kind=4),    intent (in) :: pgscav(nptsnz, maxnsg, 2)
      real(kind=4),    intent (in) :: rtw   (nptsnz)
      real(kind=4),    intent(out) :: r     (nptsnz, 25, 2)
!
! Local varariables
      integer(kind=4) :: ik
!
      real(kind=4)    :: lrt, tsqi
      real(kind=4)    :: h2a, h2b, h4, h6, h8, h10, h12a, h12b, h15a, h15b

!  initialize rates
      r = 0.0

!  evaluate necessary henry's law coefficients
      do ik = 1, nptsnz
         if (q(ik, 1) <= 0.0) cycle
         tsqi = ti(ik) * ti(ik)
         h2a = 10.0 ** (-10.65 + 1410.0 * ti(ik))
         h2b = 10.0 ** (-4.84 + 870.0 * ti(ik))
         h4 = 10.0 ** (-11.16 + 1040.0 * ti(ik))
         h6 = 10.0 ** (-11.15 + 2990.0 * ti(ik))
         h8 = 10.0 ** (-12.17 + 3780.0 * ti(ik))
!  use recent data (lind and kok, jgr, vol. 91, pp 7889-7895,
!  june 1986) for rooh solubility
         h10 = exp(5607 * ti(ik) - 13.41) * 1.0e-6
         h12a = 10.0 ** (-9.50 + 1600.0 * ti(ik))
         h12b = 10.0 ** (-13.914 + 2730.0 * ti(ik))
         h15a = 10.0 ** (-10.66 + 760.0 * ti(ik) + 58000.0 * tsqi)
         h15b = 10.0 ** (-14.25 + 5190.0 * ti(ik) - 850000.0 * tsqi)

!  note: pgscav inputs for part. by q(:,:)  + gases by q(:,2) are used!
!  evaluate mass transfer/equilibrium rate constants(ie. non-special)

         lrt = q(ik, 1) * rtw(ik)
         r(ik, 2, 1) = pgscav(ik, 1, 1)
         r(ik, 3, 1) = r(ik, 2, 1) / (h2a * h2b * lrt)
!  limit forward rate for ozone to reduce round-off error
         r(ik, 4, 1) = min(pgscav(ik, 11, 1), 1.0e-8)
         r(ik, 5, 1) = r(ik, 4, 1) / (h4 * lrt)
         r(ik, 6, 1) = pgscav(ik, 2, 1)
         r(ik, 7, 1) = r(ik, 6, 1) / (h6 * lrt)
         r(ik, 8, 1) = pgscav(ik, 7, 1)
         r(ik, 9, 1) = r(ik, 8, 1) / (h8 * lrt)
         r(ik, 10, 1) = pgscav(ik, 3, 1)
         r(ik, 11, 1) = r(ik, 10, 1) / (h10 * lrt)
         r(ik, 12, 1) = pgscav(ik, 8, 1)
         r(ik, 13, 1) = r(ik, 12, 1) / (h12a * h12b * lrt)
         r(ik, 15, 1) = min(pgscav(ik, 12, 1), 1.0e-9)
         r(ik, 16, 1) = r(ik, 15, 1) / (h15a * h15b * lrt)
         r(ik, 23, 1) = pgscav(ik, 5, 1)
         r(ik, 24, 1) = pgscav(ik, 6, 1)
         r(ik, 25, 1) = pgscav(ik, 9, 1)

      end do

      call mach_incld_upaqr2(aq, ti, r, q, nptsnz)
       
!     section for q(:,2)>0.

!  ICE scavenging rates for particles and HNO3 & NH3
      do ik = 1, nptsnz
         if (q(ik, 2) > 0.0) then
            lrt = q(ik, 2) * rtw(ik)
            tsqi = ti(ik) * ti(ik)
            h15a = 10.0 ** (-10.66 + 760.0 * ti(ik) + 58000.0 * tsqi)
            h15b = 10.0 ** (-14.25 + 5190.0 * ti(ik) - 850000.0 * tsqi)
            r(ik, 8, 2) = pgscav(ik, 7, 2)
            r(ik, 12, 2) = pgscav(ik, 8, 2)
            r(ik, 15, 2) = 1.0e-9
            r(ik, 16, 2) = r(ik, 15, 2) / (h15a * h15b * lrt)
         end if
      end do
!
      return
   end subroutine mach_incld_upaqr1
   
 end module mach_incld_upaqr_mod
