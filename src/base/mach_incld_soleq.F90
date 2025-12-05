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
! Fichier/File   : mach_incld_soleq.ftn90
! Creation       : S. Menard, W. Gong, S. Gravel, B. Pabla, GEM-MACH, June 2008
! Description    : Central operator aqueous phase chemical mechanism
!
! Extra info     : See chemical reactions below
!
! Arguments  IN
!               q(:,1)  --> cloudwater conc (g-w/m3 air)
!               q(:,2)  --> Ice/Snow conc (g-w/m3 air)
!               ti      --> Inverse temperature (1/K)
!               rt      --> TEMP*RGAS/WDEN
!              pgscav   --> Gas/part scavenging (or diffusion)
!               nptsnz  --> Total number of grids to integrate
!
!            INOUT
!               gaz_conc--> Gas/Part species conc (ppm)
!               aq      --> Aq. species conc in cloudwater and ice/snow (m/l)
!               b       --> Aqueous variable coefficients
!
!=============================================================================
!
!!if_on
subroutine mach_incld_soleq(gaz_conc, aq, b, r, nptsnz, iaq, ncw)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
!!if_off
   use mach_incld_headers_mod, only: mach_incld_funeq, &
                                     mach_incld_findh, mach_incld_concmp
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nptsnz
   integer(kind=4), intent   (in) :: ncw
   integer(kind=4), intent   (in) :: iaq     (nptsnz)
   real(kind=4),    intent(inout) :: b       (nptsnz, 5, 2)
   real(kind=4),    intent(inout) :: gaz_conc(nptsnz, maxnsg)
   real(kind=4),    intent(inout) :: aq      (nptsnz, maxnsaq, 2)
   real(kind=4),    intent   (in) :: r       (nptsnz, 25, 2)
!!if_off
!
! Local variables
!
   integer(kind=4)                       :: ik, ig, iq
   real(kind=4)                          :: tot1, tot2, tot3, dch1, resid
   real(kind=4), dimension(ncw, 12)      :: c
   real(kind=4), dimension(ncw, 5)       :: tcmp
   real(kind=4), dimension(ncw, 25)      :: r1
   real(kind=4), dimension(ncw, 5)       :: b1
   real(kind=4), dimension(ncw, maxnsaq) :: aq1
   real(kind=4), dimension(ncw, maxnsg)  :: gnew

!    CENTRAL OPERATOR AQUEOUS PHASE CHEMICAL MECHANISM
!
!                 SPECIES LIST
!       AQUEOUS                    GAS/PART.
!       -------                    ---------
!       1.  HSO3-                  1.  SO2G
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
!              EQUILIBRIUM  REACTION LIST
!              --------------------------
!
!  1.     1.SO4P           -->    B1*SO4= +   B3*H+   (WITH K= 0)
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
! 14.     1.DUST           -->    B5*FEMN +   B4*HCO3 + B4*CAT1+ (K = 0)
! 15.     1.CO2G           -->    B1*HCO3 +   B1*H+
! 16.     1.HCO3 +  1.H+   -->    B2*CO2G
! 17.     1.H+   +  1.OH-  -->   1.00H2OA
! 18.     1.H2OA           -->   1.00H+   +  1.00OH-
! 19.     1.SO4P2          -->    B1*SO4= +   B1*H+ + B1*NH4+
! 20.     1.SO4P3          -->    B1*SO4= + B3*NH4+
! 21.     1.NNO3P          -->    B1*NO3- + B1*NH4+
!
!   SPECIES LIST USED BY EQUILIBRIUM SOLVER
!
! 01    SO2 GAS                 SO2G
! 02    NITRIC ACID GAS         HNO3
! 03    NH3 GAS                 NH3G
! 04    CO2 GAS                 CO2G
! 05    AQ PROTONS              H+
! 06    AQ BISULFITE            HSO3
! 07    AQ NITRATE              NO3-
! 08    AQ AMMONIUM             NH4+
! 09    AQ BICARBONATE          HCO3
! 10    AQ HYDROXIDE            OH-
! 11    AQ SULFATE              SO4=
! 12    AQ CATIONS              CAT1
!
   do ig = 1, maxnsg
      gnew(:, ig) = pack(gaz_conc(:, ig), iaq > 1)
   end do

   do iq = 1, maxnsaq
      aq1(:, iq) = pack(aq(:, iq, 1), iaq > 1)
   end do

   do iq = 1, 25
      r1(:, iq) = pack(r(:, iq, 1), iaq > 1)
   end do

   do iq = 1, 5
      b1(:, iq) = pack(b(:, iq, 1), iaq > 1)
   end do

!  pack local array of aq & gas species for equilibrium solution
   do ik = 1, ncw
      c(ik, 1) = gnew(ik, 1)
      c(ik, 2) = gnew(ik, 7)
      c(ik, 3) = gnew(ik, 8)
      c(ik, 4) = gnew(ik, 12)
      c(ik, 5) = aq1(ik, 9)
      c(ik, 6) = aq1(ik, 1)
      c(ik, 7) = aq1(ik, 5)
      c(ik, 8) = aq1(ik, 6)
      c(ik, 9) = aq1(ik, 8)
      c(ik, 10) = aq1(ik, 10)
      c(ik, 11) = aq1(ik, 4)
      c(ik, 12) = aq1(ik, 7)
   end do

   do ik = 1, ncw
!  o3 equilibrium
      aq1(ik, 12) = b1(ik, 1) * r1(ik, 4) * gnew(ik, 11) / r1(ik, 5)
      tot1 = gnew(ik, 11) + b1(ik, 2) * aq1(ik, 12)
      aq1(ik, 12) = ((r1(ik, 4) / (r1(ik, 4) + r1(ik, 5))) / b1(ik, 2)) * tot1
      gnew(ik, 11) = tot1 - b1(ik, 2) * aq1(ik, 12)
!  h2o2 equilibrium
      tot2 = gnew(ik, 2) + b1(ik, 2) * aq1(ik, 2)
      aq1(ik, 2) = ((r1(ik, 6) / (r1(ik, 6) + r1(ik, 7))) / b1(ik, 2)) * tot2
      gnew(ik, 2) = tot2 - b1(ik, 2) * aq1(ik, 2)
!  rooh equilibrium
      tot3 = gnew(ik, 3) + b1(ik, 2) * aq1(ik, 3)
      aq1(ik, 3) = ((r1(ik, 10) / (r1(ik, 10) + r1(ik, 11))) / b1(ik, 2)) * tot3
      gnew(ik, 3) = tot3 - b1(ik, 2) * aq1(ik, 3)
   end do

!  initial component totals
   call mach_incld_funeq(tcmp, c, b1, ncw)

!  find proton concentration
!  Add [HCO3-] for initial H+ estimate to avoid unrealistically large H+
!  estimate when [SO4=], [NO3-], and [NH4+] are very small
!  use initial estimate of H+ (charge balance of major ions)
!  for [HCO3-] evaluation:

   do ik = 1, ncw
      c(ik, 5) = c(ik, 5) + r1(ik, 15) * tcmp(ik, 4) / b1(ik, 2) /   &
                  (r1(ik, 15) + r1(ik, 16) * c(ik, 5))
!  setting a lower limit for proton concentration
      c(ik, 5) = max(c(ik, 5), 3.0e-07)
   end do

   call mach_incld_findh(tcmp, c, r1, b1, ncw)

!  find other species concentrations
   call mach_incld_concmp(tcmp, c, r1, b1, ncw)

   do ik = 1, ncw
!  unpack local c array into  aq1 & g arrays
!  cloudwater
      aq1(ik, 1) = c(ik, 6)
      aq1(ik, 5) = c(ik, 7)
      aq1(ik, 6) = c(ik, 8)
      aq1(ik, 8) = c(ik, 9)
      aq1(ik, 9) = c(ik, 5)
      aq1(ik, 10) = c(ik, 10)
!  gases
      gnew(ik, 1) = c(ik, 1)
      gnew(ik, 7) = c(ik, 2)
      gnew(ik, 8) = c(ik, 3)
!  insure charge balance / adjust protons as needed
      dch1 = -2.0 * aq1(ik, 4) - aq1(ik, 1) - aq1(ik, 5) &
                    - aq1(ik, 8) - aq1(ik, 10) &
                    + aq1(ik, 6) + aq1(ik, 9) + aq1(ik, 7)
      aq1(ik, 9) = aq1(ik, 9) - dch1
   end do

   do ik = 1, ncw
!  added to avoid negative proton concentrations.
      if (aq1(ik, 9) <= 0.0) then
!  add residual to oh
         resid = 1.0e-7 - aq1(ik, 9)
         aq1(ik, 9) = 1.0e-7
         aq1(ik, 10) = resid
      else
         resid = 0.0
      end if
   end do

! Unpack the output fields 
   do ig = 1, maxnsg
      gaz_conc(:, ig) = unpack(gnew(:, ig), iaq > 1, gaz_conc(:, ig))
   end do
   do iq = 1, maxnsaq
      aq(:, iq, 1) = unpack(aq1(:, iq), iaq > 1, aq(:, iq, 1))
   end do
   do iq = 1, 5
      b(:, iq, 1) = unpack(b1(:, iq), iaq > 1, b(:, iq, 1))
   end do

   return
end subroutine mach_incld_soleq
