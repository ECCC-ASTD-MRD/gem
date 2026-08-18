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
! Fichier/File   : mach_incld_rates.ftn90
! Creation       : S. Menard, S. Gravel, GEM-MACH, June 2008
! Description    : Compute derivatives in chemeq format for chemical kinetics and
!                  mass transfer (excluding intra-hydrometeor fluxes). This version
!                  assumes daq, paq, & cdaq are only used when q > 0.
!
! Extra info     : QA  ADOM  VERSION: ADOMIIB(CLEAN)  LEVEL: 06/09/89  RATES  ENSR(AES)
!
! Arguments  IN
!               ideriv -->
!               q      --> Liquid water conc of cloudwater and Ice/snow
!               g      --> Gas/part species conc (ppm)
!               aq     --> Aq. species conc in cloudwater and ice/snow (m/l)
!               r      -->
!               b      -->
!
!            OUT
!               DG     -->
!               PG     -->
!               CDG    -->
!               DAQ    -->
!               PAQ    -->
!               CDAQ(2)-->
!
!=============================================================================
!
!!if_on
subroutine mach_incld_rates(gaz_conc, dg, pg, cdg, aq, daq, paq, cdaq, &
                            r, b, q, ideriv)
   use mach_cam_utils_mod,     only: maxnsg, maxnsaq
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent (in) :: ideriv
   real(kind=4),    intent (in) :: q       (2)
   real(kind=4),    intent (in) :: b       (5, 2)
   real(kind=4),    intent (in) :: r       (25, 2)
   real(kind=4),    intent(out) :: dg      (maxnsg)
   real(kind=4),    intent(out) :: pg      (maxnsg)
   real(kind=4),    intent(out) :: cdg     (maxnsg)
   real(kind=4),    intent (in) :: aq      (maxnsaq, 2)
   real(kind=4),    intent(out) :: daq     (maxnsaq, 2)
   real(kind=4),    intent(out) :: paq     (maxnsaq, 2)
   real(kind=4),    intent(out) :: cdaq    (maxnsaq, 2)
   real(kind=4),    intent (in) :: gaz_conc(maxnsg)
!!if_off
!
! Local variables
   integer(kind=4)              :: ngend, naqend, i
!
!  central operator aqueous phase chemical mechanism   08.23.84
!  gas kinetics
!  none in this version    set     dg = pg = zero

   ngend = maxnsg
   naqend = maxnsaq
   if (ideriv < 0) ngend = 3
   if (ideriv < 0) naqend = 4
   do i = 1, ngend
      dg(i) = 0.
      pg(i) = 0.
   end do
   do i = 1, naqend
      daq(i, 1) = 0.
      paq(i, 1) = 0.
      daq(i, 2) = 0.
      paq(i, 2) = 0.
   end do

!  cloudwater
   if (q(1) > 0.0) then
      daq(1, 1) = r(3, 1) * aq(9, 1) + r(19, 1) * aq(12, 1) +        &
                  r(20, 1) * aq(2, 1) + r(21, 1) * aq(3, 1) + r(22, 1)
      paq(1, 1) = b(1, 1) * r(2, 1) * gaz_conc(1)
      daq(2, 1) = r(7, 1) + r(20, 1) * aq(1, 1)
      paq(2, 1) = b(1, 1) * r(6, 1) * gaz_conc(2)
      daq(3, 1) = r(11, 1) + r(21, 1) * aq(1, 1)
      paq(3, 1) = b(1, 1) * r(10, 1) * gaz_conc(3)
      if (ideriv <= 0) then
         daq(4, 1) = 0.
         paq(4, 1) = (r(19, 1) * aq(12, 1) + r(20, 1) * aq(2, 1) +   &
                      r(21, 1) * aq(3, 1) + r(22, 1)) * aq(1, 1)
      else
         daq(4, 1) = 0.
         paq(4, 1) = b(1, 1) * r(1, 1) * gaz_conc(4) +               &
                     b(1, 1) * r(23, 1) * gaz_conc(5) +              &
                     b(1, 1) * r(24, 1) * gaz_conc(6)
         daq(5, 1) = r(9, 1) * aq(9, 1)
         paq(5, 1) = b(1, 1) * r(8, 1) * gaz_conc(7) +               &
                     b(1, 1) * r(25, 1) * gaz_conc(9)
         daq(6, 1) = r(13, 1) * aq(10, 1)
         paq(6, 1) = b(1, 1) * r(12, 1) * gaz_conc(8) +              &
                     b(1, 1) * r(23, 1) * gaz_conc(5) +              &
                     b(3, 1) * r(24, 1) * gaz_conc(6) +              &
                     b(1, 1) * r(25, 1) * gaz_conc(9)
         daq(7, 1) = 0.
         paq(7, 1) = b(4, 1) * r(14, 1) * gaz_conc(10)
         daq(8, 1) = r(16, 1) * aq(9, 1)
         paq(8, 1) = b(1, 1) * r(15, 1) * gaz_conc(12)
      end if
      do i = 1, naqend
         cdaq(i, 1) = paq(i, 1) - daq(i, 1) * aq(i, 1)
      end do

!  gas/part (from cloudwater)
      dg(1) = r(2, 1)
      pg(1) = b(2, 1) * r(3, 1) * aq(1, 1) * aq(9, 1)
      dg(2) = r(6, 1)
      pg(2) = b(2, 1) * r(7, 1) * aq(2, 1)
      dg(3) = r(10, 1)
      pg(3) = b(2, 1) * r(11, 1) * aq(3, 1)
      if (ideriv >= 0) then
         dg(4) = r(1, 1)
         dg(5) = r(23, 1)
         dg(6) = r(24, 1)
         dg(7) = r(8, 1)
         pg(7) = b(2, 1) * r(9, 1) * aq(5, 1) * aq(9, 1)
         dg(8) = r(12, 1)
         pg(8) = b(2, 1) * r(13, 1) * aq(6, 1) * aq(10, 1)
         dg(9) = r(25, 1)
         dg(10) = r(14, 1)
      end if
   end if

! --- ice/snow
   if (q(2) > 0.0) then
      if (ideriv > 0) then
         paq(4, 2) = b(1, 2) * r(1, 2) * gaz_conc(4) +               &
                     b(1, 2) * r(23, 2) * gaz_conc(5) +              &
                     b(1, 2) * r(24, 2) * gaz_conc(6)
         paq(5, 2) = b(1, 2) * r(8, 2) * gaz_conc(7) +               &
                     b(1, 2) * r(25, 2) * gaz_conc(9)
         paq(6, 2) = b(1, 2) * r(12, 2) * gaz_conc(8) +              &
                     b(1, 2) * r(23, 2) * gaz_conc(5) +              &
                     b(3, 2) * r(24, 2) * gaz_conc(6) +              &
                     b(1, 2) * r(25, 2) * gaz_conc(9)
         paq(7, 2) = b(4, 2) * r(14, 2) * gaz_conc(10)
         daq(8, 2) = r(16, 2) * aq(9, 2)
         paq(8, 2) = b(1, 2) * r(15, 2) * gaz_conc(12)
      end if 
      do i = 1, naqend
         cdaq(i, 2) = paq(i, 2) - daq(i, 2) * aq(i, 2)
      end do

! --- gas/part (from ice/snow)
      if (ideriv >= 0) then
         dg(4) = dg(4) + r(1, 2)
         dg(5) = dg(5) + r(23, 2)
         dg(6) = dg(6) + r(24, 2)
         dg(7) = dg(7) + r(8, 2)
         dg(8) = dg(8) + r(12, 2)
         dg(9) = dg(9) + r(25, 2)
         dg(10) = dg(10) + r(14, 2)
      end if
   end if

! --- compute net derivative for gases
   do i = 1, ngend
      cdg(i) = pg(i) - dg(i) * gaz_conc(i)
   end do

return
end subroutine mach_incld_rates
