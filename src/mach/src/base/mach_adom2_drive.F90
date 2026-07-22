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
! Fichier/File   : mach_adom2_drive.ftn90
! Creation       : P. Makar, Feb 2007 for GEM-MACH
!                  A. Kallaur, 2003 for CHRONOS on IBM supercomputer
! Description    : Gas-phase chemistry solver.
!                  This subroutine solves a set of ordinary differential equations
!                  describing gas phase chemistry of nox/voc system using the chemeq
!                  method of young and boris.
!
! Extra info     : t.r.young and j.p.boris, 1977. "a numerical technique for solving
!                  stiff ordinary differential equations associated with the chemical
!                  kinetics of reactive-flow problems". the journal of physical
!                  chemistry, vol.81, no.25. (see uprate for the species list & chemical mechanism)
!			
!		   Test for positivity on Yg and Yc1g are recoded to make the code more robust. 
!		   Positivity check for Yc1g is now done for every iteration (k = 1, niter).
!		   S. Menard (Sept 2011). See MANTIS 1944 for more information
!	           https://ulysse.cmc.ec.gc.ca/mantis/view.php?id=1944	
!
!   2017 Feb.  Change name from mach_gas_drive to mach_adom2_drive
!              to differentiate from other mech. routines (Jack C)
!              now called by "mach_adom2yb_main.ftn90"
!   2018 Mar.  Remove unused NH3 and dust species, and rearrang YG; and
!              Replace maxnsa with nsi (Deji Akingunola)
!
! Arguments:  IN
!               yg - > vector of initial gas-phase concentrations in ppm.
!               ig, jg -> x, y grid indices (used in diagnostic messages in steadyg)
!               indx -> a flag = 1 for first start up call
!                               = 0 for subsequent calls
!               zen -> solar zenith angle (deg)
!               rgs -> reaction rate
!               bgs -> product coefficients
!             OUT
!                yg    -> vector of final gas-phase concentrations for f species
!                indx -> a flag=  0  for successful completion
!                               = -1  for integration failure
!                               = -10 for integration taking too many steps
!                               = -2x for stoping before division by zero
!            MODULE Parameters
!               nsi -> number of integrated species
!               nspec -> total number of species
!               nreac ->  Number of chemical reactions in the chemical solver (114)
!               nprcf ->  Number of variable product (stoichiometric coefficients) (30)
!============================================================================
!
!!if_on
subroutine mach_adom2_drive(yg, indx, zen, ig, jg, rgs, bgs)
   use mach_pkg_adom2_mod, only: nspec, nreac, nprcf
!!if_off
   use mach_pkg_adom2_mod, only: nsi, hstart, hsub, hmax, ymings, ymin2gs
   use chm_utils_mod,      only: chm_lun_out, chm_timestep
   implicit none
!!if_on
   integer(kind=4), intent(inout) :: indx
   real(kind=4),    intent   (in) :: zen
   integer(kind=4), intent   (in) :: ig
   integer(kind=4), intent   (in) :: jg
   real(kind=4),    intent(inout) :: yg  (nspec)
   real(kind=4),    intent   (in) :: rgs (nreac)
   real(kind=4),    intent   (in) :: bgs (nprcf)
!!if_off
!
!     declare local variables
!
!     ntyp0 --|
!     ntyp1   |
!     ntyp2   |---> Used in Uprate (Rate constant types)
!     ntyp4   |
!     ntyp5 --|
!
!     nzenth -|
!     nzenm1  |
!             |---> Used in Uprate
!     nzhite  |     (Dimensions forPhotolytic Rate lookup Table)
!     nzhim1 -|
!
   integer(kind=4), parameter :: niter  = 3    ! Maximum number of corrector iterations
   integer(kind=4), parameter :: nstchk = 1    ! Check interval for stiff species
   integer(kind=4), parameter :: nstmax = 8000 ! Maximum number of iteration
!  specify average (eps) & maximum (teps) convergence criteria
   real(kind=4), parameter    :: epsng = 0.002, epssg = 0.002, tepsng = 0.01, tepssg = 0.01
   real(kind=4), parameter    :: dtmin = 1.e-10    ! Minimum dt(under any circumstances)
   real(kind=4), parameter    :: smf = 1.0e-37     ! ?
   real(kind=4)               :: ypreg(nspec)
   real(kind=4)               :: ydpreg(nsi)
   real(kind=4)               :: apreg(nsi)
   real(kind=4)               :: bpreg(nsi)
   real(kind=4)               :: yc1g(nspec)
   real(kind=4)               :: ydc1g(nsi)
   real(kind=4)               :: ac1g(nsi)
   real(kind=4)               :: bc1g(nsi)
   real(kind=4)               :: yc2g(nspec)
   real(kind=4)               :: ers(nsi)
   real(kind=4)               :: ern(nsi)
   integer(kind=4)            :: jst(nsi)
   integer(kind=4)            :: jnst(nsi)
   integer(kind=4)            :: np1g, nnormg, nstifg
   integer(kind=4)            :: i, j, k, kiter, iconv
   integer(kind=4)            :: nsteps
   real(kind=4)               :: errng, errsg, ersmsg, ersmng
   real(kind=4)               :: time, dtmax
   real(kind=4)               :: tout, xcrit
   real(kind=4)               :: dum2, xtemp, dthalf, dts, divxtemp2, xtemp2
!***
!     Declare variables used in the in-lined
!     subroutine sections of steadyg.
!***
   real(kind=4)               :: a(2, 2)
   real(kind=4)               :: prod, dest, rhs, det, deti

!************************************************************************
!                               code begins                             *
!************************************************************************
   np1g = nsi + 1
   time = 0.0
   tout = chm_timestep / 60.0
   nsteps = 0

!  select initial and maximum time step sizes
   if (indx == 0 .and. zen <= 90.) then
      dts = hsub(1)
   end if
   if (indx == 0 .and. zen > 90.) then
      dts = hsub(2)
   end if
   if (zen <= 90.0) then
      dtmax = hmax(1)
   end if
   if (zen > 90.0) then
      dtmax = hmax(2)
   end if
   if (indx == 1) then
      dts = hstart
   end if

! --- check positivity for all yg(i)
   yg = max(yg , ymin2gs)
!
!      in-line code version of suborutine:

!      subroutine steadygsc(c, igrid, jgrid)

!      called from drivesc_il in the following form:
!      call steadygsc(yg, ig, jg)

! --- o(1d)
   prod = rgs(16) * yg(5)
   dest = rgs(17) * yg(43) + rgs(18) * yg(45)
   yg(29) = prod / dest

! --- crgs1
   prod = (0.40 * rgs(85) * yg(11) +  &
            bgs(25) * rgs(89) * yg(12)) * yg(5)
   dest = rgs(92) * yg(1) +  &
            rgs(94) * yg(43) +  &
            rgs(96) * yg(15) +  &
            rgs(98) * yg(16)
   yg(39) = prod / dest

! --- crgs2
   prod = (bgs(26) * rgs(89) * yg(12) +  &
            0.20 * rgs(110) * yg(24)) * yg(5)
   dest = rgs(93) * yg(1) +  &
            rgs(95) * yg(43) +  &
            rgs(97) * yg(15) +  &
            rgs(99) * yg(16)
   yg(40) = prod / dest

! --- use the max of last concentration or 1.e-29 ppm for initial
! --- guess of [no3]
   yg(31) = max(yg(31), 1.0e-29)

! --- o(3p)
   prod = rgs(1) * yg(4) +  &
            rgs(14) * yg(31) +  &
            rgs(15) * yg(5) +  &
            rgs(18) * yg(29) * yg(45)
   dest = rgs(2) * yg(44) * yg(45) +  &
         (rgs(3) + rgs(4)) * yg(4) +  &
            rgs(86) * yg(11) +  &
            rgs(90) * yg(12) +  &
            rgs(111) * yg(24)
   yg(30) = prod / dest

! --- use the max of last concentration or 1.e-29 ppm for initial
! --- guess of [hno4]
   yg(33) = max(yg(33), 1.e-29)

! --- oh
   prod = 2.00 * rgs(17) * yg(29) * yg(43) +  &
            rgs(20) * yg(22) +  &
            rgs(23) * yg(7) +  &
         (rgs(27) * yg(3) +  &
            rgs(32) * yg(5)) * yg(25) +  &
            2.00 * rgs(37) * yg(6) +  &
            rgs(48) * yg(20) +  &
         (bgs(14) * rgs(89) * yg(12) +  &
            0.10 * rgs(110) * yg(24)) * yg(5) +  &
            0.50 * rgs(114) * yg(20) * yg(34)
   dest =   rgs(19) * yg(3) +  &
            rgs(22) * yg(4) +  &
            rgs(24) * yg(7) +  &
            rgs(25) * yg(28) +  &
            rgs(26) * yg(5) +  &
            rgs(31) * yg(33) +  &
            rgs(38) * yg(6) +  &
            rgs(43) * yg(1) +  &
            rgs(51) * yg(15) +  &
            rgs(54) * yg(16) +  &
            rgs(63) * yg(17) +  &
            rgs(65) * yg(18) +  &
            rgs(67) * yg(41) +  &
            rgs(68) * yg(42) +  &
            rgs(69) * yg(9) +  &
            rgs(70) * yg(10) +  &
            rgs(71) * yg(23) +  &
            rgs(84) * yg(11) +  &
            rgs(88) * yg(12) +  &
            rgs(100) * yg(13) +  &
            rgs(101) * yg(14) +  &
            rgs(102) * yg(19) +  &
            rgs(104) * yg(21) +  &
            rgs(109) * yg(24) +  &
            rgs(113) * yg(25) +  &
            rgs(114) * yg(20)
   yg(34) = prod / dest

! --- hno4
   prod = rgs(28) * yg(4) * yg(25)
   dest = rgs(29) + rgs(30) + rgs(31) * yg(34)
   yg(33) = prod / dest

! --- solve no3/n2o5 system
! --- no3
   rhs = -((rgs(4) * yg(30) +  &
               rgs(6) * yg(5)) * yg(4) +  &
               rgs(24) * yg(7) * yg(34))
   a(1, 2) = rgs(10)
   a(1, 1) = -(rgs(7) * yg(3) +  &
               (rgs(9) + rgs(12)) * yg(4) +  &
               rgs(13) + rgs(14) +  &
               (rgs(39) + rgs(40) * yg(45) +  &
               (rgs(41) + rgs(42)) * yg(43)) * yg(25) +  &
               rgs(52) * yg(15) +  &
               rgs(56) * yg(16) +  &
               rgs(66) * yg(18) +  &
               rgs(87) * yg(11) +  &
               rgs(91) * yg(12) +  &
               rgs(105) * yg(21) +  &
               rgs(112) * yg(24))

! --- n2o5
   a(2, 1) = rgs(9) * yg(4)
   a(2, 2) = -(rgs(10) + rgs(11) * yg(43))
   det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
   if (det == 0.0) then
      write(0, *)  '### Error in mach_adom2_drive ###'
      write(0, *)  '# subroutine steadyg was unable to solve'
      write(0, *)  '# the equations for no3 and n2o5'
      write(0, *)  '# det=0 in steadyg, i=', ig, ' j=', jg
      write(0, *)  '# a(1, 1)=', a(1, 1)
      write(0, *)  '# a(1, 2)=', a(1, 2)
      write(0, *)  '# a(2, 1)=', a(2, 1)
      write(0, *)  '# a(2, 2)=', a(2, 2)
      write(0, *)  '# rg( ) =', rgs
      write(0, *)  '# yg ( ) =', yg
      write(0, *)  '###             ABORT             ###'
      indx = -21
      return
   end if
   deti = 1.0 / det
   yg(31) = rhs * a(2, 2) * deti
   yg(32) = -rhs * a(2, 1) * deti

! --- ro2r
   prod = rgs(53) * yg(15) * yg(25) +  &
            rgs(55) * yg(16) +  &
            rgs(57) * yg(3) * yg(27) +  &
            rgs(62) * yg(17) +  &
         (1.50 * rgs(63) * yg(17) +  &
            rgs(67) * yg(41) +  &
            rgs(68) * yg(42) +  &
            rgs(69) * yg(9) +  &
            bgs(05) * rgs(70) * yg(10) +  &
            rgs(84) * yg(11) +  &
            rgs(88) * yg(12) +  &
            0.84 * rgs(100) * yg(13) +  &
            0.83 * rgs(101) * yg(14) +  &
            0.85 * rgs(104) * yg(21)) * yg(34) +  &
         (rgs(86) * yg(11) +  &
            bgs(21) * rgs(90) * yg(12)) * yg(30) +  &
            bgs(12) * rgs(89) * yg(5) * yg(12) +  &
            0.50 * rgs(114) * yg(20) * yg(34)
   dest = rgs(80) * yg(3) +  &
            rgs(81) * yg(25) +  &
            rgs(82) * yg(26) +  &
            rgs(83) * yg(27)
   yg(35) = prod / dest

! --- r2o2
   prod = (bgs(06) * rgs(70) * yg(10) +  &
            1.39 * rgs(71) * yg(23)) * yg(34) +  &
            (rgs(87) * yg(11) +  &
            rgs(91) * yg(12)) * yg(31) +  &
            (0.90 * rgs(109) * yg(34) +  &
            0.50 * rgs(111) * yg(30) +  &
            rgs(112) * yg(31)) * yg(24)
   dest = rgs(76) * yg(3) +  &
            rgs(77) * yg(25) +  &
            rgs(78) * yg(26) +  &
            rgs(79) * yg(27)
   yg(36) = prod / dest

! --- ro2n
   prod = (bgs(04) * rgs(70) * yg(10) +  &
            0.14 * rgs(104) * yg(21) +  &
            0.10 * rgs(109) * yg(24)) * yg(34)
   dest = rgs(72) * yg(3) +  &
            rgs(73) * yg(25) +  &
            rgs(74) * yg(26) +  &
            rgs(75) * yg(27)
   yg(37) = prod / dest

! --- bzo
   prod = rgs(105) * yg(21) * yg(31)
   dest = rgs(106) * yg(4) +  &
            rgs(107) * yg(25) +  &
            rgs(108)
   yg(38) = prod / dest

!---------------------------------------------------------------
!                 end steadygsc_yg  in-line
!---------------------------------------------------------------
!
!
!     initialize function
!
!      in-line version to replace subroutine:
!      call diffungsc(ng, nspec, yg, ydpreg, apreg, bpreg)
!      steadygsc_il_yg.cdk
!
!      in-line "yg" code version of subroutine:
!      subroutine diffungsc(ns, ntot, c, cdot, prod, dest)
!
!      called from mach_gas_driveil in the following form:
!      call diffungsc(ng, nspec, yg, ydpreg, apreg, bpreg)
!
!----------------------------------------------------------------------
!                          begin code
!----------------------------------------------------------------------

!....... determine chemical transformation rates .......
!
!...... compute chemical derivatives in chemeq format
!       ie.  dc/dt = cdot = prod - dest*c
!----------------------------------------------------------------------

! --- so2
   apreg(01) = 0.0
   bpreg(01) = rgs(43) * yg(34) +  &
               rgs(92) * yg(39) +  &
               rgs(93) * yg(40)

! --- so4
   apreg(02) = bpreg(01) * yg(1)
   bpreg(02) = 0.0

! --- no
   apreg(03) = (rgs(1) +  &
                  rgs(3) * yg(30) +  &
                  rgs(12) * yg(31)) * yg(4) +  &
                  rgs(13) * yg(31) +  &
                  rgs(20) * yg(22)
   bpreg(03) = rgs(5) * yg(5) +  &
               rgs(7) * yg(31) +  &
               2.00 * rgs(8) * yg(3) * yg(44) +  &
               rgs(19) * yg(34) +  &
               rgs(27) * yg(25) +  &
               rgs(57) * yg(27) +  &
               rgs(72) * yg(37) +  &
               rgs(76) * yg(36) +  &
               rgs(80) * yg(35)

! --- no2
   apreg(04) = (rgs(5) * yg(5) +  &
            2.00 * (rgs(7) * yg(31) +  &
                  rgs(8) * yg(3) * yg(44)) +  &
                  rgs(27) * yg(25) +  &
                  rgs(57) * yg(27) +  &
                  rgs(76) * yg(36) +  &
                  rgs(80) * yg(35)) * yg(3) +  &
                  rgs(10) * yg(32) +  &
                  rgs(14) * yg(31) +  &
                  rgs(23) * yg(7) +  &
               (rgs(29) +  &
                  rgs(30) +  &
                  rgs(31) * yg(34)) * yg(33) +  &
                  rgs(61) * yg(8) +  &
                  rgs(71) * yg(23) * yg(34) +  &
               (rgs(87) * yg(11) +  &
                  rgs(91) * yg(12) +  &
                  rgs(112) * yg(24)) * yg(31)
   bpreg(04) = rgs(1) +  &
               (rgs(3) +  &
               rgs(4)) * yg(30) +  &
               rgs(6) * yg(5) +  &
               rgs(9) * yg(31) +  &
               2.00 * rgs(21) * yg(43) +  &
               rgs(22) * yg(34) +  &
               rgs(28) * yg(25) +  &
               rgs(58) * yg(27) +  &
               rgs(106) * yg(38)

! --- o3
   apreg(05) = rgs(2) * yg(30) * yg(44) * yg(45)
   bpreg(05) = rgs(5) * yg(3) +  &
               rgs(6) * yg(4) +  &
               rgs(15) +  &
               rgs(16) +  &
               rgs(26) * yg(34) +  &
               rgs(32) * yg(25) +  &
               rgs(85) * yg(11) +  &
               rgs(89) * yg(12) +  &
               rgs(110) * yg(24)

! --- h2o2
   apreg(06) = (rgs(33) +  &
                  rgs(34) * yg(45) +  &
               (rgs(35) +  &
                  rgs(36)) * yg(43)) * yg(25) ** 2
   bpreg(06) = rgs(37) + rgs(38) * yg(34)

! --- hno3
   apreg(07) = (2.00 * rgs(11) * yg(32) +  &
                  rgs(21) * yg(4)) * yg(43) +  &
                  rgs(22) * yg(4) * yg(34) +  &
               ((rgs(39) + rgs(40) * yg(45)) * yg(25) +  &
               (rgs(41) + rgs(42)) * yg(25) * yg(43) +  &
                  rgs(52) * yg(15) +  &
                  rgs(56) * yg(16) +  &
                  rgs(66) * yg(18) +  &
                  rgs(105) * yg(21)) * yg(31)
   bpreg(07) = rgs(23) + rgs(24) * yg(34)

! --- pan
   apreg(08) = rgs(58) * yg(4) * yg(27)
   bpreg(08) = rgs(61)

! --- c3h8
   apreg(09) = 0.0
   bpreg(09) = rgs(69) * yg(34)

! --- alka
   apreg(10) = 0.0
   bpreg(10) = rgs(70) * yg(34)

! --- ethe
   apreg(11) = (rgs(109) * yg(34) +  &
                  0.50 * rgs(110) * yg(5) +  &
                  rgs(111) * yg(30)) * yg(24)
   bpreg(11) = rgs(84) * yg(34) +  &
               rgs(85) * yg(5) +  &
               rgs(86) * yg(30) +  &
               rgs(87) * yg(31)

! --- alke
   apreg(12) = 0.0
   bpreg(12) = rgs(88) * yg(34) +  &
               rgs(89) * yg(5) +  &
               rgs(90) * yg(30) +  &
               rgs(91) * yg(31)
! --- tolu
   apreg(13) = 0.0
   bpreg(13) = rgs(100) * yg(34)

! --- arom
   apreg(14) = 0.0
   bpreg(14) = rgs(101) * yg(34)

! --- hcho
   apreg(15) = rgs(55) * yg(16) +  &
              (rgs(57) * yg(3) +  &
               rgs(59) * yg(25) +  &
               2.00 * rgs(60) * yg(27) +  &
               rgs(75) * yg(37) +  &
               rgs(79) * yg(36) +  &
               rgs(83) * yg(35)) * yg(27) +  &
              (0.50 * rgs(63) * yg(17) +  &
               rgs(67) * yg(41) +  &
               bgs(01) * rgs(70) * yg(10) +  &
               0.16 * rgs(71) * yg(23) +  &
               1.56 * rgs(84) * yg(11) +  &
               0.11 * rgs(100) * yg(13) +  &
               bgs(24) * rgs(101) * yg(14)) * yg(34) +  &
              (rgs(85) * yg(5) +  &
               rgs(86) * yg(30) +  &
               2.00 * rgs(87) * yg(31)) * yg(11) +  &
              (bgs(08) * rgs(88) * yg(34) +  &
               bgs(10) * rgs(89) * yg(5) +  &
               bgs(18) * rgs(90) * yg(30) +  &
               bgs(08) * rgs(91) * yg(31)) * yg(12) +  &
               rgs(92) * yg(1) * yg(39) +  &
              (rgs(109) * yg(34) +  &
               rgs(110) * yg(5) +  &
               rgs(112) * yg(31)) * yg(24)
   bpreg(15) = rgs(49) +  &
               rgs(50) +  &
               rgs(51) * yg(34) +  &
               rgs(52) * yg(31) +  &
               rgs(53) * yg(25) +  &
               rgs(96) * yg(39) +  &
               rgs(97) * yg(40)

! --- ald2
   apreg(16) = rgs(62) * yg(17) +  &
              (0.50 * rgs(63) * yg(17) +  &
               rgs(68) * yg(42) +  &
               0.30 * rgs(69) * yg(9) +  &
               bgs(02) * rgs(70) * yg(10) +  &
               1.52 * rgs(71) * yg(23) +  &
               0.22 * rgs(84) * yg(11) +  &
               bgs(09) * rgs(88) * yg(12)) * yg(34) +  &
              (bgs(11) * rgs(89) * yg(5) +  &
               bgs(19) * rgs(90) * yg(30) +  &
               bgs(09) * rgs(91) * yg(31)) * yg(12) +  &
               rgs(93) * yg(1) * yg(40) +  &
              (0.20 * rgs(109) * yg(34) +  &
               0.40 * rgs(110) * yg(5) +  &
               rgs(111) * yg(30) +  &
               rgs(112) * yg(31)) * yg(24)
   bpreg(16) = rgs(54) * yg(34) +  &
               rgs(55) +  &
               rgs(56) * yg(31) +  &
               rgs(98) * yg(39) +  &
               rgs(99) * yg(40)

! --- mek
   apreg(17) = (0.50 * rgs(69) * yg(9) +  &
                bgs(03) * rgs(70) * yg(10) +  &
                0.14 * rgs(71) * yg(23)) * yg(34) +  &
               (rgs(73) * yg(25) +  &
                rgs(74) * yg(26) +  &
                rgs(75) * yg(27)) * yg(37) +  &
                bgs(17) * rgs(90) * yg(12) * yg(30)
   bpreg(17) = rgs(62) + rgs(63) * yg(34)

! --- mgly
   apreg(18) = (0.13 * rgs(100) * yg(13) +  &
                bgs(23) * rgs(101) * yg(14) +  &
                0.20 * rgs(104) * yg(21)) * yg(34) +  &
               (0.27 * rgs(109) * yg(34) +  &
                0.20 * rgs(110) * yg(5)) * yg(24)
   bpreg(18) = rgs(64) +  &
                rgs(65) * yg(34) +  &
                rgs(66) * yg(31)

! --- dial
   apreg(19) = (0.40 * rgs(100) * yg(13) +  &
                bgs(22) * rgs(101) * yg(14)) * yg(34)
   bpreg(19) = rgs(102) * yg(34) + rgs(103)

! --- rooh
   apreg(20) = (rgs(59) * yg(27) +  &
                rgs(73) * yg(37) +  &
                rgs(77) * yg(36) +  &
                rgs(81) * yg(35)) * yg(25)
   bpreg(20) = rgs(48) + rgs(114) * yg(34)

! --- cres
   apreg(21) = (0.16 * rgs(100) * yg(13) +  &
                0.17 * rgs(101) * yg(14)) * yg(34)
   bpreg(21) = 0.92 * rgs(104) * yg(34) +  &
                0.50 * rgs(105) * yg(31)

! --- hono
   apreg(22) = rgs(19) * yg(3) * yg(34) +  &
               rgs(21) * yg(4) * yg(43)
   bpreg(22) = rgs(20)

! --- rno3
   apreg(23) = rgs(72) * yg(3) * yg(37) +  &
               rgs(106) * yg(4) * yg(38)
   bpreg(23) = rgs(71) * yg(34)

! --- isop
   apreg(24) = 0.0
   bpreg(24) = rgs(109) * yg(34) +  &
               rgs(110) * yg(5) +  &
               rgs(111) * yg(30) +  &
               rgs(112) * yg(31)

! --- ho2
   apreg(25) = (rgs(25) * yg(28) +  &
                rgs(26) * yg(5) +  &
                rgs(38) * yg(6) +  &
                rgs(43) * yg(1) +  &
                0.16 * rgs(100) * yg(13) +  &
                0.17 * rgs(101) * yg(14)) * yg(34) +  &
               (rgs(29) + rgs(30)) * yg(33) +  &
                rgs(48) * yg(20) +  &
               (2.00 * rgs(49) +  &
                rgs(51) * yg(34) +  &
                rgs(52) * yg(31)) * yg(15) +  &
                rgs(55) * yg(16) +  &
               (2.00 * rgs(60) * yg(27) +  &
                rgs(75) * yg(37) +  &
                rgs(79) * yg(36) +  &
                rgs(83) * yg(35)) * yg(27) +  &
                rgs(64) * yg(18) +  &
                0.50 * rgs(74) * yg(26) * yg(37) +  &
               (rgs(80) * yg(3) +  &
                0.50 * rgs(82) * yg(26)) * yg(35) +  &
               (0.12 * rgs(85) * yg(5) +  &
                rgs(86) * yg(30)) * yg(11) +  &
               (bgs(13) * rgs(89) * yg(5) +  &
                bgs(20) * rgs(90) * yg(30)) * yg(12) +  &
                rgs(103) * yg(19) +  &
               (0.70 * rgs(109) * yg(34) +  &
                0.40 * rgs(110) * yg(5) +  &
                0.60 * rgs(111) * yg(30)) * yg(24)
   bpreg(25) = rgs(27) * yg(3) +  &
                rgs(28) * yg(4) +  &
                rgs(32) * yg(5) +  &
                2.00 * yg(25) * (rgs(33) +  &
                rgs(34) * yg(45) +  &
               (rgs(35) + rgs(36)) * yg(43)) +  &
               (rgs(39) + rgs(40) * yg(45) +  &
               (rgs(41) + rgs(42)) * yg(43)) * yg(31) +  &
                rgs(53) * yg(15) +  &
                rgs(59) * yg(27) +  &
                rgs(73) * yg(37) +  &
                rgs(77) * yg(36) +  &
                rgs(81) * yg(35) +  &
                rgs(107) * yg(38) +  &
                rgs(113) * yg(34)

! --- ro2
   apreg(26) = rgs(53) * yg(15) * yg(25) +  &
               rgs(55) * yg(16) +  &
               rgs(57) * yg(3) * yg(27) +  &
              (rgs(62) +  &
               1.50 * rgs(63) * yg(34)) * yg(17) +  &
              (rgs(67) * yg(41) +  &
               rgs(68) * yg(42) +  &
               rgs(69) * yg(9) +  &
               bgs(07) * rgs(70) * yg(10) +  &
               1.39 * rgs(71) * yg(23) +  &
               rgs(84) * yg(11) +  &
               0.84 * rgs(100) * yg(13) +  &
               0.83 * rgs(101) * yg(14) +  &
               rgs(104) * yg(21)) * yg(34) +  &
              (rgs(86) * yg(30) +  &
               rgs(87) * yg(31)) * yg(11) +  &
              (rgs(88) * yg(34) +  &
               bgs(12) * rgs(89) * yg(5) +  &
               bgs(21) * rgs(90) * yg(30) +  &
               rgs(91) * yg(31)) * yg(12) +  &
              (rgs(109) * yg(34) +  &
               0.50 * rgs(111) * yg(30) +  &
               rgs(112) * yg(31)) * yg(24) +  &
               0.50 * rgs(114) * yg(20) * yg(34)
   bpreg(26) = rgs(44) * yg(3) +  &
               rgs(45) * yg(25) +  &
               2.00 * rgs(46) * yg(26) +  &
               rgs(47) * yg(27)

! --- mco3
   apreg(27) = (rgs(54) * yg(34) +  &
                rgs(56) * yg(31)) * yg(16) +  &
                rgs(61) * yg(8) +  &
               (rgs(62) +  &
                rgs(63) * yg(34)) * yg(17) +  &
               (rgs(64) +  &
                rgs(65) * yg(34) +  &
                rgs(66) * yg(31)) * yg(18) +  &
               (rgs(102) * yg(34) +  &
                rgs(103)) * yg(19) +  &
                0.20 * rgs(109) * yg(24) * yg(34)
   bpreg(27) = rgs(57) * yg(3) +  &
               rgs(58) * yg(4) +  &
               rgs(59) * yg(25) +  &
               2.00 * rgs(60) * yg(27) +  &
               rgs(75) * yg(37) +  &
               rgs(79) * yg(36) +  &
               rgs(83) * yg(35)

!!cs_co_active
!
! --- co
    apreg(28) = rgs(49) * yg(15) +  &
                 rgs(50) * yg(15) +  &
                 rgs(51) * yg(15) * yg(34) +  &
                 rgs(52) * yg(15) * yg(31) +  &
                 rgs(55) * yg(16) +  &
                 rgs(64) * yg(18) +  &
                 rgs(65) * yg(18) * yg(34) +  &
                 rgs(66) * yg(18) * yg(31) +  &
                 rgs(85) * yg(11) * yg(5) * 0.42 +  &
                 rgs(86) * yg(11) * yg(30) +  &
                 rgs(89) * yg(12) * yg(5) * bgs(15) +  &
                 rgs(90) * yg(12) * yg(30) * bgs(16) +  &
                 rgs(100) * yg(13) * yg(34) * 0.11 +  &
                 rgs(101) * yg(14) * yg(34) * bgs(24) +  &
                 rgs(103) * yg(19)
    bpreg(28) = rgs(25) * yg(34)

! --- compute net derivative
   do i = 1, nsi
      ydpreg(i) = apreg(i) - bpreg(i) * yg(i)
   end do

!---------------------------------------------------------------
!                 end diffungsc_yg  in-line
!---------------------------------------------------------------
!
!
!      nfes = nfes + 1
!
!----------------------- begin integration step -----------------------
 18   if (nsteps /= 0) then
!
!     compute steady state species (once per step)
!
!      in-line version to replace subroutine:
!         call steadygsc(yg, ig, jg)
! steadygsc_il_yc1g.cdk
!
!***
!
!      in-line code version of suborutine:
!
!      subroutine steadygsc(c, igrid, jgrid)
!
!      called from mach_gas_driveil in the following form:
!      call steadygsc(yg, ig, jg)

! --- o(1d)
   prod = rgs(16) * yg(5)
   dest = rgs(17) * yg(43) + rgs(18) * yg(45)
   yg(29) = prod / dest

! --- crgs1
   prod = (0.40 * rgs(85) * yg(11) +  &
           bgs(25) * rgs(89) * yg(12)) * yg(5)
   dest = rgs(92) * yg(1) +  &
           rgs(94) * yg(43) +  &
           rgs(96) * yg(15) +  &
           rgs(98) * yg(16)
   yg(39) = prod / dest

! --- crgs2
   prod = (bgs(26) * rgs(89) * yg(12) +  &
           0.20 * rgs(110) * yg(24)) * yg(5)
   dest = rgs(93) * yg(1) +  &
           rgs(95) * yg(43) +  &
           rgs(97) * yg(15) +  &
           rgs(99) * yg(16)
   yg(40) = prod / dest

! --- use the max of last concentration or 1.e-29 ppm for initial
! --- guess of [no3]
   yg(31) = max(yg(31), 1.0e-29)

! --- o(3p)
   prod = rgs(1) * yg(4) +  &
          rgs(14) * yg(31) +  &
          rgs(15) * yg(5) +  &
          rgs(18) * yg(29) * yg(45)
   dest = rgs(2) * yg(44) * yg(45) +  &
         (rgs(3) + rgs(4)) * yg(4) +  &
          rgs(86) * yg(11) +  &
          rgs(90) * yg(12) +  &
          rgs(111) * yg(24)
   yg(30) = prod / dest

! --- use the max of last concentration or 1.e-29 ppm for initial
! --- guess of [hno4]
   yg(33) = max(yg(33), 1.e-29)

! --- oh
   prod = 2.00 * rgs(17) * yg(29) * yg(43) +  &
          rgs(20) * yg(22) +  &
          rgs(23) * yg(7) +  &
         (rgs(27) * yg(3) +  &
          rgs(32) * yg(5)) * yg(25) +  &
          2.00 * rgs(37) * yg(6) +  &
          rgs(48) * yg(20) +  &
         (bgs(14) * rgs(89) * yg(12) +  &
          0.10 * rgs(110) * yg(24)) * yg(5) +  &
          0.50 * rgs(114) * yg(20) * yg(34)
   dest = rgs(19) * yg(3) +  &
          rgs(22) * yg(4) +  &
          rgs(24) * yg(7) +  &
          rgs(25) * yg(28) +  &
          rgs(26) * yg(5) +  &
          rgs(31) * yg(33) +  &
          rgs(38) * yg(6) +  &
          rgs(43) * yg(1) +  &
          rgs(51) * yg(15) +  &
          rgs(54) * yg(16) +  &
          rgs(63) * yg(17) +  &
          rgs(65) * yg(18) +  &
          rgs(67) * yg(41) +  &
          rgs(68) * yg(42) +  &
          rgs(69) * yg(9) +  &
          rgs(70) * yg(10) +  &
          rgs(71) * yg(23) +  &
          rgs(84) * yg(11) +  &
          rgs(88) * yg(12) +  &
          rgs(100) * yg(13) +  &
          rgs(101) * yg(14) +  &
          rgs(102) * yg(19) +  &
          rgs(104) * yg(21) +  &
          rgs(109) * yg(24) +  &
          rgs(113) * yg(25) +  &
          rgs(114) * yg(20)
   yg(34) = prod / dest

! --- hno4
   prod = rgs(28) * yg(4) * yg(25)
   dest = rgs(29) + rgs(30) + rgs(31) * yg(34)
   yg(33) = prod / dest

! --- solve no3/n2o5 system
! --- no3
   rhs = -((rgs(4) * yg(30) +  &
               rgs(6) * yg(5)) * yg(4) +  &
               rgs(24) * yg(7) * yg(34))
   a(1, 2) = rgs(10)
   a(1, 1) = -(rgs(7) * yg(3) +  &
              (rgs(9) + rgs(12)) * yg(4) +  &
               rgs(13) + rgs(14) +  &
              (rgs(39) + rgs(40) * yg(45) +  &
              (rgs(41) + rgs(42)) * yg(43)) * yg(25) +  &
               rgs(52) * yg(15) +  &
               rgs(56) * yg(16) +  &
               rgs(66) * yg(18) +  &
               rgs(87) * yg(11) +  &
               rgs(91) * yg(12) +  &
               rgs(105) * yg(21) +  &
               rgs(112) * yg(24))

! --- n2o5
   a(2, 1) = rgs(9) * yg(4)
   a(2, 2) = -(rgs(10) + rgs(11) * yg(43))
   det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
   if (det == 0.0) then
      write(0, *)  '### Error in mach_mach_gas_driveil_omp ###'
      write(0, *)  '# subroutine steadyg was unable to solve'
      write(0, *)  '# the equations for no3 and n2o5'
      write(0, *)  '# det=0 in steadyg, i=', ig, ' j=', jg
      write(0, *)  '# a(1, 1)=', a(1, 1)
      write(0, *)  '# a(1, 2)=', a(1, 2)
      write(0, *)  '# a(2, 1)=', a(2, 1)
      write(0, *)  '# a(2, 2)=', a(2, 2)
      write(0, *)  '# rg( ) =', rgs
      write(0, *)  '# yg ( ) =', yg
      write(0, *)  '###             ABORT            ###'

!        only in in-line code
!        note: return is from drive to gaschem, not steadyg to drive!!!
      indx = -21
      return

   end if
   deti = 1.0 / det
   yg(31) = rhs * a(2, 2) * deti
   yg(32) = -rhs * a(2, 1) * deti

! --- ro2r
   prod = rgs(53) * yg(15) * yg(25) +  &
          rgs(55) * yg(16) +  &
          rgs(57) * yg(3) * yg(27) +  &
          rgs(62) * yg(17) +  &
         (1.50 * rgs(63) * yg(17) +  &
          rgs(67) * yg(41) +  &
          rgs(68) * yg(42) +  &
          rgs(69) * yg(9) +  &
          bgs(05) * rgs(70) * yg(10) +  &
          rgs(84) * yg(11) +  &
          rgs(88) * yg(12) +  &
          0.84 * rgs(100) * yg(13) +  &
          0.83 * rgs(101) * yg(14) +  &
          0.85 * rgs(104) * yg(21)) * yg(34) +  &
         (rgs(86) * yg(11) +  &
          bgs(21) * rgs(90) * yg(12)) * yg(30) +  &
          bgs(12) * rgs(89) * yg(5) * yg(12) +  &
          0.50 * rgs(114) * yg(20) * yg(34)
   dest = rgs(80) * yg(3) +  &
          rgs(81) * yg(25) +  &
          rgs(82) * yg(26) +  &
          rgs(83) * yg(27)
   yg(35) = prod / dest

! --- r2o2
   prod = (bgs(06) * rgs(70) * yg(10) +  &
           1.39 * rgs(71) * yg(23)) * yg(34) +  &
          (rgs(87) * yg(11) +  &
           rgs(91) * yg(12)) * yg(31) +  &
          (0.90 * rgs(109) * yg(34) +  &
           0.50 * rgs(111) * yg(30) +  &
           rgs(112) * yg(31)) * yg(24)
   dest = rgs(76) * yg(3) +  &
          rgs(77) * yg(25) +  &
          rgs(78) * yg(26) +  &
          rgs(79) * yg(27)
   yg(36) = prod / dest

! --- ro2n
   prod = (bgs(04) * rgs(70) * yg(10) +  &
           0.14 * rgs(104) * yg(21) +  &
           0.10 * rgs(109) * yg(24)) * yg(34)
   dest = rgs(72) * yg(3) +  &
          rgs(73) * yg(25) +  &
          rgs(74) * yg(26) +  &
          rgs(75) * yg(27)
   yg(37) = prod / dest

! --- bzo
   prod = rgs(105) * yg(21) * yg(31)
   dest = rgs(106) * yg(4) +  &
          rgs(107) * yg(25) +  &
          rgs(108)
   yg(38) = prod / dest

!---------------------------------------------------------------
!                 end steadygsc_yg  in-line
!---------------------------------------------------------------
!
!
!     evaluate function for predictor step
!
!      in-line version to replace subroutine:
!         call diffungsc(ng, nspec, yg, ydpreg, apreg, bpreg)
! steadygsc_il_yg.cdk
!
!***
!
!      in-line "yg" code version of subroutine:
!      subroutine diffungsc(ns, ntot, c, cdot, prod, dest)
!
!      called from mach_gas_driveil in the following form:
!      call diffungsc(ng, nspec, yg, ydpreg, apreg, bpreg)
!
!----------------------------------------------------------------------
!                          begin code
!----------------------------------------------------------------------

! --- check for positivity
   do i = 1, nsi
      yg(i) = max(yg(i), ymin2gs)
   end do

!----------------------------------------------------------------------
!
!....... determine chemical transformation rates .......
!
!...... compute chemical derivatives in chemeq format
!       ie.  dc/dt = cdot = prod - dest*c
!
!----------------------------------------------------------------------

!     so2
   apreg(01) = 0.0
   bpreg(01) = rgs(43) * yg(34) +  &
               rgs(92) * yg(39) +  &
               rgs(93) * yg(40)

! --- so4
   apreg(02) = bpreg(01) * yg(1)
   bpreg(02) = 0.0

! --- no
   apreg(03) = (rgs(1) +  &
                rgs(3) * yg(30) +  &
                rgs(12) * yg(31)) * yg(4) +  &
                rgs(13) * yg(31) +  &
                rgs(20) * yg(22)
   bpreg(03) = rgs(5) * yg(5) +  &
               rgs(7) * yg(31) +  &
               2.00 * rgs(8) * yg(3) * yg(44) +  &
               rgs(19) * yg(34) +  &
               rgs(27) * yg(25) +  &
               rgs(57) * yg(27) +  &
               rgs(72) * yg(37) +  &
               rgs(76) * yg(36) +  &
               rgs(80) * yg(35)

! --- no2
   apreg(04) = (rgs(5) * yg(5) +  &
          2.00 * (rgs(7) * yg(31) +  &
                rgs(8) * yg(3) * yg(44)) +  &
                rgs(27) * yg(25) +  &
                rgs(57) * yg(27) +  &
                rgs(76) * yg(36) +  &
                rgs(80) * yg(35)) * yg(3) +  &
                rgs(10) * yg(32) +  &
                rgs(14) * yg(31) +  &
                rgs(23) * yg(7) +  &
               (rgs(29) +  &
                rgs(30) +  &
                rgs(31) * yg(34)) * yg(33) +  &
                rgs(61) * yg(8) +  &
                rgs(71) * yg(23) * yg(34) +  &
               (rgs(87) * yg(11) +  &
                rgs(91) * yg(12) +  &
                rgs(112) * yg(24)) * yg(31)
   bpreg(04) = rgs(1) +  &
              (rgs(3) +  &
               rgs(4)) * yg(30) +  &
               rgs(6) * yg(5) +  &
               rgs(9) * yg(31) +  &
               2.00 * rgs(21) * yg(43) +  &
               rgs(22) * yg(34) +  &
               rgs(28) * yg(25) +  &
               rgs(58) * yg(27) +  &
               rgs(106) * yg(38)

! --- o3
   apreg(05) = rgs(2) * yg(30) * yg(44) * yg(45)
   bpreg(05) = rgs(5) * yg(3) +  &
               rgs(6) * yg(4) +  &
               rgs(15) +  &
               rgs(16) +  &
               rgs(26) * yg(34) +  &
               rgs(32) * yg(25) +  &
               rgs(85) * yg(11) +  &
               rgs(89) * yg(12) +  &
               rgs(110) * yg(24)

! --- h2o2
   apreg(06) = (rgs(33) +  &
                rgs(34) * yg(45) +  &
               (rgs(35) +  &
                rgs(36)) * yg(43)) * yg(25) ** 2
   bpreg(06) = rgs(37) + rgs(38) * yg(34)

! --- hno3
   apreg(07) = (2.00 * rgs(11) * yg(32) +  &
                rgs(21) * yg(4)) * yg(43) +  &
                rgs(22) * yg(4) * yg(34) +  &
              ((rgs(39) + rgs(40) * yg(45)) * yg(25) +  &
               (rgs(41) + rgs(42)) * yg(25) * yg(43) +  &
                rgs(52) * yg(15) +  &
                rgs(56) * yg(16) +  &
                rgs(66) * yg(18) +  &
                rgs(105) * yg(21)) * yg(31)
   bpreg(07) = rgs(23) + rgs(24) * yg(34)

! --- pan
   apreg(08) = rgs(58) * yg(4) * yg(27)
   bpreg(08) = rgs(61)

! --- c3h8
   apreg(09) = 0.0
   bpreg(09) = rgs(69) * yg(34)

! --- alka
   apreg(10) = 0.0
   bpreg(10) = rgs(70) * yg(34)

! --- ethe
   apreg(11) = (rgs(109) * yg(34) +  &
                0.50 * rgs(110) * yg(5) +  &
!     *             rgs(111)*yg(29))*yg(24)
                rgs(111) * yg(30)) * yg(24)
   bpreg(11) = rgs(84) * yg(34) +  &
               rgs(85) * yg(5) +  &
               rgs(86) * yg(30) +  &
               rgs(87) * yg(31)

! --- alke
   apreg(12) = 0.0
   bpreg(12) = rgs(88) * yg(34) +  &
               rgs(89) * yg(5) +  &
               rgs(90) * yg(30) +  &
               rgs(91) * yg(31)

! --- tolu
   apreg(13) = 0.0
   bpreg(13) = rgs(100) * yg(34)

! --- arom
   apreg(14) = 0.0
   bpreg(14) = rgs(101) * yg(34)

! --- hcho
   apreg(15) = rgs(55) * yg(16) +  &
              (rgs(57) * yg(3) +  &
               rgs(59) * yg(25) +  &
               2.00 * rgs(60) * yg(27) +  &
               rgs(75) * yg(37) +  &
               rgs(79) * yg(36) +  &
               rgs(83) * yg(35)) * yg(27) +  &
              (0.50 * rgs(63) * yg(17) +  &
               rgs(67) * yg(41) +  &
               bgs(01) * rgs(70) * yg(10) +  &
               0.16 * rgs(71) * yg(23) +  &
               1.56 * rgs(84) * yg(11) +  &
               0.11 * rgs(100) * yg(13) +  &
               bgs(24) * rgs(101) * yg(14)) * yg(34) +  &
              (rgs(85) * yg(5) +  &
               rgs(86) * yg(30) +  &
               2.00 * rgs(87) * yg(31)) * yg(11) +  &
              (bgs(08) * rgs(88) * yg(34) +  &
               bgs(10) * rgs(89) * yg(5) +  &
               bgs(18) * rgs(90) * yg(30) +  &
               bgs(08) * rgs(91) * yg(31)) * yg(12) +  &
               rgs(92) * yg(1) * yg(39) +  &
              (rgs(109) * yg(34) +  &
               rgs(110) * yg(5) +  &
               rgs(112) * yg(31)) * yg(24)
   bpreg(15) = rgs(49) +  &
               rgs(50) +  &
               rgs(51) * yg(34) +  &
               rgs(52) * yg(31) +  &
               rgs(53) * yg(25) +  &
               rgs(96) * yg(39) +  &
               rgs(97) * yg(40)

! --- ald2
   apreg(16) = rgs(62) * yg(17) +  &
              (0.50 * rgs(63) * yg(17) +  &
               rgs(68) * yg(42) +  &
               0.30 * rgs(69) * yg(9) +  &
               bgs(02) * rgs(70) * yg(10) +  &
               1.52 * rgs(71) * yg(23) +  &
               0.22 * rgs(84) * yg(11) +  &
               bgs(09) * rgs(88) * yg(12)) * yg(34) +  &
              (bgs(11) * rgs(89) * yg(5) +  &
               bgs(19) * rgs(90) * yg(30) +  &
               bgs(09) * rgs(91) * yg(31)) * yg(12) +  &
               rgs(93) * yg(1) * yg(40) +  &
              (0.20 * rgs(109) * yg(34) +  &
               0.40 * rgs(110) * yg(5) +  &
               rgs(111) * yg(30) +  &
               rgs(112) * yg(31)) * yg(24)
   bpreg(16) = rgs(54) * yg(34) +  &
               rgs(55) +  &
               rgs(56) * yg(31) +  &
               rgs(98) * yg(39) +  &
               rgs(99) * yg(40)

! --- mek
   apreg(17) = (0.50 * rgs(69) * yg(9) +  &
                bgs(03) * rgs(70) * yg(10) +  &
                0.14 * rgs(71) * yg(23)) * yg(34) +  &
               (rgs(73) * yg(25) +  &
                rgs(74) * yg(26) +  &
                rgs(75) * yg(27)) * yg(37) +  &
                bgs(17) * rgs(90) * yg(12) * yg(30)
   bpreg(17) = rgs(62) + rgs(63) * yg(34)

! --- mgly
   apreg(18) = (0.13 * rgs(100) * yg(13) +  &
                bgs(23) * rgs(101) * yg(14) +  &
                0.20 * rgs(104) * yg(21)) * yg(34) +  &
               (0.27 * rgs(109) * yg(34) +  &
                0.20 * rgs(110) * yg(5)) * yg(24)
   bpreg(18) = rgs(64) +  &
                rgs(65) * yg(34) +  &
                rgs(66) * yg(31)

! --- dial
   apreg(19) = (0.40 * rgs(100) * yg(13) +  &
                bgs(22) * rgs(101) * yg(14)) * yg(34)
   bpreg(19) = rgs(102) * yg(34) + rgs(103)

! --- rooh
   apreg(20) = (rgs(59) * yg(27) +  &
                rgs(73) * yg(37) +  &
                rgs(77) * yg(36) +  &
                rgs(81) * yg(35)) * yg(25)
   bpreg(20) = rgs(48) + rgs(114) * yg(34)

! --- cres
   apreg(21) = (0.16 * rgs(100) * yg(13) +  &
                0.17 * rgs(101) * yg(14)) * yg(34)
   bpreg(21) = 0.92 * rgs(104) * yg(34) +  &
                0.50 * rgs(105) * yg(31)

! --- hono
   apreg(22) = rgs(19) * yg(3) * yg(34) +  &
               rgs(21) * yg(4) * yg(43)
   bpreg(22) = rgs(20)

! --- rno3
   apreg(23) = rgs(72) * yg(3) * yg(37) +  &
               rgs(106) * yg(4) * yg(38)
   bpreg(23) = rgs(71) * yg(34)

! --- isop
   apreg(24) = 0.0
   bpreg(24) = rgs(109) * yg(34) +  &
               rgs(110) * yg(5) +  &
               rgs(111) * yg(30) +  &
               rgs(112) * yg(31)

! --- ho2
   apreg(25) = (rgs(25) * yg(28) +  &
                rgs(26) * yg(5) +  &
                rgs(38) * yg(6) +  &
                rgs(43) * yg(1) +  &
                0.16 * rgs(100) * yg(13) +  &
                0.17 * rgs(101) * yg(14)) * yg(34) +  &
               (rgs(29) + rgs(30)) * yg(33) +  &
                rgs(48) * yg(20) +  &
               (2.00 * rgs(49) +  &
                rgs(51) * yg(34) +  &
                rgs(52) * yg(31)) * yg(15) +  &
                rgs(55) * yg(16) +  &
               (2.00 * rgs(60) * yg(27) +  &
                rgs(75) * yg(37) +  &
                rgs(79) * yg(36) +  &
                rgs(83) * yg(35)) * yg(27) +  &
                rgs(64) * yg(18) +  &
                0.50 * rgs(74) * yg(26) * yg(37) +  &
               (rgs(80) * yg(3) +  &
                0.50 * rgs(82) * yg(26)) * yg(35) +  &
               (0.12 * rgs(85) * yg(5) +  &
                rgs(86) * yg(30)) * yg(11) +  &
               (bgs(13) * rgs(89) * yg(5) +  &
                bgs(20) * rgs(90) * yg(30)) * yg(12) +  &
                rgs(103) * yg(19) +  &
               (0.70 * rgs(109) * yg(34) +  &
                0.40 * rgs(110) * yg(5) +  &
                0.60 * rgs(111) * yg(30)) * yg(24)
   bpreg(25) = rgs(27) * yg(3) +  &
                rgs(28) * yg(4) +  &
                rgs(32) * yg(5) +  &
                2.00 * yg(25) * (rgs(33) +  &
                rgs(34) * yg(45) +  &
               (rgs(35) + rgs(36)) * yg(43)) +  &
               (rgs(39) + rgs(40) * yg(45) +  &
               (rgs(41) + rgs(42)) * yg(43)) * yg(31) +  &
                rgs(53) * yg(15) +  &
                rgs(59) * yg(27) +  &
                rgs(73) * yg(37) +  &
                rgs(77) * yg(36) +  &
                rgs(81) * yg(35) +  &
                rgs(107) * yg(38) +  &
                rgs(113) * yg(34)

! --- ro2
   apreg(26) = rgs(53) * yg(15) * yg(25) +  &
               rgs(55) * yg(16) +  &
               rgs(57) * yg(3) * yg(27) +  &
              (rgs(62) +  &
               1.50 * rgs(63) * yg(34)) * yg(17) +  &
              (rgs(67) * yg(41) +  &
               rgs(68) * yg(42) +  &
               rgs(69) * yg(9) +  &
               bgs(07) * rgs(70) * yg(10) +  &
               1.39 * rgs(71) * yg(23) +  &
               rgs(84) * yg(11) +  &
               0.84 * rgs(100) * yg(13) +  &
               0.83 * rgs(101) * yg(14) +  &
               rgs(104) * yg(21)) * yg(34) +  &
              (rgs(86) * yg(30) +  &
               rgs(87) * yg(31)) * yg(11) +  &
              (rgs(88) * yg(34) +  &
               bgs(12) * rgs(89) * yg(5) +  &
               bgs(21) * rgs(90) * yg(30) +  &
               rgs(91) * yg(31)) * yg(12) +  &
              (rgs(109) * yg(34) +  &
               0.50 * rgs(111) * yg(30) +  &
               rgs(112) * yg(31)) * yg(24) +  &
               0.50 * rgs(114) * yg(20) * yg(34)
   bpreg(26) = rgs(44) * yg(3) +  &
               rgs(45) * yg(25) +  &
               2.00 * rgs(46) * yg(26) +  &
               rgs(47) * yg(27)

! --- mco3
   apreg(27) = (rgs(54) * yg(34) +  &
                rgs(56) * yg(31)) * yg(16) +  &
                rgs(61) * yg(8) +  &
               (rgs(62) +  &
                rgs(63) * yg(34)) * yg(17) +  &
               (rgs(64) +  &
                rgs(65) * yg(34) +  &
                rgs(66) * yg(31)) * yg(18) +  &
               (rgs(102) * yg(34) +  &
                rgs(103)) * yg(19) +  &
                0.20 * rgs(109) * yg(24) * yg(34)
   bpreg(27) = rgs(57) * yg(3) +  &
               rgs(58) * yg(4) +  &
               rgs(59) * yg(25) +  &
               2.00 * rgs(60) * yg(27) +  &
               rgs(75) * yg(37) +  &
               rgs(79) * yg(36) +  &
               rgs(83) * yg(35)

! --- co
    apreg(28) = rgs(49) * yg(15) +  &
                 rgs(50) * yg(15) +  &
                 rgs(51) * yg(15) * yg(34) +  &
                 rgs(52) * yg(15) * yg(31) +  &
                 rgs(55) * yg(16) +  &
                 rgs(64) * yg(18) +  &
                 rgs(65) * yg(18) * yg(34) +  &
                 rgs(66) * yg(18) * yg(31) +  &
                 rgs(85) * yg(11) * yg(5) * 0.42 +  &
                 rgs(86) * yg(11) * yg(30) +  &
                 rgs(89) * yg(12) * yg(5) * bgs(15) +  &
                 rgs(90) * yg(12) * yg(30) * bgs(16) +  &
                 rgs(100) * yg(13) * yg(34) * 0.11 +  &
                 rgs(101) * yg(14) * yg(34) * bgs(24) +  &
                 rgs(103) * yg(19)
    bpreg(28) = rgs(25) * yg(34)

! --- compute net derivative
   do i = 1, nsi
      ydpreg(i) = apreg(i) - bpreg(i) * yg(i)
   end do

!---------------------------------------------------------------
!                 end diffungsc_yg  in-line
!---------------------------------------------------------------

   end if

!     test for stiffness every nstchk steps
   if (nsteps == 0 .or. mod(nsteps, nstchk) == 0) then
      nnormg = 0
      nstifg = 0
      xcrit = 0.9 / dts
      do i = 1, nsi
         if (abs(bpreg(i)) > xcrit) then
            nstifg = nstifg + 1
            jst(nstifg) = i
         else
            nnormg = nnormg + 1
            jnst(nnormg) = i
         end if
      end do
   end if

!     increment step counter
   nsteps = nsteps + 1
   if (nsteps > nstmax) then
      indx = -10
      return
   end if

! integrate with predictor eqns (compute ypreg from yg)
! species loop
 23   continue
   dum2 = dts + smf

!     non-stiff species
   do j = 1, nnormg
      i = jnst(j)
      ypreg(i) = yg(i) + dts * ydpreg(i)
   end do

!     stiff species
   do j = 1, nstifg
      i = jst(j)
      xtemp = 2.0 / (bpreg(i) + smf)

      xtemp2 = xtemp + dum2
      if (xtemp2 == 0.0) then
         divxtemp2 = max(abs(xtemp2),smf)
         xtemp2 = sign(divxtemp2,xtemp2)
      end if

      ypreg(i)=(yg(i)*(xtemp-dts)+xtemp*dts*apreg(i))/(xtemp2)
   end do

!     integrate with corrector step (compute yc2 from yc1 (or ypreg)) --
   do i = 1, nsi
      yc1g(i) = ypreg(i)
   end do
   do i = np1g, nspec
      yc1g(i) = yg(i)
   end do
   dthalf = 0.5 * dts

!     iteration loop
   do k = 1, niter
      ersmsg = 0.0
      ersmng = 0.0

! --- check positivity for yc1g(1:nsi)
      do i = 1, nsi
         yc1g(i) = max(yc1g(i), ymin2gs)
      end do
!
!     evaluate function for each corrector iteration
!
!      in-line version to replace subroutines:
!         call steadygsc(yc1g, ig, jg)
!         call diffungsc(ng, nspec, yc1g, ydc1g, ac1g, bc1g)
! steadygsc_il_yc1g.cdk
!
!***
!
!      in-line code version of suborutine:
!
!      subroutine steadygsc(c, igrid, jgrid)
!
!      called from mach_gas_driveil in the following form:
!      call steadygsc(yc1g, ig, jg)
!

! --- o(1d)
      prod = rgs(16) * yc1g(5)
      dest = rgs(17) * yc1g(43) + rgs(18) * yc1g(45)
      yc1g(29) = prod / dest

! --- crgs1
      prod = (0.40 * rgs(85) * yc1g(11) +  &
            bgs(25) * rgs(89) * yc1g(12)) * yc1g(5)
      dest = rgs(92) * yc1g(1) +  &
            rgs(94) * yc1g(43) +  &
            rgs(96) * yc1g(15) +  &
            rgs(98) * yc1g(16)
      yc1g(39) = prod / dest

! --- crgs2
      prod = (bgs(26) * rgs(89) * yc1g(12) +  &
            0.20 * rgs(110) * yc1g(24)) * yc1g(5)
      dest = rgs(93) * yc1g(1) +  &
            rgs(95) * yc1g(43) +  &
            rgs(97) * yc1g(15) +  &
            rgs(99) * yc1g(16)
      yc1g(40) = prod / dest

! --- use the max of last concentration or 1.e-29 ppm for initial
! --- guess of [no3]
      yc1g(31) = max(yc1g(31), 1.0e-29)

! --- o(3p)
      prod = rgs(1) * yc1g(4) +  &
            rgs(14) * yc1g(31) +  &
            rgs(15) * yc1g(5) +  &
            rgs(18) * yc1g(29) * yc1g(45)
      dest = rgs(2) * yc1g(44) * yc1g(45) +  &
            (rgs(3) + rgs(4)) * yc1g(4) +  &
            rgs(86) * yc1g(11) +  &
            rgs(90) * yc1g(12) +  &
            rgs(111) * yc1g(24)
      yc1g(30) = prod / dest

! --- use the max of last concentration or 1.e-29 ppm for initial
! --- guess of [hno4]
      yc1g(33) = max(yc1g(33), 1.e-29)

! --- oh
      prod = 2.00 * rgs(17) * yc1g(29) * yc1g(43) +  &
            rgs(20) * yc1g(22) +  &
            rgs(23) * yc1g(7) +  &
            (rgs(27) * yc1g(3) +  &
            rgs(32) * yc1g(5)) * yc1g(25) +  &
            2.00 * rgs(37) * yc1g(6) +  &
            rgs(48) * yc1g(20) +  &
            (bgs(14) * rgs(89) * yc1g(12) +  &
            0.10 * rgs(110) * yc1g(24)) * yc1g(5) +  &
            0.50 * rgs(114) * yc1g(20) * yc1g(34)
      dest = rgs(19) * yc1g(3) +  &
            rgs(22) * yc1g(4) +  &
            rgs(24) * yc1g(7) +  &
            rgs(25) * yc1g(28) +  &
            rgs(26) * yc1g(5) +  &
            rgs(31) * yc1g(33) +  &
            rgs(38) * yc1g(6) +  &
            rgs(43) * yc1g(1) +  &
            rgs(51) * yc1g(15) +  &
            rgs(54) * yc1g(16) +  &
            rgs(63) * yc1g(17) +  &
            rgs(65) * yc1g(18) +  &
            rgs(67) * yc1g(41) +  &
            rgs(68) * yc1g(42) +  &
            rgs(69) * yc1g(9) +  &
            rgs(70) * yc1g(10) +  &
            rgs(71) * yc1g(23) +  &
            rgs(84) * yc1g(11) +  &
            rgs(88) * yc1g(12) +  &
            rgs(100) * yc1g(13) +  &
            rgs(101) * yc1g(14) +  &
            rgs(102) * yc1g(19) +  &
            rgs(104) * yc1g(21) +  &
            rgs(109) * yc1g(24) +  &
            rgs(113) * yc1g(25) +  &
            rgs(114) * yc1g(20)
      yc1g(34) = prod / dest

! --- hno4
      prod = rgs(28) * yc1g(4) * yc1g(25)
      dest = rgs(29) + rgs(30) + rgs(31) * yc1g(34)
      yc1g(33) = prod / dest

! --- solve no3/n2o5 system
! --- no3
      rhs = -((rgs(4) * yc1g(30) +  &
                  rgs(6) * yc1g(5)) * yc1g(4) +  &
                  rgs(24) * yc1g(7) * yc1g(34))
      a(1, 2) = rgs(10)
      a(1, 1) = -(rgs(7) * yc1g(3) +  &
               (rgs(9) + rgs(12)) * yc1g(4) +  &
               rgs(13) + rgs(14) +  &
               (rgs(39) + rgs(40) * yc1g(45) +  &
               (rgs(41) + rgs(42)) * yc1g(43)) * yc1g(25) +  &
               rgs(52) * yc1g(15) +  &
               rgs(56) * yc1g(16) +  &
               rgs(66) * yc1g(18) +  &
               rgs(87) * yc1g(11) +  &
               rgs(91) * yc1g(12) +  &
               rgs(105) * yc1g(21) +  &
               rgs(112) * yc1g(24))

! --- n2o5
      a(2, 1) = rgs(9) * yc1g(4)
      a(2, 2) = -(rgs(10) + rgs(11) * yc1g(43))
      det = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
      if (det == 0.0) then
         write(0, *)  '### Error in mach_mach_gas_driveil_omp ###'
         write(0, *)  '# subroutine steadyg was unable to solve'
         write(0, *)  '# the equations for no3 and n2o5'
         write(0, *)  '# det=0 in steadyg, i=', ig, ' j=', jg
         write(0, *)  '# a(1, 1)   =', a(1, 1)
         write(0, *)  '# a(1, 2)   =', a(1, 2)
         write(0, *)  '# a(2, 1)   =', a(2, 1)
         write(0, *)  '# a(2, 2)   =', a(2, 2)
         write(0, *)  '# rg( )    =', rgs
         write(0, *)  '# yc1g ( ) =', yc1g
         write(0, *)  '###             Abort            ###'

!        only in in-line code
!        note: return is from drive to gaschem, not steadyg to drive!!!
         indx = -22
         return

      end if
      deti = 1.0 / det
      yc1g(31) = rhs * a(2, 2) * deti
      yc1g(32) = -rhs * a(2, 1) * deti

! --- ro2r
      prod = rgs(53) * yc1g(15) * yc1g(25) +  &
            rgs(55) * yc1g(16) +  &
            rgs(57) * yc1g(3) * yc1g(27) +  &
            rgs(62) * yc1g(17) +  &
            (1.50 * rgs(63) * yc1g(17) +  &
            rgs(67) * yc1g(41) +  &
            rgs(68) * yc1g(42) +  &
            rgs(69) * yc1g(9) +  &
            bgs(05) * rgs(70) * yc1g(10) +  &
            rgs(84) * yc1g(11) +  &
            rgs(88) * yc1g(12) +  &
            0.84 * rgs(100) * yc1g(13) +  &
            0.83 * rgs(101) * yc1g(14) +  &
            0.85 * rgs(104) * yc1g(21)) * yc1g(34) +  &
            (rgs(86) * yc1g(11) +  &
            bgs(21) * rgs(90) * yc1g(12)) * yc1g(30) +  &
            bgs(12) * rgs(89) * yc1g(5) * yc1g(12) +  &
            0.50 * rgs(114) * yc1g(20) * yc1g(34)
      dest = rgs(80) * yc1g(3) +  &
            rgs(81) * yc1g(25) +  &
            rgs(82) * yc1g(26) +  &
            rgs(83) * yc1g(27)
      yc1g(35) = prod / dest

! --- r2o2
      prod = (bgs(06) * rgs(70) * yc1g(10) +  &
            1.39 * rgs(71) * yc1g(23)) * yc1g(34) +  &
            (rgs(87) * yc1g(11) +  &
            rgs(91) * yc1g(12)) * yc1g(31) +  &
            (0.90 * rgs(109) * yc1g(34) +  &
            0.50 * rgs(111) * yc1g(30) +  &
            rgs(112) * yc1g(31)) * yc1g(24)
      dest = rgs(76) * yc1g(3) +  &
            rgs(77) * yc1g(25) +  &
            rgs(78) * yc1g(26) +  &
            rgs(79) * yc1g(27)
      yc1g(36) = prod / dest

! --- ro2n
      prod = (bgs(04) * rgs(70) * yc1g(10) +  &
            0.14 * rgs(104) * yc1g(21) +  &
            0.10 * rgs(109) * yc1g(24)) * yc1g(34)
      dest = rgs(72) * yc1g(3) +  &
            rgs(73) * yc1g(25) +  &
            rgs(74) * yc1g(26) +  &
            rgs(75) * yc1g(27)
      yc1g(37) = prod / dest

! --- bzo
      prod = rgs(105) * yc1g(21) * yc1g(31)
      dest = rgs(106) * yc1g(4) +  &
            rgs(107) * yc1g(25) +  &
            rgs(108)
      yc1g(38) = prod / dest

!---------------------------------------------------------------
!                 end steadygsc_yc1g  in-line
!---------------------------------------------------------------
!
! steadygsc_il_yc1g.cdk
!
!***
!
!      in-line "yc1g" code version of suborutine:
!      subroutine diffungsc(ns, ntot, c, cdot, prod, dest)
!
!      called from mach_gas_driveil in the following form:
!      call diffungsc(ng, nspec, yc1g, ydc1g, ac1g, bc1g)
!
!----------------------------------------------------------------------
!                          begin code
!----------------------------------------------------------------------
!
!....... determine chemical transformation rates .......
!
!...... compute chemical derivatives in chemeq format
!       ie.  dc/dt = cdot = prod - dest*c
!
!----------------------------------------------------------------------

!     so2
      ac1g(01) = 0.0
      bc1g(01) = rgs(43) * yc1g(34) +  &
               rgs(92) * yc1g(39) +  &
               rgs(93) * yc1g(40)

! --- so4
      ac1g(02) = bc1g(01) * yc1g(1)
      bc1g(02) = 0.0

! --- no
      ac1g(03) = (rgs(1) +  &
                  rgs(3) * yc1g(30) +  &
                  rgs(12) * yc1g(31)) * yc1g(4) +  &
                  rgs(13) * yc1g(31) +  &
                  rgs(20) * yc1g(22)
      bc1g(03) = rgs(5) * yc1g(5) +  &
               rgs(7) * yc1g(31) +  &
               2.00 * rgs(8) * yc1g(3) * yc1g(44) +  &
               rgs(19) * yc1g(34) +  &
               rgs(27) * yc1g(25) +  &
               rgs(57) * yc1g(27) +  &
               rgs(72) * yc1g(37) +  &
               rgs(76) * yc1g(36) +  &
               rgs(80) * yc1g(35)

! --- no2
      ac1g(04) = (rgs(5) * yc1g(5) +  &
            2.00 * (rgs(7) * yc1g(31) +  &
                  rgs(8) * yc1g(3) * yc1g(44)) +  &
                  rgs(27) * yc1g(25) +  &
                  rgs(57) * yc1g(27) +  &
                  rgs(76) * yc1g(36) +  &
                  rgs(80) * yc1g(35)) * yc1g(3) +  &
                  rgs(10) * yc1g(32) +  &
                  rgs(14) * yc1g(31) +  &
                  rgs(23) * yc1g(7) +  &
               (rgs(29) +  &
                  rgs(30) +  &
                  rgs(31) * yc1g(34)) * yc1g(33) +  &
                  rgs(61) * yc1g(8) +  &
                  rgs(71) * yc1g(23) * yc1g(34) +  &
               (rgs(87) * yc1g(11) +  &
                  rgs(91) * yc1g(12) +  &
                  rgs(112) * yc1g(24)) * yc1g(31)
      bc1g(04) = rgs(1) +  &
               (rgs(3) +  &
               rgs(4)) * yc1g(30) +  &
               rgs(6) * yc1g(5) +  &
               rgs(9) * yc1g(31) +  &
               2.00 * rgs(21) * yc1g(43) +  &
               rgs(22) * yc1g(34) +  &
               rgs(28) * yc1g(25) +  &
               rgs(58) * yc1g(27) +  &
               rgs(106) * yc1g(38)

! --- o3
      ac1g(05) = rgs(2) * yc1g(30) * yc1g(44) * yc1g(45)
      bc1g(05) = rgs(5) * yc1g(3) +  &
               rgs(6) * yc1g(4) +  &
               rgs(15) +  &
               rgs(16) +  &
               rgs(26) * yc1g(34) +  &
               rgs(32) * yc1g(25) +  &
               rgs(85) * yc1g(11) +  &
               rgs(89) * yc1g(12) +  &
               rgs(110) * yc1g(24)

! --- h2o2
      ac1g(06) = (rgs(33) +  &
                  rgs(34) * yc1g(45) +  &
               (rgs(35) +  &
                  rgs(36)) * yc1g(43)) * yc1g(25) ** 2
      bc1g(06) = rgs(37) + rgs(38) * yc1g(34)

! --- hno3
      ac1g(07) = (2.00 * rgs(11) * yc1g(32) +  &
                  rgs(21) * yc1g(4)) * yc1g(43) +  &
                  rgs(22) * yc1g(4) * yc1g(34) +  &
               ((rgs(39) + rgs(40) * yc1g(45)) * yc1g(25) +  &
               (rgs(41) + rgs(42)) * yc1g(25) * yc1g(43) +  &
                  rgs(52) * yc1g(15) +  &
                  rgs(56) * yc1g(16) +  &
                  rgs(66) * yc1g(18) +  &
                  rgs(105) * yc1g(21)) * yc1g(31)
      bc1g(07) = rgs(23) + rgs(24) * yc1g(34)

! --- pan
      ac1g(08) = rgs(58) * yc1g(4) * yc1g(27)
      bc1g(08) = rgs(61)

! --- c3h8
      ac1g(09) = 0.0
      bc1g(09) = rgs(69) * yc1g(34)

! --- alka
      ac1g(10) = 0.0
      bc1g(10) = rgs(70) * yc1g(34)

! --- ethe
      ac1g(11) = (rgs(109) * yc1g(34) +  &
                  0.50 * rgs(110) * yc1g(5) +  &
                  rgs(111) * yc1g(30)) * yc1g(24)
      bc1g(11) = rgs(84) * yc1g(34) +  &
               rgs(85) * yc1g(5) +  &
               rgs(86) * yc1g(30) +  &
               rgs(87) * yc1g(31)

! --- alke
      ac1g(12) = 0.0
      bc1g(12) = rgs(88) * yc1g(34) +  &
               rgs(89) * yc1g(5) +  &
               rgs(90) * yc1g(30) +  &
               rgs(91) * yc1g(31)

! --- tolu
      ac1g(13) = 0.0
      bc1g(13) = rgs(100) * yc1g(34)

! --- arom
      ac1g(14) = 0.0
      bc1g(14) = rgs(101) * yc1g(34)

! --- hcho
      ac1g(15) = rgs(55) * yc1g(16) +  &
               (rgs(57) * yc1g(3) +  &
               rgs(59) * yc1g(25) +  &
               2.00 * rgs(60) * yc1g(27) +  &
               rgs(75) * yc1g(37) +  &
               rgs(79) * yc1g(36) +  &
               rgs(83) * yc1g(35)) * yc1g(27) +  &
               (0.50 * rgs(63) * yc1g(17) +  &
               rgs(67) * yc1g(41) +  &
               bgs(01) * rgs(70) * yc1g(10) +  &
               0.16 * rgs(71) * yc1g(23) +  &
               1.56 * rgs(84) * yc1g(11) +  &
               0.11 * rgs(100) * yc1g(13) +  &
               bgs(24) * rgs(101) * yc1g(14)) * yc1g(34) +  &
               (rgs(85) * yc1g(5) +  &
               rgs(86) * yc1g(30) +  &
               2.00 * rgs(87) * yc1g(31)) * yc1g(11) +  &
               (bgs(08) * rgs(88) * yc1g(34) +  &
               bgs(10) * rgs(89) * yc1g(5) +  &
               bgs(18) * rgs(90) * yc1g(30) +  &
               bgs(08) * rgs(91) * yc1g(31)) * yc1g(12) +  &
               rgs(92) * yc1g(1) * yc1g(39) +  &
               (rgs(109) * yc1g(34) +  &
               rgs(110) * yc1g(5) +  &
               rgs(112) * yc1g(31)) * yc1g(24)
      bc1g(15) = rgs(49) +  &
               rgs(50) +  &
               rgs(51) * yc1g(34) +  &
               rgs(52) * yc1g(31) +  &
               rgs(53) * yc1g(25) +  &
               rgs(96) * yc1g(39) +  &
               rgs(97) * yc1g(40)

! --- ald2
      ac1g(16) = rgs(62) * yc1g(17) +  &
                  (0.50 * rgs(63) * yc1g(17) +  &
                  rgs(68) * yc1g(42) +  &
                  0.30 * rgs(69) * yc1g(9) +  &
                  bgs(02) * rgs(70) * yc1g(10) +  &
                  1.52 * rgs(71) * yc1g(23) +  &
                  0.22 * rgs(84) * yc1g(11) +  &
                  bgs(09) * rgs(88) * yc1g(12)) * yc1g(34) +  &
                  (bgs(11) * rgs(89) * yc1g(5) +  &
                  bgs(19) * rgs(90) * yc1g(30) +  &
                  bgs(09) * rgs(91) * yc1g(31)) * yc1g(12) +  &
                  rgs(93) * yc1g(1) * yc1g(40) +  &
                  (0.20 * rgs(109) * yc1g(34) +  &
                  0.40 * rgs(110) * yc1g(5) +  &
                  rgs(111) * yc1g(30) +  &
                  rgs(112) * yc1g(31)) * yc1g(24)
      bc1g(16) = rgs(54) * yc1g(34) +  &
                  rgs(55) +  &
                  rgs(56) * yc1g(31) +  &
                  rgs(98) * yc1g(39) +  &
                  rgs(99) * yc1g(40)

! --- mek
      ac1g(17) = (0.50 * rgs(69) * yc1g(9) +  &
                  bgs(03) * rgs(70) * yc1g(10) +  &
                  0.14 * rgs(71) * yc1g(23)) * yc1g(34) +  &
               (rgs(73) * yc1g(25) +  &
                  rgs(74) * yc1g(26) +  &
                  rgs(75) * yc1g(27)) * yc1g(37) +  &
                  bgs(17) * rgs(90) * yc1g(12) * yc1g(30)
      bc1g(17) = rgs(62) + rgs(63) * yc1g(34)

! --- mgly
      ac1g(18) = (0.13 * rgs(100) * yc1g(13) +  &
                  bgs(23) * rgs(101) * yc1g(14) +  &
                  0.20 * rgs(104) * yc1g(21)) * yc1g(34) +  &
               (0.27 * rgs(109) * yc1g(34) +  &
                  0.20 * rgs(110) * yc1g(5)) * yc1g(24)
      bc1g(18) = rgs(64) +  &
               rgs(65) * yc1g(34) +  &
               rgs(66) * yc1g(31)

! --- dial
      ac1g(19) = (0.40 * rgs(100) * yc1g(13) +  &
                  bgs(22) * rgs(101) * yc1g(14)) * yc1g(34)
      bc1g(19) = rgs(102) * yc1g(34) + rgs(103)

! --- rooh
      ac1g(20) = (rgs(59) * yc1g(27) +  &
                  rgs(73) * yc1g(37) +  &
                  rgs(77) * yc1g(36) +  &
                  rgs(81) * yc1g(35)) * yc1g(25)
      bc1g(20) = rgs(48) + rgs(114) * yc1g(34)

! --- cres
      ac1g(21) = (0.16 * rgs(100) * yc1g(13) +  &
                  0.17 * rgs(101) * yc1g(14)) * yc1g(34)
      bc1g(21) = 0.92 * rgs(104) * yc1g(34) +  &
               0.50 * rgs(105) * yc1g(31)

! --- hono
      ac1g(22) = rgs(19) * yc1g(3) * yc1g(34) +  &
               rgs(21) * yc1g(4) * yc1g(43)
      bc1g(22) = rgs(20)

! --- rno3
      ac1g(23) = rgs(72) * yc1g(3) * yc1g(37) +  &
               rgs(106) * yc1g(4) * yc1g(38)
      bc1g(23) = rgs(71) * yc1g(34)

! --- isop
      ac1g(24) = 0.0
      bc1g(24) = rgs(109) * yc1g(34) +  &
               rgs(110) * yc1g(5) +  &
               rgs(111) * yc1g(30) +  &
               rgs(112) * yc1g(31)

! --- ho2
      ac1g(25) = (rgs(25) * yc1g(28) +  &
                  rgs(26) * yc1g(5) +  &
                  rgs(38) * yc1g(6) +  &
                  rgs(43) * yc1g(1) +  &
                  0.16 * rgs(100) * yc1g(13) +  &
                  0.17 * rgs(101) * yc1g(14)) * yc1g(34) +  &
               (rgs(29) + rgs(30)) * yc1g(33) +  &
                  rgs(48) * yc1g(20) +  &
               (2.00 * rgs(49) +  &
                  rgs(51) * yc1g(34) +  &
                  rgs(52) * yc1g(31)) * yc1g(15) +  &
                  rgs(55) * yc1g(16) +  &
               (2.00 * rgs(60) * yc1g(27) +  &
                  rgs(75) * yc1g(37) +  &
                  rgs(79) * yc1g(36) +  &
                  rgs(83) * yc1g(35)) * yc1g(27) +  &
                  rgs(64) * yc1g(18) +  &
                  0.50 * rgs(74) * yc1g(26) * yc1g(37) +  &
               (rgs(80) * yc1g(3) +  &
                  0.50 * rgs(82) * yc1g(26)) * yc1g(35) +  &
               (0.12 * rgs(85) * yc1g(5) +  &
                  rgs(86) * yc1g(30)) * yc1g(11) +  &
               (bgs(13) * rgs(89) * yc1g(5) +  &
                  bgs(20) * rgs(90) * yc1g(30)) * yc1g(12) +  &
                  rgs(103) * yc1g(19) +  &
               (0.70 * rgs(109) * yc1g(34) +  &
                  0.40 * rgs(110) * yc1g(5) +  &
                  0.60 * rgs(111) * yc1g(30)) * yc1g(24)
      bc1g(25) = rgs(27) * yc1g(3) +  &
               rgs(28) * yc1g(4) +  &
               rgs(32) * yc1g(5) +  &
               2.00 * yc1g(25) * (rgs(33) +  &
               rgs(34) * yc1g(45) +  &
               (rgs(35) + rgs(36)) * yc1g(43)) +  &
               (rgs(39) + rgs(40) * yc1g(45) +  &
               (rgs(41) + rgs(42)) * yc1g(43)) * yc1g(31) +  &
               rgs(53) * yc1g(15) +  &
               rgs(59) * yc1g(27) +  &
               rgs(73) * yc1g(37) +  &
               rgs(77) * yc1g(36) +  &
               rgs(81) * yc1g(35) +  &
               rgs(107) * yc1g(38) +  &
               rgs(113) * yc1g(34)

! --- ro2
      ac1g(26) = rgs(53) * yc1g(15) * yc1g(25) +  &
               rgs(55) * yc1g(16) +  &
               rgs(57) * yc1g(3) * yc1g(27) +  &
               (rgs(62) +  &
               1.50 * rgs(63) * yc1g(34)) * yc1g(17) +  &
               (rgs(67) * yc1g(41) +  &
               rgs(68) * yc1g(42) +  &
               rgs(69) * yc1g(9) +  &
               bgs(07) * rgs(70) * yc1g(10) +  &
               1.39 * rgs(71) * yc1g(23) +  &
               rgs(84) * yc1g(11) +  &
               0.84 * rgs(100) * yc1g(13) +  &
               0.83 * rgs(101) * yc1g(14) +  &
               rgs(104) * yc1g(21)) * yc1g(34) +  &
               (rgs(86) * yc1g(30) +  &
               rgs(87) * yc1g(31)) * yc1g(11) +  &
               (rgs(88) * yc1g(34) +  &
               bgs(12) * rgs(89) * yc1g(5) +  &
               bgs(21) * rgs(90) * yc1g(30) +  &
               rgs(91) * yc1g(31)) * yc1g(12) +  &
               (rgs(109) * yc1g(34) +  &
               0.50 * rgs(111) * yc1g(30) +  &
               rgs(112) * yc1g(31)) * yc1g(24) +  &
               0.50 * rgs(114) * yc1g(20) * yc1g(34)
      bc1g(26) = rgs(44) * yc1g(3) +  &
               rgs(45) * yc1g(25) +  &
               2.00 * rgs(46) * yc1g(26) +  &
               rgs(47) * yc1g(27)

! -- mco3
      ac1g(27) = (rgs(54) * yc1g(34) +  &
                  rgs(56) * yc1g(31)) * yc1g(16) +  &
                  rgs(61) * yc1g(8) +  &
               (rgs(62) +  &
                  rgs(63) * yc1g(34)) * yc1g(17) +  &
               (rgs(64) +  &
                  rgs(65) * yc1g(34) +  &
                  rgs(66) * yc1g(31)) * yc1g(18) +  &
               (rgs(102) * yc1g(34) +  &
                  rgs(103)) * yc1g(19) +  &
                  0.20 * rgs(109) * yc1g(24) * yc1g(34)
      bc1g(27) = rgs(57) * yc1g(3) +  &
               rgs(58) * yc1g(4) +  &
               rgs(59) * yc1g(25) +  &
               2.00 * rgs(60) * yc1g(27) +  &
               rgs(75) * yc1g(37) +  &
               rgs(79) * yc1g(36) +  &
               rgs(83) * yc1g(35)

! --- co
      ac1g(28) = rgs(49) * yc1g(15) +  &
                  rgs(50) * yc1g(15) +  &
                  rgs(51) * yc1g(15) * yc1g(34) +  &
                  rgs(52) * yc1g(15) * yc1g(31) +  &
                  rgs(55) * yc1g(16) +  &
                  rgs(64) * yc1g(18) +  &
                  rgs(65) * yc1g(18) * yc1g(34) +  &
                  rgs(66) * yc1g(18) * yc1g(31) +  &
                  rgs(85) * yc1g(11) * yc1g(5) * 0.42 +  &
                  rgs(86) * yc1g(11) * yc1g(30) +  &
                  rgs(89) * yc1g(12) * yc1g(5) * bgs(15) +  &
                  rgs(90) * yc1g(12) * yc1g(30) * bgs(16) +  &
                  rgs(100) * yc1g(13) * yc1g(34) * 0.11 +  &
                  rgs(101) * yc1g(14) * yc1g(34) * bgs(24) +  &
                  rgs(103) * yc1g(19)

      bc1g(28) = rgs(25) * yc1g(34)

! --- compute net derivative
      do i = 1, nsi
         ydc1g(i) = ac1g(i) - bc1g(i) * yc1g(i)
      end do

!---------------------------------------------------------------
!                 end diffungsc_yc1g  in-line
!---------------------------------------------------------------

!         nfes = nfes + 1
      errng = -1.0
      errsg = -1.0

!     species loop
      do i = 1, nsi
         ers(i) = 0.0
         ern(i) = 0.0
      end do

!     non-stiff species
      do j = 1, nnormg
         i = jnst(j)
         yc2g(i) = yg(i) + dthalf * (ydc1g(i) + ydpreg(i))
      end do
      do j = 1, nnormg
         i = jnst(j)
         ern(i) = abs((yc2g(i) - yc1g(i)) /  &
                  (max(min(yc1g(i), yc2g(i)), ymings) + smf))
      end do

!     stiff species
      do j = 1, nstifg
         i = jst(j)
         xtemp = (bc1g(i) + 2.0 * smf + bpreg(i)) /  &
                  ((bc1g(i) + smf) * (bpreg(i) + smf))

         xtemp2 = xtemp + dum2
         if (xtemp2 == 0.0) then
            divxtemp2 = max(abs(xtemp2),smf)
            xtemp2 = sign(divxtemp2,xtemp2)
         end if

         yc2g(i) = (dthalf * xtemp * (ac1g(i) + apreg(i)) +  &
                  yg(i) * (xtemp - dts)) /  &
                  (xtemp2)
      end do


! sylvain test 
!	do i=1,nsi
!	if (yc2g(i) .lt.0 ) then
!	print *,'ATTENTION YC2G < 0'
!	end if
!	end do

      do j = 1, nstifg
         i = jst(j)
         ers(i) = abs((yc2g(i) - yc1g(i)) /  &
                  (max(min(yc1g(i), yc2g(i)), ymings) + smf))
      end do

!     determine average and maximum errors for non-stiff species
      do j = 1, nnormg
         i = jnst(j)
         ersmng = ersmng + ern(i)
         errng = max(errng, ern(i))
      end do

!     determine average and maximum errors for stiff species
      do j = 1, nstifg
         i = jst(j)
         ersmsg = ersmsg + ers(i)
         errsg = max(errsg, ers(i))
      end do

!     compute average error
      if (nnormg > 0) ersmng = ersmng / nnormg
      if (nstifg > 0) ersmsg = ersmsg / nstifg

!     convergence test
!     error control method:
!         max error < teps  and average error < eps
      kiter = k
      iconv = 1
      if (nnormg /= 0) then
         if (errng > tepsng) iconv = -1
         if (ersmng > epsng) iconv = -2
      end if
      if (nstifg /= 0) then
         if (errsg > tepssg) iconv = -3
         if (ersmsg > epssg) iconv = -4
      end if

!     exit iteration loop if convergence obtained
      if (iconv > 0) go to 80

!     use solution from this iteration for initial in next iteration
      do i = 1, nsi
         yc1g(i) = yc2g(i)
      end do
   end do ! k = 1, niter

!***********************************************************************
!     no convergence after niter iterations
!     select smaller time step
!     retry step beginning at predictor eqns.
   dts = max(dts * 0.50, dtmin)
   if (dts <= dtmin) go to 110
   go to 23

!     case: convergence obtained
 80   continue
   do i = 1, nsi
      yg(i) = yc2g(i)
   end do

!     check for negative concentrations
   do i = 1, nsi
      yg(i) = max(yg(i), ymin2gs)
   end do

!     increment time; exit if tout is reached
   time = time + dts
   if (time >= tout) go to 100

!     select time step for next step
   if (kiter > 2) go to 90
   if (dts >= 0.001) go to 83
   if (kiter == 1) dts = 10.0 * dts
   if (kiter == 2) dts = 3.0 * dts
   go to 90
 83   if (dts >= 0.10) go to 85
 84   if (kiter == 1) dts = 3.00 * dts
   if (kiter == 2) dts = 1.50 * dts
   go to 90
 85   if (zen > 90.0) go to 84
   if (kiter == 1) dts = 1.35 * dts
   if (kiter == 2) dts = 1.15 * dts
 90   dts = min(dts, (tout - time), dtmax)

!     return to top of integration loop for the next step
   go to 18

!     integration is completed
 100  continue
   indx = 0

!     conserve nitrogen
!
!***
!     begin cnoxsc in-line
!
!     second call to cnoxsc in-lined and
!     loop unrolled
!      call cnoxsc(yg, k2)
! paul commented out second call to cnox
!     t2 = 0.0
!     t2 = yg(3)+yg(4)+yg(7)+yg(8)+yg(22)+yg(23)
!     d      = t2 - t1
!     tnono2 = yg(3) + yg(4)
!     fno    = yg(3)/tnono2
!     fno2   = yg(4)/tnono2
!     dno    = fno*d
!     dno2   = fno2*d
!     yg(3)  = max((yg(3)-dno), 1.e-09)
!     yg(4)  = max((yg(4)-dno2), 1.e-09)
!     zo3    = yg(5) - d
!     yg(5)  = max(zo3, 1.e-06)
!
!     end cnoxsc in-line

   return
 110  continue

   write(chm_lun_out, *) 'nsi, jst, jnst ' , nsi, jst, jnst
!     integration fails
   indx = -1
   do i = 1, nsi
      yg(i) = yc2g(i)
   end do
   return
end subroutine mach_adom2_drive
