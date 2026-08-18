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
! Fichier/File   : mach_hetv_v.ftn90
! Creation       : P. Makar,  V. Bouchet,  S. Gravel,  B. Pabla,  S. Menard
! Description    : This subroutine solves the heterogeneous system for NH4-NO3-SO4,
!                  vectorizing over the gridpoint dimension and utilizing the
!                  systems of equations set out in Nenes and Pandis 1998 ISORROPIA code.
!
! Extra info     : Modified to vectorize over a 1-D array (npts dimension),
!                  reshape'd from the model 2-D array vertical slab.
!                  - (S. Gravel and A. Akingunola, August 2016)
!
!                : Data arrays for equilibrium constant calculations:
!
!                      data react /
!                   1'1:  HNO3(g_eq) <=> HNO3(aq)',
!                   2'2:  HNO3(aq) <=> H+ + NO3-',
!                   3'3:  NH3(g_eq) <=> NH3(aq)',
!                   4'4:  NH3(aq) + H2O <=> NH4+ + OH-',
!                   5'5:  HSO4- <=> SO4= + H+',
!                   6'6:  NH4HSO4(s) <=> NH4+ + HSO4-',
!                   7'7:  (NH4)2SO4(s) <=> 2 NH4+ + SO4=',
!                   8'8:  (NH4)3H(SO4)2 <=> 3 NH4+ + HSO4- + SO4=',
!                   9'9:  NH4NO3(s) <=> NH3(g_eq) + HNO3(g_eq)',
!                   1'10:  H2O <=> H+ + OH-'/
!                The Equilibrium constants are calculated using the
!                formulae described in Pilinis and Seinfeld,  1987,  from the
!                following table of thermodynamic data.  Sources:  "Lange":  Lange's
!                handbook of chemistry,  14th edition,  1992.  "CRC":  1996 CRC
!                Handbook of Chemistry and Physics.  A&P:  Ansari and Pandis,  1999.
!                "Kim (I)":  Kim et al.(I),  Aerosol Sci and Techn. 1993.
!
!                   Species          Delta Hf0     Delta Gf0      Cp0        Source
!                                  kJ/mol        kJ/mol        J/K/mol
!                HNO3(g)           -134.3         -73.94        53.34       Lange
!                HNO3(aq)          -207.36       -111.34       -86.6        Lange
!                H+                   0             0            0          Lange
!                NO3-              -207.4        -111.3        -86.6        CRC
!                NH3(g)             -45.940       -16.407       35.630      CRC
!                NH3(aq)            -80.29        -26.57        79.9(Kim I) Lange
!                NH4+              -132.51        -79.37        79.9        CRC
!                OH-               -229.99       -157.28      -148.5        Lange
!                HSO4-             -887.3        -755.9        -84.0        CRC
!                SO4=              -909.3        -744.5       -293.0        CRC
!                NH4NO3(s)         -365.56       -184.01       139.3        Lange
!                NH4HSO4(s)       -1026.96       -823.00       127.50       A&P
!                (NH4)2SO4(s)     -1180.85(A&P)  -901.90       187.49       Lange
!                (NH4)3H(SO4)2(s) -2207.00      -1730.00       315.00       A&P
!                H2O               -285.830      -237.129       75.291      Kim (I)
!
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_main_2cases(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                                 h_i, nh3_i, lwn_i, t_i, rh_i, rho_i,       &
                                 case_number)
!!if_off
   use chm_utils_mod,         only: chm_lun_out, chm_error_l
   use mach_hetv_headers_mod, only: mach_hetv_case2, mach_hetv_case12, mach_hetv_corrhno3
   use mach_hetv_mod,         only: tstd
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: npts
   real(kind=8),    intent(inout) :: so4_i      (npts)
   real(kind=8),    intent(inout) :: no3_i      (npts)
   real(kind=8),    intent(inout) :: nh4_i      (npts)
   real(kind=8),    intent(inout) :: hso4_i     (npts)
   real(kind=8),    intent(inout) :: hno3_i     (npts)
   real(kind=8),    intent(inout) :: h_i        (npts)
   real(kind=8),    intent(inout) :: nh3_i      (npts)
   real(kind=8),    intent(inout) :: lwn_i      (npts)
   real(kind=8),    intent   (in) :: t_i        (npts)
   real(kind=8),    intent   (in) :: rh_i       (npts)
   real(kind=8),    intent   (in) :: rho_i      (npts)
   real(kind=4),    intent  (out) :: case_number(npts)
!!if_off
!
! Local variables
!
   logical(kind=4),  parameter   :: diag = .false.
   integer(kind=4),  parameter   :: nsp  = 15       !species in eq rate calculations
   integer(kind=4),  parameter   :: nr   = 10        !number of eq reactions
   integer(kind=4)               :: i
   integer(kind=4)               :: ncas(npts)
   integer(kind=4)               :: n0, ne2,  ne12
   real(kind=8),  dimension(nsp) :: dh0 = (/ -134.3d0,   -207.36d0,     0.0d0,   -207.4d0,   -45.940d0,       &
                                             -80.29d0,  -132.51d0,  -229.99d0,  -887.3d0,  -909.3d0,          &
                                             -365.56d0,  -1026.96d0,  -1180.85d0,  -2207.00d0,  -285.830d0/)
   real(kind=8),  dimension(nsp) :: dg0 = (/ -73.94d0,  -111.34d0,    0.0d0,   -111.3d0,   -16.407d0,         &
                                             -26.57d0,  -79.37d0,  -157.28d0,  -755.9d0,  -744.5d0,           &
                                             -184.01d0,  -823.00d0,  -901.90d0,  -1730.00d0,  -237.129d0/)
   real(kind=8),  dimension(nsp) :: cp0 = (/ 53.34d0,  -86.6d0,    0.0d0,  -86.6d0,   35.630d0,               &
                                             79.9d0,   79.9d0,  -148.5d0,  -84.0d0,  -293.0d0,                &
                                             139.3d0,  127.50d0,  187.49d0,  315.00d0,  75.291d0/)

   real(kind=8)                       :: k0(nr), p1(nr), p2(nr)
   real(kind=8), parameter            :: r = 8.314d0

   real(kind=8)                       :: rt_inv
   real(kind=8),     dimension(npts)  :: tats, ta, ts, tn
   logical(kind=4)                    :: local_print
! 
   local_print = ((chm_lun_out > 0) .and. diag)
!
!  Start of het_v
!  (0) Calculate DG0/(RT),  DH0/(RT),  DCp0/R  for each reaction.
!  1:  HNO3(g_eq) <=> HNO3(aq)
   rt_inv = 1.0d3 / (r * tstd)
   k0(1) = exp(-(dg0(2) - dg0(1)) * rt_inv)
   p1(1) = (-(dh0(2) - dh0(1)) * rt_inv)
   p2(1) = -(cp0(2) - cp0(1)) / r
!  2:  HNO3(aq) <=> H+ + NO3-
   k0(2) = exp(-(dg0(3) + dg0(4) - dg0(2)) * rt_inv)
   p1(2) = (-(dh0(3) + dh0(4) - dh0(2)) * rt_inv)
   p2(2) = -(cp0(3) + cp0(4) - cp0(2)) / r
!  3:  NH3(g_eq) <=> NH3(aq)
   k0(3) = exp(-(dg0(6) - dg0(5)) * rt_inv)
   p1(3) = (-(dh0(6) - dh0(5)) * rt_inv)
   p2(3) = -(cp0(6) - cp0(5)) / r
!  4:  NH3(aq) + H2O <=> NH4+ + OH-
   k0(4) = exp(-(dg0(7) + dg0(8) - dg0(15) - dg0(6)) * rt_inv)
   p1(4) = (-(dh0(7) + dh0(8) - dh0(15) - dh0(6)) * rt_inv)
   p2(4) = -(cp0(7) + cp0(8) - cp0(15) - cp0(6)) / r
!  5:  HSO4- <=> SO4= + H+
   k0(5) = exp(-(dg0(10) + dg0(3) - dg0(9)) * rt_inv)
   p1(5) = (-(dh0(10) + dh0(3) - dh0(9)) * rt_inv)
   p2(5) = -(cp0(10) + cp0(3) - cp0(9)) / r
!  6:  NH4HSO4(s) <=> NH4+ + HSO4-
   k0(6) = exp(-(dg0(7) + dg0(9) - dg0(12)) * rt_inv)
   p1(6) = (-(dh0(7) + dh0(9) - dh0(12)) * rt_inv)
   p2(6) = -(cp0(7) + cp0(9) - cp0(12)) / r
!  7:  (NH4)2SO4(s) <=> 2 NH4+ + SO4=
   k0(7) = exp(-(2.0d0 * dg0(7) + dg0(10) - dg0(13)) * rt_inv)
   p1(7) = (-(2.0d0 * dh0(7) + dh0(10) - dh0(13)) * rt_inv)
   p2(7) = -(2.0d0 * cp0(7) + cp0(10) - cp0(13)) / r
!  8:  (NH4)3H(SO4)2 <=> 3 NH4+ + HSO4- + SO4=
   k0(8) = exp(-(3.0d0 * dg0(7) + dg0(9) + dg0(10) - dg0(14)) * rt_inv)
   p1(8) = (-(3.0d0 * dh0(7) + dh0(9) + dh0(10) - dh0(14)) * rt_inv)
   p2(8) = -(3.0d0 * cp0(7) + cp0(9) + cp0(10) - cp0(14)) / r
!  9:  NH4NO3(s) <=> NH3(g_eq) + HNO3(g_eq)
   k0(9) = exp(-(dg0(5) + dg0(1) - dg0(11)) * rt_inv)
   p1(9) = (-(dh0(5) + dh0(1) - dh0(11)) * rt_inv)
   p2(9) = -(cp0(5) + cp0(1) - cp0(11)) / r
! 10:  H2O <=> H+ + OH-
   k0(10) = exp(-(dg0(3) + dg0(8) - dg0(15)) * rt_inv)
   p1(10) = (-(dh0(3) + dh0(8) - dh0(15)) * rt_inv)
   p2(10) = -(cp0(3) + cp0(8) - cp0(15)) / r

   do i = 1, npts
!  Initialize output water to 0 (zero)
      lwn_i   = 0.0d0
!  Determine (total ammonia)/(total sulphate)
      ta(i)   = nh4_i(i) + nh3_i(i)
      ts(i)   = so4_i(i) + hso4_i(i)
      tn(i)   = no3_i(i) + hno3_i(i)
      if (ts(i) > 0.0d0 ) then
          tats(i) = ta(i) / ts(i)
      else
          tats(i) = 0.0d0
      end if
   end do

!
!
!  ########################################################
!   GATHER POINTS INTO CONTIGUOUS ARRAYS:  SORTING SECTION
!  #################  TA/TS < 1:  #########################
!
   n0   = 0
   ne2  = 0    !case2 counter
   ne12 = 0    !case12 counter

!           tats <  1.0,  
!     Case 2(a)
   do i = 1, npts
      if (tats(i) == 0.0d0) then
         ncas(i) = 0
         case_number(i) = 0.0
         n0 = n0 + 1
      else if (tats(i) > 0.0d0 .and. tats(i) < 1.0d0) then
         ncas(i) = 2
         case_number(i) = 2.1
         ne2  = ne2  + 1
      else if (tats(i) >= 1.0d0 .and. tats(i) < 1.5d0) then
         ncas(i) = 2
         case_number(i) = 2.2
         ne2  = ne2  + 1
      else if (tats(i) >= 1.5d0 .and. tats(i) < 2.0d0) then
         ncas(i) = 12
         case_number(i) = 12.1
         ne12 = ne12 + 1
      else if (tats(i) >= 2.0d0) then
         ncas(i) = 12
         case_number(i) = 12.2
         ne12 = ne12 + 1
      end if
   end do


!  ####################################################
!  Check:  calculate difference between sorted array
!  transformed back into original order and original
!  array (all values should be zero):
   if ((ne2 + ne12 + n0) /= npts) then

      write(0, *) '### Error in mach_hetv_main ###'
      write(0, *) '# number of gridpoints does not match total'
      write(0, *) '# in iso_v.  '
      write(0, *) '# ne2 + ne12: ', ne2 + ne12, ' npts: ', npts
      write(0, *) '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if
!
!  Case 2: 
   if (ne2 >= 1) then
      call mach_hetv_case2(npts, nr, ne2, so4_i, no3_i, nh4_i, hso4_i,       &
                           hno3_i, h_i, nh3_i, lwn_i, t_i, rh_i, ts, ta, tn, &
                           k0, p1, p2, ncas)
      if (chm_error_l) return
!
!  For the case where TA/TS < 2,  check to see if HNO3 dissolution is
!  possible and correct the variables if so:
      call mach_hetv_corrhno3(npts, nr, ne2, 2, so4_i, no3_i, nh4_i, hso4_i, &
                              hno3_i, h_i, lwn_i, t_i, rh_i, rho_i,          &
                              k0, p1, p2, ncas)
      
   end if

!  Case 12: 
   if (ne12 >= 1) then

      call mach_hetv_case12(npts, nr, ne12, so4_i, no3_i, nh4_i, hso4_i, hno3_i,   &
                            h_i, nh3_i, lwn_i, t_i, rh_i, rho_i, ta, ts, tn,       &
                            k0, p1, p2, ncas)

   end if

   if (local_print) then
      if (n0 >= 1) then
         write(chm_lun_out, *)'number of points in case 0 is', n0
      end if
      if (ne2 >= 1) then
         write(chm_lun_out, *)'number of points in case 2 is', ne2 
      end if
      if (ne12 >= 1) then
         write(chm_lun_out, *)'number of points in case 12 is', ne12
      end if
   end if  !diagnostic section
return
end
