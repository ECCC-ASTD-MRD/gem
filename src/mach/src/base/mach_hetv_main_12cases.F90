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
subroutine mach_hetv_main_12cases(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i,  &
                                  h_i, nh3_i, amsul_i, ambis_i, amnit_i,      & 
                                  leto_i, lwn_i, t_i, rh_i, rho_i, case_number)
!!if_off
   use chm_utils_mod,         only: chm_lun_out, chm_error_l
   use mach_hetv_headers_mod, only: mach_hetv_case1, mach_hetv_case2, mach_hetv_case3, mach_hetv_case4,    &
                                    mach_hetv_case5, mach_hetv_case6, mach_hetv_case7, mach_hetv_case8,    &
                                    mach_hetv_case9, mach_hetv_case10, mach_hetv_case11, mach_hetv_case12, &
                                    mach_hetv_corrhno3
   use mach_hetv_mod,         only: tstd, rg
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
   real(kind=8),    intent  (out) :: amsul_i    (npts)
   real(kind=8),    intent  (out) :: ambis_i    (npts)
   real(kind=8),    intent  (out) :: amnit_i    (npts)
   real(kind=8),    intent  (out) :: leto_i     (npts)
   real(kind=8),    intent  (out) :: lwn_i      (npts)
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
   integer(kind=4)               :: i,   j, netc, ic1, ic
   integer(kind=4)               :: ns ,  neco
   integer(kind=4)               :: n0, ne1,  ne2, ne3, ne4, ne5, ne6, ne7, ne8, ne9, ne10
   integer(kind=4)               :: ne11, ne12
   integer(kind=4), dimension(npts) :: ncas, ncor
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
   real(kind=8)                       :: tx1, tx2, tcf

   integer(kind=4),  dimension(npts)  :: locat
   real(kind=8),     dimension(npts)  :: tats, drh_ambis, ta, ts, tn, drh_amnit, &
                                         drh_leto, drh_amsul, mdrh_leto_ambis,   &
                                         mdrh_leto_amsul, mdrh_amnit_amsul
   real(kind=8),     dimension(npts)  :: hno3_c, nh3_c, amsul_c, amnit_c, kamnit, &
                                         ta_c, ts_c, tn_c, zeros                                   
   real(kind=8)                       :: rt_inv
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
!  Set temperatures of crossovers in DRH.  TX1 is the
!  point where the ammonium sulphate and letovicite
!  curves intersect.  TX2 is the point where the
!  ammonium sulphate and ammonium nitrate curves intersect.
!  These values are based on equating the appropriate
!  pair of DRH functions above,  and solving for the temperature
!  where the DRH functions intersect.  They should therefore
!  be updated if the DRH values are updated.

   tx1 = 2.107059277D+02   ! (NH4)2SO4 and (NH4)3H(SO4)2
   tx2 = 2.712042017D+02   ! NH4NO3 and (NH4)2SO4

   netc = 0    !net counter
   n0  = 0
   ns  = 1
   ne1 = 0
   ne2 = 0
   ne3 = 0
   ne4 = 0
   ne5 = 0
   ne6 = 0
   ne7 = 0
   ne8 = 0
   ne9 = 0
   ne10= 0
   ne11= 0
   ne12= 0
   ncor= 0

   do i = 1, npts
!  Initialize output amsul, amnit, ambis, leto, and water to 0 (zero)
      amsul_i = 0.0d0
      ambis_i = 0.0d0
      amnit_i = 0.0d0
      leto_i  = 0.0d0
      lwn_i   = 0.0d0
!  (a) Determine (total ammonia)/(total sulphate)
      ta(i) = nh4_i(i) + nh3_i(i) + 2.0d0 * amsul_i(i) + ambis_i(i) + amnit_i(i) + 3.0d0 * leto_i(i)
      ts(i) = so4_i(i) + hso4_i(i) + amsul_i(i) + ambis_i(i) + 2.0d0 * leto_i(i)
      tn(i) = no3_i(i) + amnit_i(i) + hno3_i(i)
      tcf = (1.0d0 / t_i(i) - 1.0d0 / 298.15d0)
      if (ts(i) > 0.0d0 ) then
          tats(i) = ta(i) / ts(i)
      else
          tats(i) = 0.0d0
      end if
!  (b) Determine the relevant DRH,  MDRH values for each case.  Formulae from ISORROPIA subroutine INIT2
      drh_ambis(i) = 0.4000d0 * exp(384.0d0 * tcf)
      drh_amnit(i) = 0.6183d0 * exp(852.0d0 * tcf)
      drh_leto(i)  = 0.6900d0 * exp(186.0d0 * tcf)
      drh_amsul(i) = 0.7997d0 * exp(80.0d0  * tcf)
      mdrh_leto_ambis(i)  = max(0.01d0, min(drh_ambis(i), drh_leto(i)) - 0.022d0)    !0.3780d0
      mdrh_leto_amsul(i)  = max(0.01d0, min(drh_amsul(i), drh_leto(i)) - 0.02d0)     !0.6700d0  !0.6900d0
      mdrh_amnit_amsul(i) = max(0.01d0, min(drh_amnit(i), drh_amsul(i)) - 0.0183d0)  !0.6000d0
   end do
!
!  #################  TA/TS < 1:  ######################
!
!    case1
!          tats <  1.0,  rh <  drh_ambis          -> calcc1_v
!    case2(a)
!          tats <  1.0,  rh >= drh_ambis          -> calcc2_v
!
!  #################  1 <= TA/TS <1.5 ##################
!
!    case3
!          1.0 <= tats <  1.5,  rh <  mdrh_leto_ambis   -> calcb1ab_v
!    case4
!          1.0 <= tats <  1.5,  c   mdrh_leto_ambis <= rh < drh_ambis,
!              (wet:)  rh < drh_leto           -> calcb1b:
!                                    -> {calcb1ab_v(dry),  calcb2b_v(wet)}
!    case5
!          1.0 <= tats <  1.5,
!          drh_ambis <= rh < drh_leto           -> calcb2b_v
!    case2(b)
!          1.0 <= tats <  1.5,
!          drh_leto <=  rh < drh_amsul          -> calcb3b_v
!    case2(c)
!          1.0 <= tats <  2.0 ,  rh >= drh_amsul         -> calcb4_v
!
!  #################  1.5 <= TA/TS < 2.0  ##############
!
!    case6(a)
!          1.5 <= tats <  2.0,  rh <  mdrh_leto_ambis  -> calcb1aa_v
!    case6(b)
!          1.5 <= tats <  2.0,
!          mdrh_leto_ambis <= rh < drh_ambis,         -> calcb1b:
!                                    -> {calcb1aa_v(dry),  calcb2a("wet")}
!          (note: calcb2a = calcb1aa in terms of system solved);
!           cases 6a,  6b,  6c can be solved using same algorithm as case6a.
!    case6(c)
!          1.5 <= tats <  2.0,
!          drh_ambis <= rh <  mdrh_leto_amsul   -> calcb2a_v
!    case7
!          1.5 <= tats < 2.0,
!          mdrh_leto_amsul <= rh < drh_leto           -> calcb2a2:
!                                    -> {calcb1aa_w(dry),  calcb3a_v(wet)}
!    case8
!          1.5 <= tats <  2.0,
!          drh_leto <= rh <  drh_amsul            -> calcb3a_v
!
!
!
!  #################  TA/TS >= 2.0   ###################
!
!    case9
!          tats >= 2.0,  rh <  mdrh_amnit_amsul   -> calcd1a_v
!    case10
!          tats >= 2.0,
!          mdrh_amnit_amsul <= rh < drh_amnit           -> calcd1a_v(dry),
!                                              (calcd1a_v, calcd2_v(wet))
!    case11
!          tats >= 2.0,
!          drh_amnit <= rh < drh_amsul                   -> calcd1a_v(dry), calcd2_v(wet)
!    case12
!          tats >= 2.0,  rh >= drh_amsul           -> calcd1a_v(dry),  calcd3_v(wet)
!
!  ########################################################
!   GATHER POINTS INTO CONTIGUOUS ARRAYS:  SORTING SECTION
!  #################  TA/TS < 1:  #########################
!
!       tats <  1.0,  rh <  drh_ambis          -> calcc1_v

   do i = 1, npts

      if (tats(i) == 0.0d0) then
         netc = netc + 1
         n0 = n0 + 1
         ncas(i) = 0
         locat(netc) = i
         case_number(i) = 0.0
!  Case 1:
      else if (tats(i) < 1.0d0 .and. rh_i(i) < drh_ambis(i)) then
         netc = netc + 1
         ne1 = ne1 + 1
         ncas(i)= 1
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 1.0
!           tats <  1.0,  rh >= drh_ambis          -> calcc2_v
!     Case 2(a)
      else if (tats(i) < 1.0d0 .and. (rh_i(i) >= drh_ambis(i))) then
         netc = netc + 1
         ne2 = ne2 + 1
         ncas(i)= 2
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 2.1

!    1.0 <= tats <  1.5,
!          drh_leto <=  rh < drh_amsul          -> calcb3b_v
!  Note:  Domain extended to all RH greater than
!  drh_leto.  Same result as combined calcb3b,
!  calcb4,  in same region.  Also avoids crossover
!  problems at temp<211K:  assumption there
!  is that no (NH4)2SO4 forms even below its DRH if
!  1 < TA/TS < 1.5.
!  Same algorithm as case 2a used; use same counters
      else if (tats(i) >= 1.0d0 .and. tats(i) < 1.5d0 .and. (rh_i(i) >= drh_leto(i))) then
         netc = netc + 1
         ne2 = ne2 + 1
         ncas(i)= 2
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 2.2
!  #################  1 <= TA/TS <1.5 ##################
!    1.0 <= tats <  1.5,  rh <  mdrh_leto_ambis   -> calcb1ab_v
      else if (tats(i) >= 1.0d0 .and. tats(i) < 1.5d0 .and. rh_i(i) < mdrh_leto_ambis(i)) then
         netc = netc + 1
         ne3 = ne3 + 1
         ncas(i)= 3
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 3.0
!    1.0 <= tats <  1.5,  c   mdrh_leto_ambis <= rh < drh_ambis,
!              (wet:)  rh < drh_leto           -> calcb1b:
!                                    -> {calcb1ab_v(dry),  calcb2b_v(wet)}

      else if (tats(i) >= 1.0d0 .and. tats(i) < 1.5d0 .and. rh_i(i) >= mdrh_leto_ambis(i) .and. rh_i(i) < drh_ambis(i)) then
         netc = netc + 1
         ne4 = ne4 + 1
         ncas(i)= 4
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 4.0
!    1.0 <= tats <  1.5,
!          drh_ambis <= rh < drh_leto           -> calcb2b_v
      else if (tats(i) >= 1.0d0 .and. tats(i) < 1.5d0 .and. rh_i(i) >= drh_ambis(i) .and. rh_i(i) < drh_leto(i) ) then
         netc = netc + 1
         ne5 = ne5 + 1
         ncas(i)= 5
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 5.0
!  #################  1.5 <= TA/TS < 2.0  ##############
!    1.5 <= tats <  2.0,  rh <  mdrh_leto_ambis  -> calcb1aa_v
!  Case 6,  all together:
      else if (tats(i) >= 1.5d0 .and. tats(i) < 2.0d0 .and. rh_i(i) < mdrh_leto_amsul(i)) then
         netc = netc + 1
         ne6 = ne6 + 1
         ncas(i)= 6
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 6.0
!    1.5 <= tats < 2.0,
!   mdrh_leto_amsul <= rh < drh_leto           -> calcb2a2:
!                                    -> {calcb1aa_w(dry),  calcb3a_v(wet)}

      else if (tats(i) >= 1.5d0 .and. tats(i) < 2.0d0 .and. rh_i(i) >= mdrh_leto_amsul(i) .and. rh_i(i) < drh_leto(i)) then
         netc = netc + 1
         ne7 = ne7 + 1
         ncas(i)= 7
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 7.0
!    1.5 <= tats <  2.0,
!        drh_leto <= rh <  drh_amsul            -> calcb3a_v
      else if (tats(i) >= 1.5d0 .and. tats(i) < 2.0d0 .and. rh_i(i) >= drh_leto(i) .and. rh_i(i) < drh_amsul(i)) then
         netc = netc + 1
         ne8 = ne8 + 1
         ncas(i)= 8
         locat(netc) = i
         ncor(i)= 1
         case_number(i) = 8.0
!  #################  TA/TS >= 2.0   ###################
!           tats >= 2.0,  rh <  mdrh_amnit_amsul   -> calcd1a_v

      else if (tats(i) >= 2.0d0 .and. rh_i(i) < mdrh_amnit_amsul(i) ) then
         netc = netc + 1
         ne9 = ne9 + 1
         ncas(i)= 9
         locat(netc) = i
         case_number(i) = 9.0
!           tats >= 2.0,
! mdrh_amnit_amsul <= rh < drh_amnit           -> calcd1a_v(dry),
!                                              (calcd1a_v, calcd2_v(wet))
      else if (tats(i) >= 2.0d0 .and. rh_i(i) >= mdrh_amnit_amsul(i) .and. rh_i(i) < drh_amnit(i) ) then
         netc = netc + 1
         ne10 = ne10 + 1
         ncas(i)= 10
         locat(netc) = i
         case_number(i) = 10.0
!           tats >= 2.0,
! drh_amnit <= rh < drh_amsul                   -> calcd1a_v(dry), calcd2_v(wet)
      else if (tats(i) >= 2.0d0 .and. rh_i(i) >= drh_amnit(i) .and. rh_i(i) < drh_amsul(i) ) then
         netc = netc + 1
         ne11 = ne11 + 1
         ncas(i)= 11
         locat(netc) = i
         case_number(i) = 11.0
!  Correction:  case 12 used for the following, rather than case 12,
!  to improve gas-phase ammonia concentration calculation. P.A. Makar, Nov. 2008
      else if (tats(i) >= 1.5d0 .and. tats(i) < 2.0d0 .and. &
          ((rh_i(i) >= drh_amsul(i) .and. t_i(i) > tx1) .or. &
           (rh_i(i) >= drh_leto(i) .and. t_i(i) <= tx1))  ) then
         netc = netc + 1
         ne12 = ne12 + 1
         ncas(i)= 12
         locat(netc) = i
         case_number(i) = 12.1
!           tats >= 2.0,  rh >= drh_amsul       -> calcd1a_v(dry),  calcd3_v(wet)
      else if (tats(i) >= 2.0d0 .and. ((rh_i(i) >= drh_amsul(i) .and. t_i(i) > tx2) .or.  &
         (rh_i(i) >= drh_amnit(i) .and. t_i(i) <= tx2))) then
         netc = netc + 1
         ne12 = ne12 + 1
         ncas(i)= 12
         locat(netc) = i
         case_number(i) =12.2
      end if
   end do

!  ####################################################
!  Check:  calculate difference between sorted array
!  transformed back into original order and original
!  array (all values should be zero):
   if (netc /= npts) then

      write(0, *) '### Error in mach_hetv_main_12cases ###'
      write(0, *) '# number of gridpoints does not match total'
      write(0, *) '# in iso_v.  '
      write(0, *) '# netc: ', netc, ' npts: ', npts

!  Find the gridpoints that missed being classified, and write out their concentrations,  etc.,  before stopping.
      do i = 1, npts
         ic1 = 0
         do j = 1, netc
            if (locat(j) == i) ic1 = ic1 + 1
         end do
         if (ic1 == 0) then
         write (0, *)'# *************************'
         write (0, *)'# point: ', i, ' was not classified.'
         write (0, *)'# input data for this point: '
         write (0, 2271)'# so4', so4_i(i)
         write (0, 2271)'# nh4', nh4_i(i)
         write (0, 2271)'# no3', no3_i(i)
         write (0, 2271)'# hso4', hso4_i(i)
         write (0, 2271)'# amsul', amsul_i(i)
         write (0, 2271)'# amnit', amnit_i(i)
         write (0, 2271)'# ambis', ambis_i(i)
         write (0, 2271)'# leto', leto_i(i)
         write (0, 2271)'# hno3', hno3_i(i)
         write (0, 2271)'# nh3', nh3_i(i)
         write (0, 2271)'# h', h_i(i)
         write (0, 2271)'# rh', rh_i(i)
         write (0, 2271)'# t', t_i(i)
         write (0, 2271)'# rho', rho_i(i)
         write (0, 2271)'# mdrh_leto_ambis', mdrh_leto_ambis(i)
         write (0, 2271)'# drh_ambis', drh_ambis(i)
         write (0, 2271)'# mdrh_leto_amsul', mdrh_leto_amsul(i)
         write (0, 2271)'# drh_leto', drh_leto(i)
         write (0, 2271)'# mdrh_amnit_amsul', mdrh_amnit_amsul(i)
         write (0, 2271)'# drh_amnit', drh_amnit(i)
         write (0, 2271)'# drh_amsul', drh_amsul(i)
         write (0, *)'# *************************'
         end if
 2271           format (2x, a16, 3x, 1pe15.8)
      end do
!  check to make sure that double classifications have not occurred:
      do i = 1, npts
         ic = 0
         do j = 1, netc
            if (locat(j) == i) ic = ic + 1
         end do
         if (ic > 1) then
            write(0, *)    '# double assignment for gridpoint: ', i
            write(0, *)    '# values: '
            write(0, 2271) '# so4', so4_i(i)
            write(0, 2271) '# nh4', nh4_i(i)
            write(0, 2271) '# no3', no3_i(i)
            write(0, 2271) '# hso4', hso4_i(i)
            write(0, 2271) '# amsul', amsul_i(i)
            write(0, 2271) '# amnit', amnit_i(i)
            write(0, 2271) '# ambis', ambis_i(i)
            write(0, 2271) '# leto', leto_i(i)
            write(0, 2271) '# hno3', hno3_i(i)
            write(0, 2271) '# nh3', nh3_i(i)
            write(0, 2271) '# h', h_i(i)
            write(0, 2271) '# rh', rh_i(i)
            write(0, 2271) '# t', t_i(i)
            write(0, 2271) '# rho', rho_i(i)
            write(0, 2271) '# mdrh_leto_ambis', mdrh_leto_ambis(i)
            write(0, 2271) '# mdrh_leto_amsul', mdrh_leto_amsul(i)
            write(0, 2271) '# mdrh_amnit_amsul', mdrh_amnit_amsul(i)
            write(0, 2271) '# drh_ambis', drh_ambis(i)
            write(0, 2271) '# drh_leto', drh_leto(i)
            write(0, 2271) '# drh_amnit', drh_amnit(i)
            write(0, 2271) '# drh_amsul', drh_amsul(i)
            write(0, *)    '###         ABORT         ###'
            chm_error_l = .true.
            return
         end if
      end do
      write(0, *)    '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if

   if (local_print) then

      write(chm_lun_out, 202)'point no.', 'ta/ts', 'rh'
      do i = 1, npts
         write(chm_lun_out, 203)i, tats(i), rh_i(i)
      end do

      if (n0 >= 1) then
         write(chm_lun_out, *)'number of points in case 0 is', n0 
      end if
      if (ne1 >= 1) then
         write(chm_lun_out, *)'number of points in case 1 is', ne1 
      end if
      if (ne2 >= 1) then
         write(chm_lun_out, *)'number of points in case 2 is', ne2
      end if
      if (ne3 >= 1) then
         write(chm_lun_out, *)'number of points in case 3 is', ne3
      end if
      if (ne4 >= 1) then
         write(chm_lun_out, *)'number of points in case 4 is', ne4
      end if
      if (ne5 >= 1) then
         write(chm_lun_out, *)'number of points in case 5 is', ne5
      end if
      if (ne6 >= 1) then
         write(chm_lun_out, *)'number of points in case 6 is', ne6
      end if
      if (ne7 >= 1) then
         write(chm_lun_out, *)'number of points in case 7 is', ne7 
      end if
      if (ne8 >= 1) then
         write(chm_lun_out, *)'number of points in case 8 is', ne8
      end if
      if (ne9 >= 1) then
         write(chm_lun_out, *)'number of points in case 9 is', ne9
      end if
      if (ne10 >= 1) then
         write(chm_lun_out, *)'number of points in case 10 is', ne10
      end if
      if (ne11 >= 1) then
         write(chm_lun_out, *)'number of points in case 11 is', ne11
      end if
      if (ne12 >= 1) then
         write(chm_lun_out, *)'number of points in case 12 is', ne12
      end if
   end if  !diagnostic section

!  for each of the subsystems identified,  call the appropriate
!  (vectorized over gridpoint) subroutine:

!  Case 1:  tats <  1.0,  rh <  drh_ambis
   if (ne1 > 0) then

      call mach_hetv_case1(npts, nr, ne1, so4_i, no3_i, nh4_i, hso4_i,    &
                           hno3_i, h_i, nh3_i, ambis_i, lwn_i, t_i, rh_i, &
                           ts, ta, tn, k0, p1, p2, ncas)

   end if

!  Case 2:
!  2(a):  tats <  1.0,  rh >= drh_ambis
!  2(b):  1.0 <= tats <  1.5, ,  drh_leto <=  rh < drh_amsul
!  2(c):  1.0 <= tats <  2.0 ,  rh >= drh_amsul
   if (ne2 > 0) then

      call mach_hetv_case2(npts, nr, ne2, so4_i, no3_i, nh4_i, hso4_i,       &
                           hno3_i, h_i, nh3_i, lwn_i, t_i, rh_i, ts, ta, tn, &
                           k0, p1, p2, ncas)
      if (chm_error_l) return
   end if

!  Case 3: 1.0 <= tats <  1.5,  rh <  mdrh_leto_ambis
   if (ne3 > 0) then

      call mach_hetv_case3(npts, so4_i, no3_i, nh4_i, hso4_i, hno3_i, h_i, &
                           nh3_i, ambis_i, leto_i, lwn_i, ts, ta, tn, ncas)
   end if

!  Case 4: 1.0 <= tats <  1.5,  c   mdrh_leto_ambis <= rh < drh_ambis
!
   if (ne4 > 0) then

      call mach_hetv_case4(npts, nr, ne4, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, ambis_i, leto_i, mdrh_leto_ambis,       &
                           drh_ambis, lwn_i, t_i, rh_i, ts, ta, tn,            &
                           k0, p1, p2, ncas)
      if (chm_error_l) return
   end if

!  Case 5: 1.0 <= tats <  1.5,  drh_ambis <= rh < drh_leto
   if (ne5 > 0) then

      call mach_hetv_case5(npts, nr, ne5, 5, so4_i, no3_i, nh4_i, hso4_i, &
                           hno3_i, h_i, nh3_i, leto_i, lwn_i, t_i, rh_i,  &
                           ts, ta, tn, k0, p1, p2, ncas)
   end if

!  Case 6:
!  6(a):  1.5 <= tats <  2.0,  rh <  mdrh_leto_ambis
!  6(b):  1.5 <= tats <  2.0,  mdrh_leto_ambis <= rh < drh_ambis
!  6(c):  1.5 <= tats <  2.0,  drh_ambis <= rh <  mdrh_leto_amsul
   if (ne6 > 0) then

      call mach_hetv_case6(npts, so4_i, no3_i, nh4_i, hso4_i,   &
                           hno3_i, h_i, nh3_i, amsul_i, leto_i, &
                           lwn_i, ts, ta, tn, ncas)
   end if

!  Case 7: 1.5 <= tats < 2.0,  mdrh_leto_amsul <= rh < drh_leto
   if (ne7 > 0) then

      call mach_hetv_case7(npts, nr, ne7, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, amsul_i, leto_i, mdrh_leto_amsul,       &
                           drh_leto, lwn_i, t_i, rh_i, ts, ta, tn,             &
                           k0, p1, p2, ncas)
      if (chm_error_l) return
   end if

!  Case 8: 1.5 <= tats <  2.0,  drh_leto <= rh <  drh_amsul
   if (ne8 > 0) then

      call mach_hetv_case8(npts, nr, ne8, 8, so4_i, no3_i, nh4_i, hso4_i, &
                           hno3_i, h_i, nh3_i, amsul_i, lwn_i, t_i, rh_i, &
                           ts, ta, tn, k0, p1, p2, ncas)
   end if

   if (ne9 > 0) then
      zeros = 0.0d0
!  Calculate rates for Equilibrium constants
!  9:  NH4NO3(s) <=> NH3(g_eq) + HNO3(g_eq)
      where (ncas == 9)
        kamnit = k0(9) * exp(p1(9) * (tstd / t_i - 1.0d0) +   &
                 p2(9) * (1.0d0 + log(tstd / t_i) - tstd / t_i))
!  Original units are atm^2:  Gas is in molal units already; convert rate constant to (kg^3 air)/mole:
        kamnit = kamnit / ((rg * t_i * rho_i) * (rg * t_i * rho_i))
      elsewhere
        kamnit = 0.0d0
      end where
      hno3_c  = pack(hno3_i,  ncas == 9, zeros)
      nh3_c   = pack(nh3_i,   ncas == 9, zeros)
      amsul_c = pack(amsul_i, ncas == 9, zeros)
      amnit_c = pack(amnit_i, ncas == 9, zeros)
      ta_c    = pack(ta,      ncas == 9, zeros)
      ts_c    = pack(ts,      ncas == 9, zeros)
      tn_c    = pack(tn,      ncas == 9, zeros)
      kamnit  = pack(kamnit,  ncas == 9, zeros)
      call mach_hetv_case9(npts, ne9, hno3_c, nh3_c, amsul_c, amnit_c, &
                           kamnit, ta_c, ts_c, tn_c)
!
      hno3_i  = unpack(hno3_c,  ncas == 9, hno3_i)
      nh3_i   = unpack(nh3_c,   ncas == 9, nh3_i)
      amsul_i = unpack(amsul_c, ncas == 9, amsul_i)
      amnit_i = unpack(amnit_c, ncas == 9, amnit_i)
   end if

!  Case 11: tats >= 2.0,  drh_amnit <= rh < drh_amsul
   if (ne11 > 0) then

      call mach_hetv_case11(npts, nr, ne11, 11, so4_i, no3_i, nh4_i, hso4_i, &
                            hno3_i, h_i,  nh3_i, amsul_i, lwn_i, t_i, rh_i,  &
                            rho_i, ts, ta, tn, k0, p1, p2, ncas)
      if (chm_error_l) return

   end if

!  Case 10:  tats >= 2.0,  mdrh_amnit_amsul <= rh < drh_amnit
   if (ne10 > 0) then

      call mach_hetv_case10(npts, nr, ne10, so4_i, no3_i, nh4_i, hso4_i, &
                            hno3_i, h_i, nh3_i, amsul_i, amnit_i,        &
                            mdrh_amnit_amsul, drh_amnit, lwn_i, t_i,     &
                            rh_i, rho_i, ts, ta, tn, k0, p1, p2, ncas)
      if (chm_error_l) return

   end if
!  Case 12: tats >= 2.0,  rh >= drh_amsul
   if (ne12 > 0) then

      call mach_hetv_case12(npts, nr, ne12, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                            h_i, nh3_i, lwn_i, t_i, rh_i, rho_i, ta, ts, tn,     &
                            k0, p1, p2, ncas)

   end if

!  For all cases where TA/TS < 2,  check to see if HNO3 dissolution is
!  possible and correct the variables if so:
   neco = ne1 + ne2 + ne3 + ne4 + ne5 + ne6 + ne7 + ne8
   if (neco > 0) then
      call mach_hetv_corrhno3(npts, nr, neco, 1, so4_i, no3_i, nh4_i, hso4_i, hno3_i,  &
                              h_i, lwn_i, t_i, rh_i, rho_i, k0, p1, p2, ncor)
   end if

202   format (3(2x, a9, 2x))
203   format (2x, i9, 2x, 2(2x, 1pe9.2, 2x))
302   format (3(2x, a9, 2x))
304   format (2x, i9, 2x, 2(2x, 1pe24.17, 2x))
992   format (i4, 2x, 10(1x, 1pe8.1, 1x))

  return
end subroutine mach_hetv_main_12cases
