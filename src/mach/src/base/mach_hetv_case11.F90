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
! Fichier/File   : mach_hetv_case11.ftn90
! Creation       : P. Makar, V. Bouchet, A. Nenes, S. Gravel, B. Pabla, S .Menard
! Description    : Heterogeneous chemistry solver for case11.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne11 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case 11's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne11.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                  TA/TS >= 2.0,  drh_amnit <= rh < drh_amsul
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!                  (1)
!                   (NH4)2SO4 <=> 2 NH4 + SO4,  k(NH4)2SO4
!                  (2)
!                    HNO3 <=> H + NO3, kHNO3
!                    NH3 + H2O <=> NH4 + OH, kNH3
!                    H2O <=> H + OH, kH2O
!                  (3)
!                    H + NH4 = NO3 + 2 SO4 + OH = 0
!                  (4)
!                     NH3 + H2O <=> NH4 + kH2O H2O/H, kNH3
!                  (5)
!                     HSO4 <=> H + SO4, kHSO4
!
!                  The system of equations to be solved in this case are as follows:
!
!                  (1)
!                    NH4available = NH4NO3dry + (NH3
!                    (delNH42SO4)^3 + NH4available*(delNH42SO4)^2
!                                 + 1/4 (NH4available)^2 * delNH42SO4
!                                                  -    k(NH4)2SO4/4 = 0
!                    SO4 = SO4dry + delNH42SO4,
!                    NH4 = NH4available + 2 delNH42SO4
!                    NH3 = NH3dry - delNH3
!                    (NH4)2SO4 = (NH4)2SO4dry - delNH42SO4
!                  (2)
!                    HNO3 = HNO3dry - delNO3,
!                    NO3 = NH4NO3dry +  delNO3,
!                    kHNO3 HNO3 = (NH4)(NO3)kH2O/(NH3 kNH3)
!                  (3)
!                    H2 + (NH4 - NO3 - 2 SO4) H - kH2O = 0
!                  (4)
!                    NH3 = (NH4)(OH)/(kNH3 H2O),
!                    OH/H2O = kH2O/H,
!                  (5)
!                    kHSO4 (HSO4) =
!                         (H-delHSO4)(SO4 - delHSO4)
!                    Hnew = H -delHSO4,
!                    SO4new = SO4 -delHSO4
!                    HSO4 = HSO4in + delHSO4
!
!                   The solution to the system of equations:
!
!                  (1)
!                    delNH42SO4 = largest positive real root of cubic equation or zero.
!                    Note: SO4dry  = 0
!                    SO4 = delNH42SO4,
!                    NH4 = NH4NO3dry + delNH3 + 2 delNH42SO4
!                    NH3 = NH3dry - delNH3
!                    (NH4)2SO4 = (NH4)2SO4dry - delNH42SO4
!                  (2)
!                    delNO3 = [kHNO3 kNH3 (NH3)(HNO3dry)
!                                    -  kH2O(NH4) (NH4NO3dry)]
!                       /   [kH2O (NH4) + kHNO3 kNH3 NH3]
!                    HNO3 = HNO3dry - delNO3,
!                    NO3 = NH4NO3dry +  delNO3.
!                  (3)
!                    H = 0.5 {-(NH4 - NO3 - 2 SO4) + sqrt( (NH4 - NO3 - 2 SO4)^2 + 4 kH2O)
!                  (4)
!                    (NH3)(H)kNH3 - (NH4)kH2O = 0
!                  (5)
!                    delHSO4 = smallest positive real root of
!                    delHSO4^2 - (H + SO4 + kHSO4) delHSO4
!                      + [(H)(SO4)] = 0, or zero.
!                    Hnew = H -delHSO4,
!                    SO4new = SO4 -delHSO4
!                    Variable, Range
!                    delNH3, 0 to NH3dry
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case11(npts, nr, ne11, nc, so4_i, no3_i, nh4_i, hso4_i, &
                            hno3_i, h_i, nh3_i, amsul_i, lwn_i, t_i, rh_i,   &
                            rho_i, ts_i, ta_i, tn_i, k0, p1, p2, ncas)
!!if_off
   use chm_utils_mod,         only: chm_error_l
   use mach_hetv_headers_mod, only: mach_hetv_case9, mach_hetv_water, &
                                    mach_hetv_activity, mach_hetv_poly3v, &
                                    mach_hetv_corrhno3
   use mach_hetv_mod,         only: lwmin, tstd, small, smrt, iter, itero, ndiv, &
                                    eps, eps2, rg
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne11
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: ncas   (npts)
   real(kind=8),    intent   (in) :: k0     (nr)
   real(kind=8),    intent   (in) :: p1     (nr)
   real(kind=8),    intent   (in) :: p2     (nr)
   real(kind=8),    intent(inout) :: so4_i  (npts)
   real(kind=8),    intent(inout) :: no3_i  (npts)
   real(kind=8),    intent(inout) :: nh4_i  (npts)
   real(kind=8),    intent(inout) :: hso4_i (npts)
   real(kind=8),    intent(inout) :: hno3_i (npts)
   real(kind=8),    intent(inout) :: h_i    (npts)
   real(kind=8),    intent(inout) :: nh3_i  (npts)
   real(kind=8),    intent(inout) :: amsul_i(npts)
   real(kind=8),    intent(inout) :: lwn_i  (npts)
   real(kind=8),    intent   (in) :: t_i    (npts)
   real(kind=8),    intent   (in) :: rh_i   (npts)
   real(kind=8),    intent   (in) :: rho_i  (npts)
   real(kind=8),    intent   (in) :: ta_i   (npts)
   real(kind=8),    intent   (in) :: ts_i   (npts)
   real(kind=8),    intent   (in) :: tn_i   (npts)
!!if_off
!
!  Local variables:
!
   integer(kind=4) :: i, jo, rooteval, rootsrch, jsum, iwe
   real(kind=8) :: dndiv
   real(kind=8), parameter :: h2o = 1000.D0 / 18.01528D0    !moles H2O/kg H2o
   real(kind=8) :: delouter, del, a, b, c, d, v, root
   real(kind=8), dimension(ne11) :: so4, no3, nh4, hso4, amsul, hno3, nh3, &
                                    h, lwn, t, aw, rho, ta, ts, tn,        &
                                    nh3dry, hno3dry, amnitdry, amsuldry
   real(kind=8), dimension(npts) :: awu, law
   real(kind=8), dimension(npts) :: a2, b2, c2, delamsul2
   real(kind=8), dimension(ne11) :: g_h_no3, g_h_no3_o, g_h_hso4, g_h_hso4_o, &
                                    g_nh42_so4, g_nh42_so4_o,                 &
                                    g_nh4_no3, g_nh4_no3_o,                   &
                                    g_h2_so4, g_h2_so4_o,                     &
                                    lwo, delnh3_lo, delnh3_hi
   real(kind=8), dimension(ne11) :: choice, rtbis, dx, diff_lo, diff, delno3, delamsul,       &
                                    delnh3, nh4available, difflo, diff_1, locmin, avdiff
   real(kind=8), dimension(ne11) :: khso4, kh2o, knh3, khno3, kamsul, kamnit
   real(kind=8) :: khso4_m, knh3_m, khno3_m, kamsul_m, work1, work2, work3, lwo_lwn
   integer, dimension(ne11)      :: ind2, rflag, rootloc, bflag, wflag, rootmin, islv, islv2
! For HNO3 correction:
   integer(kind=4) :: nh4_counter
   integer(kind=4), dimension(ne11) :: ncas11_cor
!
!   so4   = pack(so4_i,   ncas == nc)
!   nh4   = pack(nh4_i,   ncas == nc)
!   no3   = pack(no3_i,   ncas == nc)
!   hso4  = pack(hso4_i,  ncas == nc)
   hno3  = pack(hno3_i,  ncas == nc)
!   h     = pack(h_i,     ncas == nc)
   nh3   = pack(nh3_i,   ncas == nc)
!   amsul = pack(amsul_i, ncas == nc)
   ta    = pack(ta_i,    ncas == nc)
   ts    = pack(ts_i,    ncas == nc)
   tn    = pack(tn_i,    ncas == nc)
   t     = pack(t_i,     ncas == nc)
   aw    = pack(rh_i,    ncas == nc)
   rho   = pack(rho_i,   ncas == nc)

!  NOTE:  THIS OPTION NOT ACTIVE IN THIS VERSION
!  If the results of the dry calculation are such that there is no NH4NO3(s),
!  then deliquescence will not occur; the correct solution to the system of
!  equations is the dry calculation. Although this will be calculated correctly
!  in the code that follows, its more efficient to reorganize the input arrays
!  at this stage, and only perform the calculations for those cases in which
!  the dry NH4NO3 is non-zero.

   do i = 1, ne11
!  Calculate rates for Equilibrium constants
      work1 = tstd / t(i) - 1.0d0
      work2 = 1.0d0 + log(tstd / t(i)) - tstd / t(i)
      work3 = rg * t(i) * rho(i)
!  1:  HNO3(g_eq) <=> HNO3(aq)
!  2:  HNO3(aq) <=> H+ + NO3-
!  Net  HNO3(g_eq) <=> H+ + NO3- : khno3=k1*k2
      khno3(i) = (k0(1) * exp(p1(1) * work1 + p2(1) * work2)) * &
                 (k0(2) * exp(p1(2) * work1 + p2(2) * work2)) * work3
!  3:  NH3(g_eq) <=> NH3(aq)
!  4:  NH3(aq) + H2O <=> NH4+ + OH-
!  Net NH3(g_eq) + H2O <=> NH4+ + OH-:  knh3=k3*k4
      knh3(i) = (k0(3) * exp(p1(3) * work1 + p2(3) * work2)) * &
                (k0(4) * exp(p1(4) * work1 + p2(4) * work2)) * work3
!  5:  HSO4- <=> SO4= + H+
      khso4(i) = k0(5) * exp(p1(5) * work1 + p2(5) * work2)
!  7:  (NH4)2SO4(s) <=> 2 NH4+ + SO4=
      kamsul(i) = k0(7) * exp(p1(7) * work1 + p2(7) * work2)
!  9:  NH4NO3(s) <=> NH3(g_eq) + HNO3(g_eq)
      kamnit(i) = k0(9) * exp(p1(9) * work1 + p2(9) * work2)
!  Original units are atm^2:  Gas is in molal units already; convert rate constant to (kg^3 air)/mole:
      kamnit(i) = kamnit(i) / (work3**2)
!
!  Calculate OH concentration for activity coefficient calculations
      kh2o(i) = k0(10) * exp(p1(10) * work1 + p2(10) * work2)
!
!  polynomials for electrolytes only good between aw=0.002 and aw=0.98
      awu(i) = min(max(aw(i), 2.0d-03), 0.98d0)
      law(i) = dble(log(real(awu(i))))    !double precision log unnec.
!
!  Initialize dry particle gases and components
      hno3dry(i)  = 0.0d0
      nh3dry(i)   = 0.0d0
      amsuldry(i) = 0.0d0
      amnitdry(i) = 0.0d0
   end do
!  Calculate dry particle gases and components using Case 9.  Note that the "dry" values of nh3, hno3,
!  and nh4no3 are the only parts of this solution that are used here.
   call mach_hetv_case9(ne11, ne11, hno3dry, nh3dry, amsuldry, amnitdry, &
                        kamnit, ta, ts, tn)

!  Set value of subdivision for root search:
   dndiv = 1.0D0 / dble(ndiv)
!  Initial guess for components:  SO4 = 0, HSO4 = 0,
!  NH4 = NH4NO3dry, NO3 = NH4NO3dry, H = 0
   do i = 1, ne11
      so4(i) = 0.D0
      hso4(i) = 0.D0
      nh4(i) = amnitdry(i)
      h(i) = 0.D0
      no3(i) = amnitdry(i)
      amsul(i) = 0.D0
!
      lwn(i) = 1.D0
      wflag(i) = 0
      delnh3_hi(i) = 1.D0
      delnh3_lo(i) = 1.D0

!  Set value of subdivision for root search:
      dx(i) = (nh3dry(i) + 2.D0 * amsuldry(i) + amnitdry(i)) * dndiv
    end do

!  Looping point for outer iteration (subdivision
!  of function interval and search for roots):
   rooteval = 0
   rootsrch = 0
   delouter = 1.0d0

!  Repeat the calculations, for all ndiv+1 initial function evaluations:
   do while ((rooteval <= ndiv + 1) .or. (rootsrch <= 3) .or. &
             (delouter > eps .and. rootsrch < iter))

      do i = 1, ne11
!   Initial guess for liquid water content (will be replaced
!   in calculations to follow:
         lwo(i) = 1.D0
!
!   Reset concentrations to starting liquid water content
         nh3(i)   = nh3(i)  * lwn(i)
         hno3(i)  = hno3(i) * lwn(i)
         so4(i)   = 0.D0  !so4(i)   * lwn(i)
         h(i)     = 0.D0  !h(i)     * lwn(i)
         no3(i)   = amnitdry(i) * lwn(i) !no3(i)   * lwn(i)
         nh4(i)   = amnitdry(i) * lwn(i) !nh4(i)   * lwn(i)
         hso4(i)  = 0.D0  !hso4(i)  * lwn(i)
         amsul(i) = 0.D0  !amsul(i) * lwn(i)
         nh3dry(i)  = nh3dry(i)  * lwn(i)
         hno3dry(i) = hno3dry(i) * lwn(i)
         amnitdry(i) = amnitdry(i) * lwn(i)
         amsuldry(i) = amsuldry(i) * lwn(i)
         delnh3(i)   = delnh3(i) * lwn(i)

!  Zero the 'old value' activity coefficient arrays
         g_h_hso4_o(i)   = 0.D0
         g_h2_so4_o(i)   = 0.D0
         g_h_no3_o(i)    = 0.D0
         g_nh4_no3_o(i)  = 0.D0
         g_nh42_so4_o(i) = 0.D0
      end do
!
      rooteval = rooteval + 1

!  First root evaluation is for lower bound,
!   delnh3 = - (2.D0*amsuldry(i)+amnitdry(i))
      if (rooteval == 1) then
         do i = 1, ne11
            delnh3(i) = - (2.D0 * amsuldry(i) + amnitdry(i))
         end do
      end if

      if (rooteval == 2) then
!  Store results of first root evaluation:
         do i = 1, ne11
            difflo(i) = diff_lo(i)
            diff_1(i) = diff_lo(i)
!
!  Set the indirect address for the location
!  of the sign change (rootloc) and the counter for the number of
!  sign changes (rflag) for the given gridpoint to zero:
            rootloc(i) = 0
            rflag(i)   = 0
            rootmin(i) = 0
            locmin(i)  = 1.D30
         end do
      end if

!  Subsequent root evaluations (but before the binary search
!  refinement):  determine the value of the variable to be
!  used for the next function evaluation:
      if (rooteval > 1 .and. rooteval <=  ndiv + 1) then
         do i = 1, ne11
            delnh3(i) = - (2.D0 * amsuldry(i) + amnitdry(i)) + dx(i) * dble(rooteval - 1)
         end do
      end if
!
      if (rooteval > 2 .and. rooteval <=  ndiv + 2) then
!  Determine whether the boundaries of the current sub interval
!  constitute the minimum for the entire interval.  If so,
!  store the right-hand index:
         do i = 1, ne11
            avdiff(i) = 0.5D0 * abs(diff(i) + diff_1(i))
            if (avdiff(i) <=  locmin(i)) then
               rootmin(i) = rooteval - 1
               locmin(i) = avdiff(i)
            end if
         end do

!  Determine whether a sign change occurred in the last
!  interval:
         do i = 1, ne11
            if (sign(1.d0, diff_1(i)) * sign(1.d0, diff(i)) <  0.d0) then
               rflag(i) = rflag(i) + 1
            end if
         end do

!  If it was the first sign change detected for the given gridpoint, save its location:

         do i = 1, ne11
            if (rflag(i) == 1) then
               rootloc(i) = rooteval - 1
               rflag(i) = 2
            end if
         end do

!  Update the variable and evaluated function arrays:
         do i = 1, ne11
            diff_1(i) = diff(i)
         end do
      end if

!  The next section determines the boundaries to be used in the
!  binary search refinement of the given root.  Sometimes both limits
!  of the search are set to the same value (upper or lower boundary
!  of the entire interval), if no sign change occurred in the entire interval.
      if (rooteval == ndiv + 2) then

!  Set the "boundary has been chosen" flag to "off":
         do i = 1, ne11
            bflag(i) = 0
         end do

! following block changed (see paul's e-mail dated Jan 15/03)
!  Use the lower limit of the variable space, if that lower limit
!  corresponds to a root:
         do i = 1, ne11
            if (abs(difflo(i)) <=  eps .and. rootloc(i) == 0) then
               delnh3_lo(i) = - (2.D0 * amsuldry(i) + amnitdry(i))
               delnh3_hi(i) = - (2.D0 * amsuldry(i) + amnitdry(i))
               bflag(i)  = 1
            end if
         end do

!  If a root was found in the variable space for the function at the
!  given gridpoint, use the stored index to give the boundaries:

         do i = 1, ne11
            if (rootloc(i) /= 0 .and. bflag(i) == 0) then
               delnh3_lo(i) = - (2.D0 * amsuldry(i) + amnitdry(i)) + dx(i) * dble(rootloc(i) - 2)
               delnh3_hi(i) = delnh3_lo(i) + dx(i)
               bflag(i) = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  but the upper limit is a root, use that value for the root:
         do i = 1, ne11
            if (bflag(i) == 0 .and. abs(diff(i)) <=  eps .and. rootloc(i) == 0) then
               delnh3_lo(i) = delnh3(i)
               delnh3_hi(i) = delnh3(i)
               bflag(i)  = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  and all the values of the root evaluations were positive, then
!  K[NH3] (NH3)(H) - (NH4) K[H2O] > 0 always - this implies that
!  there will always be higher NH3 concentrations than is necessary
!  to satisfy the NH3 and water equilibrium equations - while
!  local minima (not roots) may exist, the chemically correct course
!  in this case is to assume that NH3 is at the minimum
!  value.  This in turn requires that the variable, del_nh3 (the amount
!  of NH3 lost from the dry phase when NH4NO3 deliquesces) is set to the
!  maximum value; NH3 = NH3dry - del_nh3.  In order to ensure that this
!  is done, both boundaries of the root search are set to the maximum value
!  of del_nh3:

         do i = 1, ne11
            if (bflag(i) == 0 .and. diff(i) > 0.0d0) then
               delnh3_lo(i) = delnh3(i)
               delnh3_hi(i) = delnh3(i)
               bflag(i) = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  and all the values of the root evaluations were negative, then
!  K[NH3] (NH3)(H) - (NH4) K[H2O] < 0 always - this implies that
!  there will always be lower NH3 concentrations than is necessary
!  to satisfy the NH3 and water equilibrium equations - while
!  local minima (not roots) may exist, the chemically correct course
!  in this case is to assume that NH3 is at the maximum
!  value.  This in turn requires that the variable, del_nh3 (the amount
!  of NH3 lost from the dry phase when NH4NO3 deliquesces) is set to the
!  minimum value; NH3 = NH3dry - del_nh3.  In order to ensure that this
!  is done, both boundaries of the root search are set to the minimum value
!  of del_nh3:

         do i = 1, ne11
            if (bflag(i) == 0 .and. diff(i) <  0.0d0) then
               delnh3_lo(i) = - (2.D0 * amsuldry(i) + amnitdry(i))
               delnh3_hi(i) = - (2.D0 * amsuldry(i) + amnitdry(i))
               bflag(i) = 1
            end if
         end do

!  Check to make sure that some boundaries for the binary search were found:
         jsum = sum(bflag)
         if (jsum < ne11) then
            write(*, *) '> WARNING <'
            write(*, *) '> Root interval not set in Case 1'
            write(*, *) '> Condition occurred at the following gridpoint:'
            do i = 1, ne11
               if (bflag(i) == 0) write(*, *) '> ', i
            end do
         end if
      end if

!  If the initial subdivision of the variable space has been completed,
!  start the bisection search refinement of the chosen function interval

      if (rooteval >= ndiv + 2) rootsrch = rootsrch + 1

!  Reevaluate the function on the lhs of the interval:
      if (rootsrch == 1) then
         do i = 1, ne11
            delnh3(i) = delnh3_lo(i)
         end do
      end if

!  Reevaluate the function on the rhs of the interval:
      if (rootsrch == 2) then
         do i = 1, ne11
            delnh3(i) = delnh3_hi(i)
         end do
      end if

!  On the third function evaluation, the bisection search
!  has commenced.  The first stage of this is to orient the
!  search so that the positive value of the root equation
!  lies on the right of the interval.  Note that from this
!  point forward, the "_hi" value ceases to be the upper
!  limit of the variable - it is now the increment in the
!  variable as the function proceeds towards the root.
      if (rootsrch == 3) then
         do i = 1, ne11
!  "choice" is -1 if diff_lo<0, 1 otherwise:
            choice(i) = sign (1.E+0, real(diff_lo(i)))
            dx(i) = choice(i) * (delnh3_lo(i) - delnh3_hi(i))

!  "choice" is 1 if diff_lo<0, 0 otherwise
            choice(i) = 0.5d0 * (1.0d0 - choice(i))
            rtbis(i)  = choice(i) * delnh3_lo(i) + &
                        (1.0d0 - choice(i)) * delnh3_hi(i)
            delnh3_hi(i) = 0.5d0 * dx(i)
            delnh3_lo(i) = rtbis(i)
            delnh3(i) = delnh3_lo(i) + delnh3_hi(i)
!  rtbis is used to store the value of dx for future convergence calcs.
            rtbis(i) = max(abs(dx(i)), small)
         end do
      end if

!  All subsequent root evaluations:  determine which side of the
!  interval contains the root crossing, rearrange the bounds, and
!  subdivide the interval again:

      delouter = 0.D0
      if (rootsrch > 3) then
         do i = 1, ne11
!  "choice" is 1 if diff<0, 0 otherwise:
            choice(i) = 0.5d0 * (1.0d0 - sign (1.E+0, real(diff(i))))

!  don't update if old value of delhi/dello < eps
            if (abs(delnh3_hi(i) / rtbis(i))  >  eps) then
               delnh3_lo(i) = delnh3(i) * choice(i) + &
                              (1.0d0 - choice(i)) * delnh3_lo(i)
               delnh3_hi(i) = delnh3_hi(i) * 0.5D0
               delnh3(i)    = delnh3_lo(i) + delnh3_hi(i)
            end if
            delouter        = max(delouter, abs(delnh3_hi(i) / rtbis(i)))
         end do
      end if

!  The convergence parameter "del" and the number of bisection
!  iterations "iter" are used to exit the outer loop further down in the code.

      jo = 0
      del = 1.0D0
!  Convergence check:  Have the activity coefficients
!  changed by more than eps2 in the last iteration?  If
!  so, repeat the calculation with the new activity
!  coefficients.  Keep track of the number of iterations.
      do while (del > eps2 .and. jo < itero)
!
         jo = jo + 1

!  Calculate liquid water, based on current ion concentration.

         call mach_hetv_water(ne11, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)

!  Important note:  if the dry particle initial conditions
!  are such that there's no initial ammonium nitrate, then
!  deliquescence should not take place, the water should be
!  close to zero, and many of the equations in the original
!  system will no longer be valid.  In HETV, this may
!  manifest itself as erroneously high water values via
!  HNO3 being used as a potential solute in the ZSR calculations.
!  This is corrected here; in the event that the initial (dry)
!  ammonium nitrate concentration is zero, the water will be
!  reset to its lower limit value of 1D-20

         where (amnitdry <= 0.D0)
            lwn = 1.D-20
         end where

         do i = 1, ne11
            lwo_lwn  = lwo(i) / lwn(i)
            nh3(i)   = nh3(i)   * lwo_lwn
            hno3(i)  = hno3(i)  * lwo_lwn
            so4(i)   = so4(i)   * lwo_lwn
            h(i)     = h(i)     * lwo_lwn
            no3(i)   = no3(i)   * lwo_lwn
            nh4(i)   = nh4(i)   * lwo_lwn
            hso4(i)  = hso4(i)  * lwo_lwn
            amsul(i) = amsul(i) * lwo_lwn
            nh3dry(i)   = nh3dry(i)  * lwo_lwn
            hno3dry(i)  = hno3dry(i) * lwo_lwn
            amnitdry(i) = amnitdry(i) * lwo_lwn
            amsuldry(i) = amsuldry(i) * lwo_lwn
            delnh3(i)   = delnh3(i)   * lwo_lwn
!
!  Update liquid water content
            lwo(i) = lwn(i)
         end do

!  Calculate activity coefficients

         call mach_hetv_activity(ne11, t, h, hso4, nh4, no3, so4, g_h_hso4, &
                                 g_h2_so4, g_h_no3, g_nh4_no3, g_nh42_so4)

         del = 0.D0
         do i = 1, ne11
            del = max(del, abs((g_h_hso4(i) - g_h_hso4_o(i)) / g_h_hso4(i)))
            g_h_hso4_o(i) = g_h_hso4(i)
            del = max(del, abs((g_h2_so4(i) - g_h2_so4_o(i)) / g_h2_so4(i)))
            g_h2_so4_o(i) = g_h2_so4(i)
            del = max(del, abs((g_h_no3(i) - g_h_no3_o(i)) / g_h_no3(i)))
            g_h_no3_o(i) = g_h_no3(i)
            del = max(del, abs((g_nh4_no3(i) - g_nh4_no3_o(i)) / g_nh4_no3(i)))
            g_nh4_no3_o(i) = g_nh4_no3(i)
            del = max(del, abs((g_nh42_so4(i) - g_nh42_so4_o(i)) / g_nh42_so4(i)))
            g_nh42_so4_o(i) = g_nh42_so4(i)
          end do

!  Factor activity coefficients into equilibrium constants
!  to simplify calculations.
!  *  Conversion factors:
!  *
!  *        atmospheres -> moles/kg H2O:
!  *
!  *   atmospheres *1.01325E5(kg m^-1 s^-2 /atm) /
!  *    { 8.3144 (kg m^2 s^-2 mol^-1 K^-1) * T (K) * Lw (kg H2O/kg air)
!  *          * rho (kg air/m3 air)  }    =  moles/(kg H2O)
!  *
!  *        atmospheres^-1 -> (moles/kg H2O)^-1 :
!  *
!  *         1          8.3144 * T * Lw * rho     kg H2O
!  *   ------------- * ----------------------  = --------
!  *   atmospheres      1.01325E5                 moles
!  *

!   Solve inner system of equations.
!
!  The value of delNH3 is varied from 0 to NH3dry in a multi-stage
!  binary search.  In the first stage, system (1, 2 & 3) is solved with
!  the given value of delNH3 as a constant.  System (1) determines the
!  amount of ammonium sulphate leaving the crystalline state with
!  the addition of water (the resulting NH4 is partitioned between
!  NH4+ and NH3).  System (2) determines the amount of nitrate being
!  transferred to/from HNO3 following the addition of water.  The
!  dry phase ammonium nitrate is assumed to be totally deliquesed.
!  System (3) determines the H ion concentration required to maintain
!  a charge balance with the ion concentrations determined in (1)
!  and (2).  Note that HSO4 is not part of the system of equations
!  at this stage; HSO4 concentrations are determined using a later
!  adjustment in system (5).  The root of the (HSO4-less) equations
!  is tested using system (4):  the overall result is a search for the
!  root of the NH4 <=> NH3 equilibrium equation as a function of the
!  available NH3.  Once this system has converged, the concentration
!  of HSO4 (and adjustments to H and SO4) are calculated.

         iwe = 0
         do i = 1, ne11
! 7:  (NH4)2SO4(s) <=> 2 NH4+ + SO4=
!  Divide by the cube of the ammonium sulphate activity coefficient:
            kamsul_m = kamsul(i) / (g_nh42_so4(i)**2 * g_nh42_so4(i))
!
            nh4available(i) = amnitdry(i) + delnh3(i)
            a  = nh4available(i)                             !*lwn(i)
            b  = 0.25D0 * nh4available(i) * nh4available(i)  !*lwn(i)*lwn(i)
            c  = - 0.25D0 * kamsul_m                         !*lwn(i)*lwn(i)*lwn(i)

!  Determine the smallest positive real root of the cubic polynomial
!  delnh42so4^3
!              +  nh4available * delnh42so4^2
!                        + 0.25 * (nh4available)^2 * delnh42so4
!                                              - 0.25 * knh42so4 = 0
!
!  Note:  if the liquid water content reaches the minimum value, then
!  the condensed phase has iterated to a waterless result.  Some of the
!  parts of the poly3 solver will have low number limit difficulties
!  with this situation.  In order to avoid these problems, poly3 is called
!  only for those points with water larger than the minimum limit, and
!  the delamsul value for the low water cases is taken to be zero.

            if (lwn(i) > lwmin) then
               iwe = iwe + 1
               ind2(iwe) = i
               a2(iwe) = a
               b2(iwe) = b
               c2(iwe) = c
            end if
         end do
         if (iwe /= 0) then
            call mach_hetv_poly3v(a2, b2, c2, delamsul2, islv2, iwe)
            if (chm_error_l) return
         end if
         do i = 1, ne11
            delamsul(i) = 0.D0
            islv(i) = 1
         end do
         do i = 1, iwe
            delamsul(ind2(i)) = delamsul2(i)
            islv(ind2(i)) = islv2(i)
         end do

         do i = 1, ne11
!           delamsul(i) = delamsul(i)   !/lwn(i)

!  Note:
!  (1) If islv(i) = 1, then the equation does not have a positive
!     real root.  If this is the case, no ammonium sulphate dissolves.
            delamsul(i) = dble(1 - islv(i)) * delamsul(i)
!
!  (2) The (NH4)2SO4 that may be lost from the crystalline
!     phase calculated in poly3  (delamsul) is the amount that would
!     have to be lost from the crystalline phase in order for
!     kamsul - (nh4)^2 * so4 = 0 to be maintained.  If this amount is
!     greater than the dry phase ammonium sulphate, then the ammonium
!     sulphate is assumed to be totally dissolved; delamsul = amsuldry,
!     and amsul = 0.
!
            delamsul(i) = min(delamsul(i), amsuldry(i))

!  Determine so4, nh4, nh3, (nh4)2so4 values (rest of system 1):
            so4(i) = delamsul(i)
            nh4(i) = max(0.D0, amnitdry(i) + delnh3(i) + 2.D0 * delamsul(i))
            amsul(i) = max(0.D0, amsuldry(i) - delamsul(i))
            nh3(i) = max(0.D0, nh3dry(i) - delnh3(i))
         end do

!  Note:  if there's no (NH4)2SO4 dissociating, and no dry NH4NO3
!  to dissociate, then the gas-phase concentration of NH3 should be
!  the dry value.  i.e. if delamsul=0, indicating that
!  no root could be found for the (NH4)2SO4 <-> 2 NH4 + SO4 equation,
!  and nh4no3dry=0, reset NH3 to the dry value, and
!  ensure that NH4+ is zero.

         do i = 1, ne11
            if (delamsul(i) == 0.0d0 .and. amnitdry(i) == 0.0d0) then
               nh3(i) = nh3dry(i)
               nh4(i) = 0.0D0
            end if
         end do

!  Note that kh2o = kw in the following, and that the aw values
!  cancel out.
         do i = 1, ne11
! 1, 2:  HNO3(g_eq) <=> HNO3(aq) <=> H+ + NO3-
! k1: Original units (moles HNO3(aq)/kg H2O)/(atm HNO3(g)), but
! gas is already in moles/(kg H2O); multiply by (R T Lwn rho).
! k2: Divide by the square of the HNO3 activity coefficient
!
            khno3_m = khno3(i) * lwn(i) / (g_h_no3(i) * g_h_no3(i))
!
! 3:  NH3(g_eq) <=> NH3(aq), NH3(aq) + H2O <=> NH4+ + OH-
!  ==> NH3(g) + H2O <=> NH4+ + OH-
!
! k3: Original units (moles NH3(aq)kg H2O)/(atm NH3(g)), but
! gas is already in moles/(kg H2O); multiply by (R T Lwn rho).
! k4: Multiply by the square of the HNO3 activity coefficient,
!  and divide by (the square of the NH4NO3 activity coefficient
!  times Kh2o (=Kw).
!  i.e. k4=(NH4)(OH)/( (NH3)(H2O) ), kH2O = (H)(OH)/(H2O)
!  Therefore,
!      k4/kH2O = (NH4)/( (NH3) (H) ) .  Then multiply by
!  gamma(H)/gamma(NO3) to correct for non-equilibrium.
!  Rather than divide by kH2O here, the kH2O value
!  is carried to the denominator in the relevant equations
            knh3_m = knh3(i) * lwn(i) * (g_h_no3(i)**2) / (g_nh4_no3(i)**2)
!
!  Correction, Nov 17, 2008:  if the NH4+ = 0, then
!  the ammonia equilibrium can no longer be used, and
!  the delno3 equation is invalid.  P.A. Makar/J. Murphy
            if (nh4(i) /= 0.D0) then
!  Non-zero NH4+, use delno3 equation:
               delno3(i) = min(max((khno3(i) * knh3_m * nh3(i) * hno3dry(i) - &
                      kh2o(i) * nh4(i) * amnitdry(i)) /                  &
                      (kh2o(i) * nh4(i) + khno3_m * knh3_m * nh3(i)),    &
                      0.0d0), hno3dry(i))
!
               hno3(i) = hno3dry(i) - delno3(i)
               no3(i) = amnitdry(i) + delno3(i)
            else
!  NH4+ = 0; set HNO3 = HNO3dry, and NO3- to zero, for now.
               hno3(i) = hno3dry(i)
               no3(i) = 0.D0
            end if

!  Determine the concentration of H required to maintain
!  charge balance.  Note that c(i) is defined negative, therefore
!  the root is always real and larger in magnitude than b(i).
!  If b(i) is negative, the positive root results in a positive H;
!  if b(i) is positive, the positive root results in a positive H.
!  The negative root will result in a negative value of H, regardless
!  of the sign of b(i).  Positive root for the quadratic is therefore
!  used.
!
!  Note that a taylor expansion of the square root is used in the
!  instance that round-off error can mess up the root.
!  "kh2o" = kw in the following.  Note that gamma H+ x gamma OH-
!  = 1 is assumed.
            b = nh4(i) - 2.D0 * so4(i) - no3(i)
            c = - kh2o(i) * aw(i)

            if (b /= 0.D0) then
               d = c / (b * b)
               v = 4.d0 * d
            else
               v = 1.D+03
            end if

            if (abs(v) <=  smrt .and. b > 0.D0) then
               h(i) =  - ((((14.D0 * d + 5.D0) * d + 2.D0) *  &
                          d + 1.D0) * d + 1.D0) * c / b
            else
               h(i) = 0.5d0 * (-b + sqrt(b * b - 4.d0 * c) )
            end if
            h(i) = max(h(i), 0.D0)
!
!   Activity coefficients changed by less than eps2 during
!   the previous iteration; the first part of the problem, for
!   the lower limit of the variable (delnh3), has converged.
!   Evaluate the root (system 4).
            diff(i) = knh3_m * nh3(i) * h(i) - nh4(i) * kh2o(i)
         end do

      end do ! End inner loop

      if (rooteval == 1 .or. rootsrch == 1) then
         do i = 1, ne11
            diff_lo(i) = diff(i)
         end do
      end if
!
!  Once the first midpoint has been calculated, use the
!  convergence checks from the last iteration do decide on
!  whether or not the set of problems have converged.  If
!  not, do another iteration.  If so, continue on and fill in
!  the remaining terms in the system of equations.
   end do
!
!  Final stage is an adjustment of the so4 values in order to
!  allow for the possible existence of hso4.  Note the b(i) is
!  defined negative, so -b(i) is positive, and the (non-complex)
!  quantity within the root will be a number smaller in magnitude
!  than b(i).  The sign of the root must therefore be positive
!  in order to get the smaller of the two roots.
   do i = 1, ne11
! 5:  HSO4- <=> SO4= + H+
!  Multiply by the square of the H-HSO4 activity coefficient,
!  and divide by the cube of the H2SO4 activity coefficient:
      khso4_m = khso4(i) * g_h_hso4(i)**2 / (g_h2_so4(i) * g_h2_so4(i)**2)
!
      b = - (h(i) + so4(i) + khso4_m)
      c = h(i) * so4(i)
      root = b * b - 4.D0 * c
!
!  Modification:  need to add condition that if there's no SO4
!  to convert into HSO4, you don't bother to take the root:
      if (root >= 0.D0 .and. so4(i) > 0.D0) then
         hso4(i) = 0.5d0 * ( - b - sqrt(root) )
      else
         hso4(i) = 0.D0
      end if
      hso4(i) = max(hso4(i), 0.0D0)

      h(i) = max(h(i) - hso4(i), 0.0D0)
      so4(i) = max(so4(i) - hso4(i), 0.0D0)
   end do

   do i = 1, ne11
!  Convert values back to moles/kg air from moles/kg water by
!  multiplying by the liquid water content:
      nh3(i)   = nh3(i)   * lwn(i)
      hno3(i)  = hno3(i)  * lwn(i)
      so4(i)   = so4(i)   * lwn(i)
      h(i)     = h(i)     * lwn(i)
      no3(i)   = no3(i)   * lwn(i)
      nh4(i)   = nh4(i)   * lwn(i)
      hso4(i)  = hso4(i)  * lwn(i)
      amsul(i) = amsul(i) * lwn(i)
   end do
!
! For those cases in which the ammonium ion is zero, some
! nitric acid may still dissolve via the HNO3 <-> NO3-(aq) + H+(aq)
! equilibrium.  Those cases are extracted here from the entire
! set that are passed to this case, and are corrected for the
! HNO3 concentration accordingly.
!
   nh4_counter = 0
   ncas11_cor  = 0
   do i = 1, ne11
      if (nh4(i) == 0) then
         nh4_counter = nh4_counter + 1
         ncas11_cor(i) = 11
      end if
   end do
!
   if (nh4_counter > 0) then
      call mach_hetv_corrhno3(ne11, nr, nh4_counter, 11, so4, no3, nh4, hso4, &
                              hno3, h, lwn, t, aw, rho, k0, p1, p2, ncas11_cor)
   end if

   so4_i   = unpack(so4,   ncas == nc, so4_i)
   nh4_i   = unpack(nh4,   ncas == nc, nh4_i)
   no3_i   = unpack(no3,   ncas == nc, no3_i)
   hso4_i  = unpack(hso4,  ncas == nc, hso4_i)
   hno3_i  = unpack(hno3,  ncas == nc, hno3_i)
   nh3_i   = unpack(nh3,   ncas == nc, nh3_i)
   h_i     = unpack(h,     ncas == nc, h_i)
   lwn_i   = unpack(lwn,   ncas == nc, lwn_i)
   amsul_i = unpack(amsul, ncas == nc, amsul_i)

   return
end subroutine mach_hetv_case11
