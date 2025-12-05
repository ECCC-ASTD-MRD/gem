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
! Fichier/File   : mach_hetv_case8.ftn90
! Creation       : P. Makar,  V. Bouchet,  A. Nenes,  S. Gravel,  B. Pabla,  S. Menard
! Description    : Heterogeneous chemistry solver for case8.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne8 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case 8's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne8.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                        1.5 <= TA/TS < 2.0, drh_leto <= rh < drh_ambis
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!                  (1)
!                     HSO4 <=> H + SO4 , kHSO4
!                     TA - NH4 - 2 (NH4)2SO4 = 0
!                     TS - SO4 - HSO4 - (NH4)2SO4 = 0
!                     H + NH4 = 2 SO4 + HSO4
!                   (2)
!                      (NH4)2SO4 <=> 2 NH4 + SO4, k(NH4)2SO4
!                   The system of equations to be solved in this case are as follows:
!                      Equations
!                   (1)
!                      kHSO4 (HSO4)- (H)(SO4) = 0
!                      TA - NH4 - (NH4)2SO4 = 0,
!                      TS - SO4 - HSO4 - (NH4)2SO4 = 0,
!                      H + NH4 - 2 SO4 - HSO4 = 0.
!                   (2)
!                      k(NH4)2SO4 - (NH4)^2 * (SO4) = 0
!                   The solution to the system of equations:
!
!                   (1)
!                      Let
!                      b = [kHSO4 - TA + (NH4)2SO4 + TS],
!                      c = kHSO4 [(NH4)2SO4 - TS],
!                      SO4 = 0.5*(-b + sqrt(b**2 -4c)),
!                      NH4 = TA - 2 (NH4)2SO4
!                      H = - TA + (NH4)2SO4 + SO4 + TS
!                      HSO4 = TS - SO4 - (NH4)2SO4
!                      (2)
!                      k(NH4)2SO4 - (NH4)^2 * (SO4) = 0
!                   Variable, Range
!                      (NH4)2SO4, 0 to 2TA - 3TS
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case8(npts, nr, ne8, nc, so4_i, no3_i, nh4_i, hso4_i, &
                           hno3_i, h_i, nh3_i, amsul_i, lwn_i, t_i, rh_i,  &
                           ts_i, ta_i, tn_i, k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_water, mach_hetv_activity
   use mach_hetv_mod,         only: tstd, small, smrt, iter, itero, ndiv, eps, eps2
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne8
   integer(kind=4), intent   (in) :: nc
   integer(kind=4), intent   (in) :: ncas    (npts)
   real(kind=8),    intent   (in) :: k0      (nr)
   real(kind=8),    intent   (in) :: p1      (nr)
   real(kind=8),    intent   (in) :: p2      (nr)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
   real(kind=8),    intent(inout) :: amsul_i (npts)
   real(kind=8),    intent(inout) :: lwn_i   (npts)
   real(kind=8),    intent   (in) :: ts_i    (npts)
   real(kind=8),    intent   (in) :: ta_i    (npts)
   real(kind=8),    intent   (in) :: tn_i    (npts)
   real(kind=8),    intent   (in) :: t_i     (npts)
   real(kind=8),    intent   (in) :: rh_i    (npts)
!!if_off
!
!  Local variables:
!
   real(kind=8), dimension(ne8) :: so4, no3, nh4, hso4, hno3, h, lwn
   real(kind=8), dimension(ne8) :: nh3, amsul, lwo
   real(kind=8), dimension(ne8) :: t, aw
   real(kind=8), dimension(ne8) :: ta, ts, tn
!
   integer(kind=4), dimension(ne8) :: rflag, rootloc, bflag, rootmin
   integer(kind=4)                 :: i, jo, rooteval, rootsrch, jsum
   real(kind=8)                    :: delouter, del, b, c, d, v, dndiv
   real(kind=8), dimension(ne8)    :: g_h_hso4, g_h_hso4_o, &
                                      g_h2_so4, g_h2_so4_o, &
                                      g_nh42_so4, g_nh42_so4_o
   real(kind=8), dimension(ne8)    :: awu, law

   real(kind=8), dimension(ne8)    :: amsul_lo, amsul_hi, choice, rtbis, dx, diff, &
                                      diff_lo, difflo, diff_1, locmin, avdiff
!  arrays for reaction rates:
   real(kind=8), dimension(ne8) :: khso4, kamsul
   real(kind=8) :: work1, work2, khso4_m, kamsul_m, lwo_lwn

!   so4   = pack(so4_i,  ncas == nc)
!   nh4   = pack(nh4_i,  ncas == nc)
   no3   = pack(no3_i,  ncas == nc)
!   hso4  = pack(hso4_i, ncas == nc)
   hno3  = pack(hno3_i, ncas == nc)
!   h     = pack(h_i,    ncas == nc)
!   nh3   = pack(nh3_i,  ncas == nc)
!   amsul = pack(amsul_i, ncas == nc)
   ta    = pack(ta_i,   ncas == nc)
   ts    = pack(ts_i,   ncas == nc)
   tn    = pack(tn_i,   ncas == nc)
   aw    = pack(rh_i,   ncas == nc)
   t     = pack(t_i,    ncas == nc)
!   lwn   = pack(lwn_i,  ncas == nc)

!  Set value of subdivision for root search:
   dndiv = 1.0D0 / dble(ndiv)

   do i = 1, ne8
!  Initial guess for components:  TS split evenly between SO4 and HSO4,
!  NH4 = TA, and H = max(0, 2 SO4 + HSO4 - NH4), (NH4)3H(SO4)2 = 0
      so4(i)   = 0.5 * ts(i)
      hso4(i)  = 0.5 * ts(i)
      nh4(i)   = ta(i)
      h(i)     = max(0.d0, 2.d0 * so4(i) + hso4(i) - nh4(i))
      amsul(i) = 0.d0

      lwn(i) = 1.d0
!
!  Set value of subdivision for root search:
      dx(i) = max(0.d0, (2.d0 * ta(i) - 3.d0 * ts(i)) * dndiv)
!      
!  Calculate rates for Equilibrium constants
      work1 = tstd / t(i) - 1.D0
      work2 = 1.d0 + log(tstd / t(i)) - tstd / t(i)
!  5:  HSO4- <=> SO4= + H+
      khso4(i) = k0(5) * exp(p1(5) * work1 + p2(5) * work2)
!  7:  (NH4)2SO4(s) <=> 2 NH4+ + SO4=
      kamsul(i) = k0(7) * exp(p1(7) * work1 + p2(7) * work2)
!
!  polynomials for electrolytes only good between aw=0.002 and aw=0.98
      awu(i) = min(max(aw(i), 2.0d-03), 0.98d0)
      law(i) = dble(log(real(awu(i))))    !double precision log unnec.
   end do

!  Looping point for outer iteration 
!  (subdivision of function interval and search for roots):
   rooteval = 0
   rootsrch = 0
   delouter = 1.0d0

!  Repeat the calculations, for all ndiv+1 initial function evaluations:
   do while ((rooteval <= ndiv + 1) .or. (rootsrch <= 3) .or. &
             (delouter > eps .and. rootsrch < iter))

!   Initial guess for liquid water content (will be replaced in calculations to follow):
      do i = 1, ne8
         lwo(i) = 1.d0

!   Reset concentrations to starting liquid water content
         ta(i)    = ta(i)    * lwn(i)
         ts(i)    = ts(i)    * lwn(i)
         tn(i)    = tn(i)    * lwn(i)
         hno3(i)  = hno3(i)  * lwn(i)
         so4(i)   = 0.5d0 * ts(i)     !so4(i)   * lwn(i)
!         h(i)     = h(i)     * lwn(i)
         no3(i)   = 0.d0            !no3(i)   * lwn(i)
         nh4(i)   = ta(i)           !nh4(i)   * lwn(i)
         hso4(i)  = 0.5d0 * ts(i)     !hso4(i)  * lwn(i)
         amsul(i) = amsul(i) * lwn(i)
         h(i) = max(0.d0, 2.d0 * so4(i) + hso4(i) - nh4(i))

!  Zero the 'old value' activity coefficient arrays
         g_h_hso4_o(i)     = 0.D0
         g_h2_so4_o(i)     = 0.D0
         g_nh42_so4_o(i) = 0.d0
      end do

      rooteval = rooteval + 1

!  First root evaluation is for upper bound, amsul = maximum (total dissolution)
      if (rooteval == 1) then
         do i = 1, ne8
            amsul(i) = max((2.d0 * ta(i) - 3.d0 * ts(i) ), 0.d0)
         end do
      end if

      if (rooteval == 2) then
!  Store results of first root evaluation:
         do i = 1, ne8
            difflo(i) = diff_lo(i)
            diff_1(i) = diff_lo(i)
         end do

!  Set the indirect address for the location of the sign change (rootloc) and the counter
!  for the number of sign changes (rflag) for the given gridpoint to zero:
         do i = 1, ne8
            rootloc(i) = 0
            rflag(i)   = 0
            rootmin(i) = 0
            locmin(i)  = 1.d30
         end do
      end if

!  Subsequent root evaluations (but before the binary search
!  refinement):  determine the value of the variable to be
!  used for the next function evaluation:
      if (rooteval > 1 .and. rooteval <=   ndiv + 1) then
         do i = 1, ne8
            amsul(i) = max((2.d0 * ta(i) - 3.d0 * ts(i) ) - dx(i) * &
                           dble(rooteval - 1), 0.d0)
         end do
      end if

!  Save results of previous function evaluations,
!  and check for sign change of function:
      if (rooteval > 2 .and. rooteval <=   ndiv + 2) then
!
!  Determine whether the boundaries of the current sub interval
!  constitute the minimum for the entire interval.  If so,
!  store the right-hand index:
         do i = 1, ne8
            avdiff(i) = 0.5D0 * abs(diff(i) + diff_1(i))
            locmin(i) = min(locmin(i), avdiff(i))
         end do
         do i = 1, ne8
            if (avdiff(i) <= locmin(i)) rootmin(i) = rooteval - 1
         end do

!  Determine whether a sign change occurred in the last
!  interval:
         do i = 1, ne8
            if (sign(1.D0, diff_1(i)) * sign(1.D0, diff(i)) < 0.D0) then
               rflag(i) = rflag(i) + 1
            end if
         end do

!  If it was the first sign change detected for the given
!  gridpoint, save its location:
         do i = 1, ne8
            if (rflag(i) == 1) then
               rootloc(i) = rooteval - 1
               rflag(i) = 2
            end if
         end do

!  Update the evaluated function arrays:
         do i = 1, ne8
            diff_1(i) = diff(i)
         end do
      end if

!  The next section determines the boundaries to be used in the
!  binary search refinement of the given root.  Sometimes both limits
!  of the search are set to the same value (upper or lower boundary
!  of the entire interval), if no sign change occurred in the
!  entire interval.
      if (rooteval == ndiv + 2) then
!  Set the "boundary has been chosen" flag to "off":
         do i = 1, ne8
            bflag(i) = 0
         end do

! following block changed (see paul's e-mail dated Jan 15/03)
!  Use the lower limit of the variable space, if that lower limit
!  corresponds to a root and a sign change does not occur elsewhere:
         do i = 1, ne8
            if (abs(difflo(i)) <=   eps .and. rootloc(i) == 0) then
               amsul_lo(i) = small
               amsul_hi(i) = small
               bflag(i)  = 1
            end if
         end do

!  If a root was found in the variable space for the function at the
!  given gridpoint, use the stored index to give the boundaries:
         do i = 1, ne8
            if (rootloc(i) /= 0 .and. bflag(i) == 0) then
               amsul_lo(i) = max((2.d0 * ta(i) - 3.d0 * ts(i) ) - dx(i) * dble(rootloc(i) - 2), 0.d0)
               amsul_hi(i) = max(amsul_lo(i) - dx(i), 0.d0)
               bflag(i) = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  but the upper limit is a root, use that value for the root:
         do i = 1, ne8
            if (bflag(i) == 0 .and. abs(diff(i)) <=  eps .and. rootloc(i) == 0) then
               amsul_lo(i) = amsul(i)
               amsul_hi(i) = amsul(i)
               bflag(i)  = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  and all the values of the root evaluations were positive, then
!  (NH4)^2 (SO4)/ K[(NH4)2SO4] - 1 > 0 always - this implies that
!  there will always be higher ionic concentrations (and insufficient
!  (NH4)2SO4) to satisfy the (NH4)2SO4 formation equation - while
!  local minima (not roots) may exist, the chemically correct course
!  in this case is to assume that (NH4)2SO4 is at the maximum
!  value.  In order to ensure that this is done, both boundaries of the
!  root search are set to the upper value of (NH4)2SO4:
         do i = 1, ne8
            if (bflag(i) == 0 .and. diff(i) > 0.0d0) then
               amsul_lo(i) = amsul(i)
               amsul_hi(i) = amsul(i)
               bflag(i) = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  and all the values of the root evaluations were negative, then
!  (NH4)^2 (SO4)/ K[(NH4)2SO4] -1 < 0 always - this implies that
!  there will always be insufficient ionic concentrations (and an excess
!  of (NH4)2SO4) to satisfy the (NH4)2SO4 formation equation - while
!  local minima (not roots) may exist, the chemically correct course
!  in this case is to assume that (NH4)2SO4 is at the minimum
!  value.  In order to ensure that this is done, both boundaries of the
!  root search are set to the minimum value of (NH4)2SO4:
         do i = 1, ne8
            if (bflag(i) == 0 .and. diff(i) < 0.0d0) then
               amsul_lo(i) = small
               amsul_hi(i) = small
               bflag(i) = 1
            end if
         end do

!  Check to make sure that some boundaries for the binary search were found:
         jsum = sum(bflag)
         if (jsum < ne8) then
            write(*, *) '> WARNING <'
            write(*, *) '> Root interval not set in Case 8'
            write(*, *) '> Condition occurred at the following gridpoint:'
            do i = 1, ne8
               if (bflag(i) == 0) write(*, *) '> ', i
            end do
         end if
      end if

!  If the initial subdivision of the variable space has been completed,
!  start the bisection search refinement of the chosen function interval
      if (rooteval >=  ndiv + 2) rootsrch = rootsrch + 1

!  Reevaluate the function on the lhs of the interval:
      if (rootsrch == 1) then
         do i = 1, ne8
            amsul(i) = amsul_lo(i)
         end do
      end if

!  Reevaluate the function on the rhs of the interval:
      if (rootsrch == 2) then
         do i = 1, ne8
            amsul(i) = amsul_hi(i)
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
         do i = 1, ne8
!  "choice" is -1 if diff_lo < 0, 1 otherwise:
            choice(i) = sign(1.E+0, real(diff_lo(i)))
            dx(i) = choice(i) * (amsul_lo(i) - amsul_hi(i))
!  "choice" is 1 if diff_lo < 0, 0 otherwise
            choice(i) = 0.5D+0 * (1.D+0 - choice(i))
            rtbis(i) = choice(i) * amsul_lo(i) + (1.D+0 - choice(i)) * amsul_hi(i)
            amsul_hi(i) = 0.5d0 * dx(i)
            amsul_lo(i) = rtbis(i)
            amsul(i) = max(amsul_lo(i) + amsul_hi(i), 0.D0)
!  rtbis is used to store the value of dx for future convergence calcs.
            rtbis(i) = max(abs(dx(i)), small)
         end do
      end if

!  All subsequent root evaluations:  determine which side of the interval contains
!  the root crossing, rearrange the bounds, and subdivide the interval again:

      delouter = 0.D0
      if (rootsrch > 3) then
         do i = 1, ne8
!  "choice" is 1 if diff<0, 0 otherwise:
            choice(i)   = 0.5D+0 * (1.D+0 - sign (1.E+0, real(diff(i)) )   )
            amsul_lo(i) = amsul(i) * choice(i) + (1.D+0 - choice(i)) * amsul_lo(i)
            delouter    = max( delouter, abs(amsul_hi(i) / rtbis(i)))
            amsul_hi(i) = amsul_hi(i) * 0.5d0
            amsul(i)    = max(amsul_lo(i) + amsul_hi(i), 0.D0)
         end do
      end if

!  The convergence parameter "del" and the number of bisection iterations "iter"
!  are used to exit the outer loop further down in the code.

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
         call mach_hetv_water(ne8, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)

!  Important note:  if the dry particle initial conditions are such that there's no initial
!  letovicite, thendeliquescence should not take place, the water should be close to zero,
!  and many of the equations in the original system will no longer be valid.  In HETV, this
!  may manifest itself as erroneously high water values via acids being used as a potential
!  solute in the ZSR calculations. This is corrected here; in the event that the initial (dry)
!  letovicite concentration is zero, the water will be reset to its lower limit value of 1D-20
         do i = 1, ne8
            if (2.d0 * ts(i) - ta(i) <= 0.d0) lwn(i) = 1.d-20
         end do

!  (2)  Calculate new concentrations of ions and other phases; 
!       multiply old concentration by ratio of old to new liquid water contents
         do i = 1, ne8
            lwo_lwn  = lwo(i) / lwn(i)
            ta(i)    = ta(i)    * lwo_lwn
            ts(i)    = ts(i)    * lwo_lwn
            tn(i)    = tn(i)    * lwo_lwn
            hno3(i)  = hno3(i)  * lwo_lwn
            so4(i)   = so4(i)   * lwo_lwn
            h(i)     = h(i)     * lwo_lwn
            no3(i)   = no3(i)   * lwo_lwn
            nh4(i)   = nh4(i)   * lwo_lwn
            hso4(i)  = hso4(i)  * lwo_lwn
            amsul(i) = amsul(i) * lwo_lwn
!  Update liquid water content
            lwo(i) = lwn(i)
         end do

!  Calculate activity coefficients
         call mach_hetv_activity(ne8, t, h, hso4, nh4, no3, so4, g_h_hso4, &
                                 g_h2_so4, g_nh42_so4=g_nh42_so4)

!  Calculate maximum change in activity coefficients since the previous call
!  in this subroutine:

         del = 0.D0
         do i = 1, ne8
            del = max(del, abs( (g_h_hso4(i) - g_h_hso4_o(i)) / g_h_hso4(i)))
            g_h_hso4_o(i) = g_h_hso4(i)
            del = max(del, abs( (g_h2_so4(i) - g_h2_so4_o(i)) / g_h2_so4(i)))
            g_h2_so4_o(i) = g_h2_so4(i)
            del = max(del, abs( (g_nh42_so4(i) - g_nh42_so4_o(i)) / g_nh42_so4(i) ) )
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

         do i = 1, ne8
! 5:  HSO4- <=> SO4= + H+
!  Multiply by the square of the H-HSO4 activity coefficient,
!  and divide by the cube of the H2SO4 activity coefficient:
            khso4_m = khso4(i) * (g_h_hso4(i)**2) / (g_h2_so4(i) * g_h2_so4(i)**2)
! 7:  (NH4)2SO4(s) <=> 2 NH4+ + SO4=
!  Divide by the cube of the ammonium sulphate activity coefficient:
            kamsul_m = kamsul(i) / (g_nh42_so4(i) * g_nh42_so4(i) * g_nh42_so4(i))

!   Solve inner system of equations.
!  The value of (NH4)2SO4 is varied from 0 to 2 TA - 3 TS in a two
!  stage binary search.  In the first stage, system (1) is solved
!  with the given value of (NH4)2SO4 as a constant.  The resulting
!  SO4 and NH4 value is tested using system (2):  the overall result is a
!  search for the root of the lammonium sulphate equation as a function of the
!  available ammonium sulphate.  Note that 1.5 <= TA/TS <2 if this
!  subroutine is called; the value within the root is always positive,
!  since "c" is <= 0.
!  In order to prevent possible round off errors in the SO4 determination,
!  a Taylor series approximation is used in the event that |c| << |b|,
!  and b > 0.
            b = khso4_m - ta(i) + amsul(i) + ts(i)
            c = khso4_m * (amsul(i)  - ts(i))
            if (b /= 0.d0) then
               d = c / (b * b)
               v = 4.d0 * d
            else
               v = 1.d+03
            end if
            if (abs(v) <= smrt .and. b > 0.d0) then
               so4(i) = - ((((14.d0 * d + 5.d0) * d + 2.d0) * d + 1.d0) * &
                              d + 1.d0) * c / b
            else
               so4(i) = 0.5d0 * (-b + sqrt(b * b - 4.d0 * c))
            end if
            so4(i)  = max(so4(i), 0.D0)
            nh4(i)  = max(ta(i) - 2.d0 * amsul(i), 0.d0)
            h(i)    = max(ts(i) + so4(i) + amsul(i) - ta(i), 0.d0)
            hso4(i) = max(ts(i) - so4(i) - amsul(i), 0.d0)
      
            diff(i) = (nh4(i) * nh4(i) * so4(i)) / kamsul_m - 1.d0
         end do

      end do ! End inner loop

!   Activity coefficients changed by less than eps2 during
!   the previous iteration; the first part of the problem, for
!   the lower limit of the variable (nh4), has converged.
!   Evaluate the root for the final equation to get the lower
!   bound value of the root.
      if (rooteval == 1 .or. rootsrch == 1) then
         do i = 1, ne8
            diff_lo(i) = diff(i)
         end do
      end if

!  Once the first midpoint has been calculated, use the
!  convergence checks from the last iteration do decide on
!  whether or not the set of problems have converged.  If
!  not, do another iteration.  If so, continue on and fill in
!  the remaining terms in the system of equations.
   end do
   
   do i = 1, ne8
      hno3(i)  = tn(i)
      no3(i)   = 0.d0
      nh3(i)   = 0.d0

!  Convert values back to moles/kg air from moles/kg water by multiplying by the liquid water content:
      hno3(i)  = hno3(i)  * lwn(i)
      so4(i)   = so4(i)   * lwn(i)
      h(i)     = h(i)     * lwn(i)
      nh4(i)   = nh4(i)   * lwn(i)
      hso4(i)  = hso4(i)  * lwn(i)
      amsul(i) = amsul(i) * lwn(i)
   end do

   so4_i   = unpack(so4,   ncas == nc, so4_i)
   no3_i   = unpack(no3,   ncas == nc, no3_i)
   nh4_i   = unpack(nh4,   ncas == nc, nh4_i)
   hso4_i  = unpack(hso4,  ncas == nc, hso4_i)
   hno3_i  = unpack(hno3,  ncas == nc, hno3_i)
   h_i     = unpack(h,     ncas == nc, h_i)
   nh3_i   = unpack(nh3,   ncas == nc, nh3_i)
   amsul_i = unpack(amsul, ncas == nc, amsul_i)
   lwn_i   = unpack(lwn,   ncas == nc, lwn_i)

   return
end subroutine mach_hetv_case8