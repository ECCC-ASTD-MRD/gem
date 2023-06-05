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
! Fichier/File   : mach_hetv_case1.ftn90
! Creation       : P. Makar, V. Bouchet, A. Nenes, S. Gravel, B. Pabla, S. Menard
! Description    : Heterogeneous chemistry solver for case1.
!                  Based on algorithm's in Nenes et al's ISORROPIA solver,
!                  recoded to vectorize over the gridpoint dimension.
!                  All input arrays are of length npts (total number of gridpoints
!                  submitted from the calling code.  The subsection of each
!                  1-D array between index 1 and ne1 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case 1's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne1.
!
!                  Units on input are moles/kg air.
!
!                  Revised root search algorithm; uses subdivision of interval to
!                  find smallest root, followed by binary search (algorithm similar to ISORROPIA).
!
! Extra info :     Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                  TA/TS < 1, rh < drh_ambis
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!                  (1)
!                    HSO4 <=> H + SO4 , kHSO4
!                    NH4HSO4 = TA - NH4
!                    TS = SO4 + HSO4 + NH4HSO4
!                    H + NH4 = 2 SO4 + HSO4
!                  (2)
!                    NH4HSO4 <=> NH4 + HSO4, kNH4HSO4
!                  The system of equations to be solved in this case are as follows:
!                  Equations
!                  (1)
!                   kHSO4 (HSO4)- (H)(SO4) = 0
!                   TA - NH4 - NH4HSO4 = 0,
!                   TS - SO4 - HSO4 - NH4HSO4 = 0,
!                   H + NH4 - 2 SO4 - HSO4 = 0.
!                  (2)
!                    kNH4HSO4 - (NH4)(HSO4) = 0
!
!                  The solution to the system of equations:
!                  (1)
!                    Let
!                    b = kHSO4 - TA + TS,
!                    c = kHSO4 (TA - NH4 - TS)
!                    SO4 = 0.5*(-b + sqrt(b**2 -4c)),
!                    NH4HSO4 = TA - NH4,
!                    H = SO4 + TS - TA,
!                    HSO4 = TS - SO4 - TA + NH4
!                  (2)
!                    kNH4HSO4 - (NH4)(HSO4) = 0
!                  Variable, Range
!                    NH4, 0 to TA
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case1(npts, nr, ne1, so4_i, no3_i, nh4_i, hso4_i,    &
                           hno3_i, h_i, nh3_i, ambis_i, lwn_i, t_i, rh_i, &
                           ts_i, ta_i, tn_i, k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_water, mach_hetv_activity
   use mach_hetv_mod,         only: lwmin, tstd, small, smrt, iter, itero, ndiv, eps, eps2
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne1
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
   real(kind=8),    intent(inout) :: ambis_i (npts)
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
   real(kind=8), dimension(ne1) :: so4, no3, nh4, hso4, hno3, h, lwn
   real(kind=8), dimension(ne1) :: nh3, ambis
   real(kind=8), dimension(ne1) :: t, aw
   real(kind=8), dimension(ne1) :: ta, ts, tn
!
   integer(kind=4)              :: i, jo, rooteval, rootsrch, jsum
   real(kind=8)                 :: del, delouter, b, c, d, v
   real(kind=8)                 :: lwo_lwn, dndiv
!  Workspace arrays
   real(kind=8), dimension(ne1) :: g_h_hso4, g_h_hso4_o
   real(kind=8), dimension(ne1) :: g_h2_so4, g_h2_so4_o
   real(kind=8), dimension(ne1) :: g_nh4_hso4, g_nh4_hso4_o
   real(kind=8), dimension(ne1) :: lwo
   real(kind=8), dimension(ne1) :: awu, law
   real(kind=8), dimension(ne1) :: nh4_lo, nh4_hi
   real(kind=8), dimension(ne1) :: choice, rtbis, dx
   real(kind=8), dimension(ne1) :: diff_lo, diff
   real(kind=8), dimension(ne1) :: difflo, diff_1
!  arrays for reaction rates:
   real(kind=8), dimension(ne1) :: khso4, kambis
   real(kind=8) :: khso4_m, kambis_m, work1, work2

!  Local variables that may be used in the solution:
   Integer(kind=4), dimension(ne1) :: rflag, rootloc, bflag

!   so4   = pack(so4_i,  ncas == 1)
!   nh4   = pack(nh4_i,  ncas == 1)
   no3   = pack(no3_i,  ncas == 1)
!   hso4  = pack(hso4_i, ncas == 1)
   hno3  = pack(hno3_i, ncas == 1)
!   h     = pack(h_i,    ncas == 1)
!   ambis = pack(ambis_i, ncas == 1)
   ta    = pack(ta_i,   ncas == 1)
   ts    = pack(ts_i,   ncas == 1)
   tn    = pack(tn_i,   ncas == 1)
   aw    = pack(rh_i,   ncas == 1)
   t     = pack(t_i,    ncas == 1)
!   lwn   = pack(lwn_i, ncas == 1)

!  Set value of subdivision for root search:
   dndiv = 1.0D0 / dble(ndiv)
!
   do i = 1, ne1
!  Initial guess for components:  TS split evenly between SO4 and HSO4,
!  NH4 = TA, and H = max(0, 2 SO4 + HSO4 - NH4), NH4HSO4 = 0
      so4(i)   = 0.5d0 * ts(i)
      hso4(i)  = 0.5d0 * ts(i)
      nh4(i)   = ta(i)
      h(i)     = max(0.d0, 2.d0 * so4(i) + hso4(i) - nh4(i))
      ambis(i) = 0.d0
!
      lwn(i) = 1.d0
!
!  Set value of subdivision for root search:
      dx(i) = (ta(i) - small) * dndiv
!
!  Calculate rates for Equilibrium constants
      work1 = tstd / t(i) - 1.D0
      work2 = 1.d0 + log(tstd / t(i)) - tstd / t(i)
!  5:  HSO4- <=> SO4= + H+
      khso4(i) = k0(5) * exp(p1(5) * work1 + p2(5) * work2)
!  6:  NH4HSO4(s) <=> NH4+ + HSO4-
      kambis(i) = k0(6) * exp(p1(6) * work1 + p2(6) * work2)
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

      do i = 1, ne1
! Initial guess for liquid water content (will be replaced in calculations to follow:
         lwo(i) = 1.d0

!   Reset concentrations to starting liquid water content
         ta(i)    = ta(i)    * lwn(i)
         ts(i)    = ts(i)    * lwn(i)
         tn(i)    = tn(i)    * lwn(i)
         hno3(i)  = hno3(i)  * lwn(i)
         so4(i)   = so4(i)   * lwn(i)
         h(i)     = h(i)     * lwn(i)

         no3(i)   = no3(i)   * lwn(i)
         nh4(i)   = nh4(i)   * lwn(i)
         hso4(i)  = hso4(i)  * lwn(i)
         ambis(i) = ambis(i) * lwn(i)

!  Zero the 'old value' activity coefficient arrays
         g_h_hso4_o(i)   = 0.D0
         g_h2_so4_o(i)   = 0.D0
         g_nh4_hso4_o(i) = 0.D0
      end do

      rooteval = rooteval + 1

!  First root evaluation is for lower bound, nh4 = 0
      if (rooteval == 1) then
         do i = 1, ne1
            nh4(i) = small
         end do
      end if

      if (rooteval == 2) then
!  Store results of first root evaluation:
!  Set the indirect address for the location of the sign change (rootloc) and the
!  counter for the number of sign changes (rflag) for the given gridpoint to zero:
         do i = 1, ne1
            difflo(i) = diff_lo(i)
            diff_1(i) = diff_lo(i)
            rootloc(i) = 0
            rflag(i)   = 0
         end do
      end if

!  Subsequent root evaluations (but before the binary search refinement):  determine
!  the value of the variable to be used for the next function evaluation:
      if (rooteval > 1 .and. rooteval <= ndiv + 1) then
         do i = 1, ne1
            nh4(i) = small + dx(i) * dble(rooteval - 1)
         end do
      end if

!  Save results of previous function evaluations, and check for sign change of function:
      if (rooteval > 2 .and. rooteval <= ndiv+2) then

!  Determine whether a sign change occurred in the last interval:
         do i = 1, ne1
            if (sign(1.D0, diff_1(i))*sign(1.D0, diff(i)) < 0.D0) &
               rflag(i) = rflag(i) + 1
         end do

!  If it was the first sign change detected for the given
!  gridpoint, save its location:
         do i = 1, ne1
            if (rflag(i) == 1) then
               rootloc(i) = rooteval - 1
               rflag(i) = 2
            end if
         end do

!  Update the variable and evaluated function arrays:
         do i = 1, ne1
            diff_1(i) = diff(i)
         end do
      end if

!  The next section determines the boundaries to be used in the binary search refinement
!  of the given root.  Sometimes both limits of the search are set to the same value
!  (upper or lower boundary of the entire interval), if no sign change occurred in the
!  entire interval.
      if (rooteval ==  ndiv + 2) then

!  Set the "boundary has been chosen" flag to "off":
         do i = 1, ne1
            bflag(i) = 0
         end do

!  following block changed (see paul's e-mail dated Jan 15/03)
!  Use the lower limit of the variable space, if that lower limit
!  corresponds to a root and a sign change does not occur elsewhere:
         do i = 1, ne1
            if (abs(difflo(i)) <= eps .and. rootloc(i) == 0) then
               nh4_lo(i) = small
               nh4_hi(i) = small
               bflag(i)  = 1
            end if
         end do

!  If a root was found in the variable space for the function at the
!  given gridpoint, use the stored index to give the boundaries:
         do i = 1, ne1
            if (rootloc(i) /= 0 .and. bflag(i) ==  0) then
               nh4_lo(i) = small + dx(i) * dble(rootloc(i) - 2)
               nh4_hi(i) = nh4_lo(i) + dx(i)
               bflag(i) = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  but the upper limit is a root, use that value for the root:
         do i = 1, ne1
            if (bflag(i) == 0 .and. abs(diff(i)) <= eps .and. &
                rootloc(i) == 0) then
               nh4_lo(i) = nh4(i)
               nh4_hi(i) = nh4(i)
               bflag(i)  = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  and all the values of the root evaluations were positive, then
!  K(NH4HSO4) > (NH4)(HSO4) always - this implies that there will
!  never be enough NH4 to satisfy the ammonium bisulphate formation
!  equation - while local minima (not roots) may exist, the chemically
!  correct course in this case is to assume that NH4 is the maximum
!  value.  In order to ensure that this is done, both boundaries of the
!  root search are set to the upper value of NH4:
         do i = 1, ne1
            if (bflag(i) == 0 .and. diff(i) > 0.) then
               nh4_lo(i) = nh4(i)
               nh4_hi(i) = nh4(i)
               bflag(i) = 1
            end if
         end do

!  If no root was found in the variable space at the given gridpoint,
!  and all the values of the root evaluations were negative, then
!  K(NH4HSO4) < (NH4)(HSO4) always - this implies that there will
!  always be too much NH4 to satisfy the ammonium bisulphate formation
!  equation - while local minima (not roots) may exist, the chemically
!  correct course in this case is to assume that NH4 is the minimum
!  value.  In order to ensure that this is done, both boundaries of the
!  root search are set to the minimum value of NH4:
         do i = 1, ne1
            if (bflag(i) == 0 .and. diff(i) < 0.) then
               nh4_lo(i) = small
               nh4_hi(i) = small
               bflag(i) = 1
            end if
         end do

!  Check to make sure that some boundaries for the binary search were found:
         jsum = sum(bflag)
         if (jsum < ne1) then
            write(*, *) '> WARNING <'
            write(*, *) '> Root interval not set in Case 1'
            write(*, *) '> Condition occurred at the following gridpoint:'
            do i = 1, ne1
               if (bflag(i) == 0) write(*, *) '> ', i
            end do
         end if
      end if

!  If the initial subdivision of the variable space has been completed,
!  start the bisection search refinement of the chosen function interval
      if (rooteval >= ndiv + 2) rootsrch = rootsrch + 1

!  Reevaluate the function on the lhs of the interval:
      if (rootsrch == 1) then
         do i = 1, ne1
            nh4(i) = nh4_lo(i)
         end do
      end if

!  Reevaluate the function on the rhs of the interval:
      if (rootsrch == 2) then
         do i = 1, ne1
            nh4(i) = nh4_hi(i)
         end do
      end if

!  On the third function evaluation, the bisection search
!  has commenced.  The first stage of this is to orient the
!  search so that the positive value of the root equation
!  lies on the right of the interval.  Note that from this
!  point forward, the "_hi" value ceases to be the upper
!  limit of the variable - it is now the increment in the
!  variable as the function proceeds towards the root.
!
      if (rootsrch == 3) then
         do i = 1, ne1
!  "choice" is -1 if diff_lo<0, 1 otherwise:
            choice(i) = sign(1.E+0, real(diff_lo(i)))
            dx(i) = choice(i) * (nh4_lo(i) - nh4_hi(i))

!   "choice" is 1 if diff_lo<0, 0 otherwise
            choice(i) = 0.5D0 * (1.D0 - choice(i))
!
            rtbis(i) = choice(i) * nh4_lo(i) + (1.D0 - choice(i)) * nh4_hi(i)
            nh4_hi(i) = 0.5D0 * dx(i)
            nh4_lo(i) = rtbis(i)
            nh4(i) = nh4_lo(i) + nh4_hi(i)
!  rtbis is used to store the value of dx for future convergence calcs.
            rtbis(i) = max(abs(dx(i)), small)
         end do
      end if

!  All bisection evaluations:  determine which side of the
!  interval contains the root crossing, rearrange the bounds, and
!  subdivide the interval again:
      delouter = 0.D0
      if (rootsrch > 3) then
         do i = 1, ne1
!  "choice" is 1 if diff<0, 0 otherwise:
            choice(i) = 0.5D0 * (1.D0 - sign(1.E+0, real(diff(i))))
!
            nh4_lo(i) = nh4(i) * choice(i) + (1.D0 - choice(i)) * nh4_lo(i)
!
            delouter  = max(delouter, abs(nh4_hi(i) / rtbis(i)))
            nh4_hi(i) = nh4_hi(i) * 0.5D0
            nh4(i)    = nh4_lo(i) + nh4_hi(i)
         end do
      end if

!  The convergence parameter "del" and the number of bisection
!  iterations "iter" are used to exit the outer loop.

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

         call mach_hetv_water(ne1, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)

!  (2)  Calculate new concentrations of ions and other phases; multiply 
!       old concentration by ratio of old to new liquid water contents
         do i = 1, ne1
            lwo_lwn     = lwo(i) / lwn(i)
            ta(i)    = max(ta(i)    * lwo_lwn, lwmin)
            ts(i)    = max(ts(i)    * lwo_lwn, lwmin)
            tn(i)    = max(tn(i)    * lwo_lwn, lwmin)
            hno3(i)  = max(hno3(i)  * lwo_lwn, lwmin)
            so4(i)   = max(so4(i)   * lwo_lwn, lwmin)
            h(i)     = max(h(i)     * lwo_lwn, lwmin)
            no3(i)   = max(no3(i)   * lwo_lwn, lwmin)
            nh4(i)   = max(nh4(i)   * lwo_lwn, lwmin)
            hso4(i)  = max(hso4(i)  * lwo_lwn, lwmin)
            ambis(i) = max(ambis(i) * lwo_lwn, lwmin)
!
!  Update liquid water content
            lwo(i) = lwn(i)
         end do

!  Calculate activity coefficients
         call mach_hetv_activity(ne1, t, h, hso4, nh4, no3, so4, g_h_hso4,  &
                                 g_h2_so4, g_nh4_hso4=g_nh4_hso4)

!  Calculate maximum change in activity coefficients since the previous call
!  in this subroutine
         del = 0.D0
         do i = 1, ne1
            del = max(del, abs((g_h_hso4(i) - g_h_hso4_o(i)) / g_h_hso4(i)))
            g_h_hso4_o(i)   = g_h_hso4(i)
            del = max(del, abs((g_nh4_hso4(i) - g_nh4_hso4_o(i)) / g_nh4_hso4(i)))
            g_nh4_hso4_o(i) = g_nh4_hso4(i)
            del = max(del, abs((g_h2_so4(i) - g_h2_so4_o(i)) / g_h2_so4(i)))
            g_h2_so4_o(i)   = g_h2_so4(i)
         end do
!
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
         do i = 1, ne1
! 5:  HSO4- <=> SO4= + H+
!  Multiply by the square of the H-HSO4 activity coefficient,
!  and divide by the cube of the H2SO4 activity coefficient:
            khso4_m = khso4(i) * (g_h_hso4(i)**2) / (g_h2_so4(i) * g_h2_so4(i)**2)
!
!  Solve inner system of equations.
!  The value of NH4 is varied from 0 to TA in a two stage binary search. In the first stage,
!  system (1) is solved with the given value of NH4 as a constant.  The resulting HSO4 value is
!  tested using system (2):  the overall result is a search for the root of the ammonium
!  bisulphate equation as a function of the available NH4.  Note that TA/TS <1 if this
!  subroutine is called; the value within the root is always positive, since "c" is always negative
            b = khso4_m - ta(i) + ts(i)
            c = khso4_m * (ta(i) - nh4(i) - ts(i))
!
!  Taylor series used if round-off error possible:
            if (b /= 0.D0) then
               d = c / (b * b)
               v = 4.d0 * d
            else
               v = 1.D+03
            end if
            if (abs(v) <= smrt .and. b >  0.D0) then
               so4(i) =  - ((((14.D0 * d + 5.D0) * d + 2.D0) * &
                               d + 1.D0) * d + 1.D0) * c / b
            else
               so4(i) = 0.5d0 * (-b + sqrt(b * b - 4.d0 * c))
            end if
         end do
   
         do i = 1, ne1
            so4(i)   = max(so4(i), 0.0D0)
            ambis(i) = max(ta(i) - nh4(i), 0.D0)
            h(i)     = max(so4(i) + ts(i) - ta(i), 0.D0)
            hso4(i)  = max(ts(i) - so4(i) - ta(i) + nh4(i), 0.D0)
!      
! 6:  NH4HSO4(s) <=> NH4+ + HSO4-
!  Divide by the square of the ammonium bisulphate activity coefficient:
            kambis_m = kambis(i) / (g_nh4_hso4(i) * g_nh4_hso4(i))
            diff(i) = kambis_m - nh4(i) * hso4(i)
         end do

      end do ! End inner loop

!   Activity coefficients changed by less than eps2 during
!   the previous iteration; the first part of the problem, for
!   the lower limit of the variable (nh4), has converged.
!   Evaluate the root for the final equation to get the lower
!   bound value of the root.
      if (rooteval == 1 .or. rootsrch == 1) then
         do i = 1, ne1
            diff_lo(i) = diff(i)
         end do
      end if

!  Once the first midpoint has been calculated, use the
!  convergence checks from the last iteration to decide on
!  whether or not the set of problems have converged.  If
!  not, do another iteration.  If so, continue on and fill in
!  the remaining terms in the system of equations.
   end do
   
   do i = 1, ne1
      hno3(i)  = tn(i)
      no3(i)   = 0.D0
      nh3(i)   = 0.D0
!  Convert values back to moles/kg air from moles/kg water by
!  multiplying by the liquid water content:
      hno3(i)  = hno3(i)  * lwn(i)
      so4(i)   = so4(i)   * lwn(i)
      h(i)     = h(i)     * lwn(i)
      nh4(i)   = nh4(i)   * lwn(i)
      hso4(i)  = hso4(i)  * lwn(i)
      ambis(i) = ambis(i) * lwn(i)
   end do

   so4_i   = unpack(so4,   ncas == 1, so4_i)
   no3_i   = unpack(no3,   ncas == 1, no3_i)
   nh4_i   = unpack(nh4,   ncas == 1, nh4_i)
   hso4_i  = unpack(hso4,  ncas == 1, hso4_i)
   hno3_i  = unpack(hno3,  ncas == 1, hno3_i)
   h_i     = unpack(h,     ncas == 1, h_i)
   nh3_i   = unpack(nh3,   ncas == 1, nh3_i)
   ambis_i = unpack(ambis, ncas == 1, ambis_i)
   lwn_i   = unpack(lwn,   ncas == 1, lwn_i)

   return
end subroutine mach_hetv_case1
