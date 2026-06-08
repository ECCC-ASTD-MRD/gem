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
! Fichier/File   : mach_hetv_case12.ftn90
! Creation       : P. Makar, V. Bouchet, A. Nenes, S. Gravel, B. Pabla, S. Menard
! Description    : Heterogeneous chemistry solver for case12.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne12 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case 12's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne12.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                  TA/TS >= 2.0, rh > drh_amsul
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!                  (1)
!                     SO4 = TS
!                     TA = NH3 + NH4
!                     TN = HNO3 + NO3
!                     HNO3 <=> H + NO3, kHNO3
!                     NH3 + H2O <=> NH4 + kH2O (H2O/H), kNH3
!                  (2)
!                     H + NH4 = 2 SO4 + NO3 + kH2O (H2O / H)
!                  (3)
!                     NH3 + H2O <=> NH4 + kH2O (H2O)/(H), kNH3
!                  That is, stages 1 to 3 are for the system in which the only
!                  form of sulphur is SO4, and the concentration of OH is
!                  replaced by the corresponding water equilibrium expression.
!                  The value of HSO4 is determined AFTER a converged solution
!                  for (1) to (3) is found.  This is done using the HSO4
!                  equilibrium equation and mass conservation:
!                  (4)
!                     HSO4 <=> H + SO4, kHSO4
!
!                  The system of equations to be solved in this case are as follows:
!
!                  (1)
!                     SO4 = TS
!                     TA - NH3 - NH4 = 0
!                     TN - HNO3 -  NO3 = 0
!                     HNO3 = HNO3dry - delNO3,
!                     NO3 = NH4NO3dry +  delNO3,
!                     kHNO3 HNO3 = (NH4)(NO3)kH2O/(NH3 kNH3)
!                  (2)
!                  H2 + (NH4 - 2 SO4 - NO3) H - kH2O H2O = 0
!                  Note that activity coefficient convergence is over steps (1) and (2),
!                  and root convergence is over step (3).  Once the root has converged,
!                  step (4) is evaluated.
!                  (3)
!                  kNH3 (NH3) (H)/(kH2O NH4) - 1 = 0
!                  (4)
!                  kHSO4 (HSO4in + (HSO4) = (H-delHSO4)(SO4 - delHSO4)
!                  Hnew = H -delHSO4,
!                  SO4new = SO4 -delHSO4
!                  HSO4 = delHSO4
!
!                  The solution to the system of equations:
!
!                  (1)
!                     NH4 = TA - NH3
!                     SO4 = TS
!                     delNO3 = [kHNO3 kNH3 (NH3)(HNO3dry) -  kw(NH4) (NH4NO3dry)]
!                        /   [kw (NH4) + kHNO3 kNH3 NH3]
!                     HNO3 = HNO3dry - delNO3,
!                     NO3 = NH4NO3dry +  delNO3,
!                  (2)
!                     b = (NH4 - 2 SO4 - NO3),
!                     c = - kH2O H2O
!                     H = 0.5*(-b + sqrt(b**2 -4c)),
!                  (3)
!                     (NH3)(H)kNH3 - (NH4)kH2O = 0
!                  (4)
!                     delHSO4 = smallest positive real root of
!                     delHSO4^2 - (H + SO4 + kHSO4) delHSO4
!                        + [(H)(SO4)] = 0, or zero.
!                     Hnew = H -delHSO4,
!                     SO4new = SO4 -delHSO4
!                     HSO4 = delHSO4, if positive, otherwise zero.
!                  Variable, Range
!                     NH3, 0 to NH3dry
!                  ----------------------------------------------------------------
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case12(npts, nr, ne12, so4_i, no3_i, nh4_i, hso4_i, hno3_i,   &
                            h_i, nh3_i, lwn_i, t_i, rh_i, rho_i, ta_i, ts_i, tn_i, &
                            k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_case9, mach_hetv_water, mach_hetv_activity
   use mach_hetv_mod,         only: tstd, small, smrt, iter, itero, eps, eps2, rg
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne12
   integer(kind=4), intent   (in) :: ncas(npts)
   real(kind=8),    intent   (in) :: k0(nr)
   real(kind=8),    intent   (in) :: p1(nr)
   real(kind=8),    intent   (in) :: p2(nr)
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
   real(kind=8),    intent   (in) :: ta_i       (npts)
   real(kind=8),    intent   (in) :: ts_i       (npts)
   real(kind=8),    intent   (in) :: tn_i       (npts)
!!if_off
!
!  Local variables:
!
   real(kind=8), dimension(ne12) :: so4, no3, nh4, hso4, hno3, h, nh3, lwn
   real(kind=8), dimension(ne12) :: t, aw, rho
   real(kind=8), dimension(ne12) :: ta, ts, tn
!
   integer(kind=4)               :: i, jo, rooteval
   real(kind=8)                  :: delouter, del
!  Workspace arrays:
   real(kind=8)                  :: lwo_lwn, choice, dx, delno3
   real(kind=8)                  :: b, c, v, root, d
   real(kind=8), dimension(ne12) :: nh3dry, hno3dry, amnitdry, amsuldry
   real(kind=8), dimension(ne12) :: g_h_no3, g_h_no3_o
   real(kind=8), dimension(ne12) :: g_h_hso4, g_h_hso4_o
   real(kind=8), dimension(ne12) :: g_nh4_no3, g_nh4_no3_o
   real(kind=8), dimension(ne12) :: g_h2_so4, g_h2_so4_o
   real(kind=8), dimension(ne12) :: lwo, awu, law
   real(kind=8), dimension(ne12) :: nh3_lo, nh3_hi
   real(kind=8), dimension(ne12) :: rtbis
   real(kind=8), dimension(ne12) :: diff_lo, diff

!  arrays for reaction rates:
   real(kind=8), dimension(ne12) :: khso4, kh2o, knh3, khno3, kamnit
   real(kind=8) :: khso4_m, khno3_m, knh3_m, work1, work2, work3

!   so4  = pack(so4_i,  ncas == 12)
!   nh4  = pack(nh4_i,  ncas == 12)
!   no3  = pack(no3_i,  ncas == 12)
!   hso4 = pack(hso4_i, ncas == 12)
   hno3 = pack(hno3_i, ncas == 12)
   nh3  = pack(nh3_i,  ncas == 12)
!   h    = pack(h_i,    ncas == 12)
   ta   = pack(ta_i,   ncas == 12)
   ts   = pack(ts_i,   ncas == 12)
   tn   = pack(tn_i,   ncas == 12)
   aw   = pack(rh_i,   ncas == 12)
   t    = pack(t_i,    ncas == 12)
   rho  = pack(rho_i,  ncas == 12)

!  Initial guess for components:  SO4 = TS, HSO4 = 0,
!  NH4 = TA - NH3dry, and H = max(0, 2 SO4 + HSO4 - NH4), (NH4)2SO4 = 0
!
!  Initial guess for components:  SO4 = TS, HSO4 = 0,
!  NH4 = NH4NO3dry, NO3 = NH4NO3dry, H = 0
   do i = 1, ne12
      lwn(i) = 1.D0
!
!  Calculate rates for Equilibrium constants
      work1 = tstd / t(i) - 1.D0
      work2 = 1.d0 + log(tstd / t(i)) - tstd / t(i)
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
!  Calculate OH concentration for activity coefficient calculations
      kh2o(i) = k0(10) * exp(p1(10) * work1 + p2(10) * work2)
!  9:  NH4NO3(s) <=> NH3(g_eq) + HNO3(g_eq)
!  Original units are atm^2:  Gas is in molal units already; convert rate constant to (kg^3 air)/mole:
      kamnit(i) = (k0(9) * exp(p1(9) * work1 + p2(9) * work2)) / (work3 * work3) 
!
!  polynomials for electrolytes only good between aw=0.002 and aw=0.98
      awu(i) = min(max(aw(i), 2.0d-03), 0.98d0)
      law(i) = dble(log(real(awu(i))))    !double precision log unnec.
   end do

   call mach_hetv_case9(ne12, ne12, hno3dry, nh3dry, amsuldry, amnitdry, &
                        kamnit, ta, ts, tn)

!  Looping point for outer iteration (search):
   rooteval = 0
   delouter = 1.0d0
   do while (rooteval <= 3 .or. (delouter > eps .and. rooteval < iter))

!   Initial guess for liquid water content will be replaced in calculations to follow:
      do i = 1, ne12
         lwo(i) = 1.D0

!   Reset concentrations to starting liquid water content
         ta(i)    = ta(i)    * lwn(i)
         ts(i)    = ts(i)    * lwn(i)
         nh3(i)   = nh3(i)   * lwn(i)
         hno3(i)  = hno3(i)  * lwn(i)
         so4(i)   = 0.0D0    !so4(i)   * lwn(i)

         h(i)     = 0.0D0    !h(i)     * lwn(i)
         no3(i)   = amnitdry(i) * lwn(i)  !no3(i)   * lwn(i)
         nh4(i)   = amnitdry(i) * lwn(i)  !nh4(i)   * lwn(i)
         hso4(i)  = 0.0D0    !hso4(i)  * lwn(i)
         nh3dry(i)   = nh3dry(i)* lwn(i)
         hno3dry(i)  = hno3dry(i) * lwn(i)
         amnitdry(i) = amnitdry(i) * lwn(i)

!  Zero the 'old value' activity coefficient arrays
         g_h_hso4_o(i)  = 0.D0
         g_h2_so4_o(i)  = 0.D0
         g_h_no3_o(i)   = 0.D0
         g_nh4_no3_o(i) = 0.D0
      end do

      rooteval = rooteval + 1

!  First root evaluation is for lower bound, nh3 = 0
      if (rooteval == 1) then
         do i = 1, ne12
            nh3(i)    = 0.D0  !-(amnitdry(i)+amsuldry(i))  !0.D0
            nh3_lo(i) = nh3(i)
         end do
      end if

!  Second root evaluation is for upper bound, nh3=nh3dry
      if (rooteval == 2) then
         do i = 1, ne12
            nh3(i)    = nh3dry(i) + 2.D0 * amsuldry(i) + amnitdry(i) !nh3dry(i)
            nh3_hi(i) = nh3(i)
         end do
      end if

!  On the third root evaluation, the bisection search
!  has commenced.  The first stage of this is to orient the
!  search so that the positive value of the root equation
!  lies on the right of the interval:
      if (rooteval == 3) then
         do i = 1, ne12
!  "choice" is -1 if diff_lo<0, 1 otherwise:
            choice = sign(1.D0, diff_lo(i))
            dx = choice * (nh3_lo(i) - nh3_hi(i))
!   "choice" is 1 if diff_lo<0, 0 otherwise
            choice = 0.5D0 * (1.D0 - choice)
            rtbis(i) = choice * nh3_lo(i) + (1.D0 - choice) * nh3_hi(i)
            nh3_hi(i) = 0.5D0 * dx
            nh3_lo(i) = rtbis(i)
            nh3(i) = nh3_lo(i) + nh3_hi(i)
!  rtbis is used to store the value of dx for future convergence calcs.
            rtbis(i) = max(abs(dx), small)
         end do
      end if
!  All subsequent root evaluations:  determine which side of the
!  interval contains the root crossing, rearrange the bounds, and
!  subdivide the interval again:
      if (rooteval >  3) then
         delouter = 0.D0
         do i = 1, ne12
!  "choice" is 1 if diff<0, 0 otherwise:
            choice = 0.5D0 * (1.D0 - sign(1.D0, diff(i)))
            nh3_lo(i) = nh3(i) * choice + (1.D0 - choice) * nh3_lo(i)
            delouter  = max(delouter, abs(nh3_hi(i) / rtbis(i)))
            nh3_hi(i) = nh3_hi(i) * 0.5D0
            nh3(i)    = nh3_lo(i) + nh3_hi(i)
         end do
      end if

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
         call mach_hetv_water(ne12, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)

!  (2)  Calculate new concentrations of ions and other phases; multiply
!       old concentration by ratio of old to new liquid water contents
         do i = 1, ne12
            lwo_lwn     = lwo(i) / lwn(i)

            ta(i)       = ta(i)      * lwo_lwn
            ts(i)       = ts(i)      * lwo_lwn
            nh3(i)      = nh3(i)     * lwo_lwn
            hno3(i)     = hno3(i)    * lwo_lwn
            so4(i)      = so4(i)     * lwo_lwn

            h(i)        = h(i)       * lwo_lwn
            no3(i)      = no3(i)     * lwo_lwn
            nh4(i)      = nh4(i)     * lwo_lwn
            hso4(i)     = hso4(i)    * lwo_lwn

            nh3dry(i)   = nh3dry(i)   * lwo_lwn
            hno3dry(i)  = hno3dry(i)  * lwo_lwn
            amnitdry(i) = amnitdry(i) * lwo_lwn
!  Update liquid water content
            lwo(i) = lwn(i)
         end do
!
!  Calculate activity coefficients
         call mach_hetv_activity(ne12, t, h, hso4, nh4, no3, so4, g_h_hso4, &
                                 g_h2_so4, g_h_no3, g_nh4_no3)

!  Calculate maximum change in activity coefficients since the previous call
!  in this subroutine:

         del = 0.D0
         do i = 1, ne12
            del = max(del, abs((g_h_hso4(i) - g_h_hso4_o(i)) / g_h_hso4(i)))
            g_h_hso4_o(i) = g_h_hso4(i)
            del = max(del, abs((g_h2_so4(i) - g_h2_so4_o(i)) / g_h2_so4(i)))
            g_h2_so4_o(i) = g_h2_so4(i)
            del = max(del, abs((g_h_no3(i) - g_h_no3_o(i)) / g_h_no3(i)))
            g_h_no3_o(i) = g_h_no3(i)
            del = max(del, abs((g_nh4_no3(i) - g_nh4_no3_o(i)) / g_nh4_no3(i)))
            g_nh4_no3_o(i) = g_nh4_no3(i)
         end do
!
         do i = 1, ne12
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

! 1, 2:  HNO3(g_eq) <=> HNO3(aq) <=> H+ + NO3-
! k1: Original units (moles HNO3(aq)/kg H2O)/(atm HNO3(g)), but
! gas is already in moles/(kg H2O); multiply by (R T Lwn rho).
! k2: Divide by the square of the HNO3 activity coefficient
            khno3_m = khno3(i) * lwn(i) / (g_h_no3(i)**2)
!
! 3:  NH3(g_eq) <=> NH3(aq), NH3(aq) + H2O <=> NH4+ + OH-
!  ==> NH3(g) + H2O <=> NH4+ + OH-
!
! k3: Original units (moles NH3(aq)kg H2O)/(atm NH3(g)), but
! gas is already in moles/(kg H2O); multiply by (R T Lwn rho).
! k4: Multiply by the square of the HNO3 activity coefficient,
!  and divide by (the square of the NH4NO3 activity coefficient
!  times Kh2o (=Kw).
!  i.e. k4=(NH4)(OH)/( (NH3)(H2O) ), kw = (H)(OH)/(H2O)
!  Therefore,
!      k4/kw = (NH4)/( (NH3) (H) ) .  Then multiply by
!  gamma(H)/gamma(NO3) to correct for non-equilibrium.
!  Rather than divide by kH2O here, the kH2O value
!  is carried to the denominator in the relevant equations
            knh3_m = knh3(i) * lwn(i) * g_h_no3(i)**2 / (g_nh4_no3(i)**2)
!
!  Solve inner system of equations.
!  The value of NH3 is varied from 0 to NH3dry in a two stage
!  binary search.  In the first stage, system (1&2) is solved with the
!  given value of NH3 as a constant.  The resulting H and NH4 values
!  are tested using system (3):  the overall result is a search for the
!  root of the NH4 <=> NH3 equilibrium equation as a function of the
!  available NH3.  Note that c is defined negative, therefore
!  the root is always real and larger in magnitude than b.
!  If b is negative, the positive root results in a positive H;
!  if b is positive, the positive root results in a positive H.
!  The negative root will result in a negative value of H, regardless
!  of the sign of b.  Positive root is therefore used.

!  Note that kh2o = kw in the following, and that the aw values cancel out.
            so4(i) = ts(i)
            nh4(i) = max(0.D0, ta(i) - nh3(i))
!            nh3(i) = max(0.D0, nh3dry(i)-nh3(i))
!  Correction, Nov17, 2008:  if the NH4+ = 0, then
!  the ammonia equilibrium can no longer be used in the system
!  of equations, and the delno3 equation is invalid. P. Makar/J. Murphy
            if (nh4(i) /= 0.D0) then
               delno3 = min(max(( khno3_m * knh3_m * nh3(i) * hno3dry(i) &
                            - kh2o(i) * nh4(i) * amnitdry(i)) / &
                        (kh2o(i)*nh4(i) + khno3_m * knh3_m * nh3(i)), 0.D0), hno3dry(i))
               hno3(i) = hno3dry(i) - delno3
               no3(i) = amnitdry(i) + delno3
            else
!  NH4+ = 0; set HNO3 = HNO3dry, and NO3- = 0.
               hno3(i) = hno3dry(i)
               no3(i) = 0.D0
            end if

!  Note that a taylor expansion of the square root is used in the
!  instance that round-off error can mess up the root for the
!  hydrogen equation.
!  Determine the concentration of H required to maintain
!  charge balance.  Note that c is defined negative, therefore
!  the root is always real and larger in magnitude than b.
!  If b is negative, the positive root results in a positive H;
!  if b is positive, the positive root results in a positive H.
!  The negative root will result in a negative value of H, regardless
!  of the sign of b.  Positive root for the quadratic is therefore used.
!  Note that a taylor expansion of the square root is used in the
!  instance that round-off error can mess up the root.
!  "kh2o" = kw in the following.  Note that gamma H+ x gamma OH- = 1 is assumed.
            b = nh4(i) - 2.D0 * so4(i) - no3(i)
            c =  - kh2o(i) * aw(i)
            if (b /= 0.D0) then
               d = c / (b * b)
               v = 4.d0 * d
            else
               v = 1.D03
            end if
            if (abs(v) <= smrt .and. b > 0.D0) then
               h(i) = - ((((14.D0 * d + 5.D0) * d + 2.D0) * &
                            d + 1.D0) * d + 1.D0) * c / b
            else
               h(i) = 0.5d0 * (-b + sqrt(b * b - 4.d0 * c) )
            end if
            h(i) = max(h(i), 0.D0)
!
!   Evaluate the root for the final equation to get the lower
!   bound value of the root.
            diff(i) = knh3_m * nh3(i) * h(i) - nh4(i) * kh2o(i)
         end do
         
      end do ! End inner loop
!
!   Activity coefficients changed by less than eps2 during
!   the previous iteration; the first part of the problem, for
!   the lower limit of the variable (nh4), has converged.
      if (rooteval == 1) then
         do i = 1, ne12
            diff_lo(i) = diff(i)
         end do
      end if

!  Repeat the calculations, for the upper boundary of the interval
!  in which the root is expected, and the first midpoint evaluation:
   end do

   do i = 1, ne12
! 5:  HSO4- <=> SO4= + H+
!  Multiply by the square of the H-HSO4 activity coefficient,
!  and divide by the cube of the H2SO4 activity coefficient:
      khso4_m = khso4(i) * (g_h_hso4(i)**2) / (g_h2_so4(i) * g_h2_so4(i)**2)

!  Final stage is an adjustment of the so4 values in order to
!  allow for the possible existence of hso4.  Note the b is
!  defined negative, so -b is positive, and the (non-complex)
!  quantity within the root will be a number smaller in magnitude
!  than b.  The sign of the root must therefore be positive
!  in order to get the smaller of the two roots.
      b = - (h(i) + so4(i) + khso4_m)
      c =    h(i) * so4(i)
      root = b * b - 4.D0 * c
!   Modification:  need to add condition that if there's no SO4
!  to convert into HSO4, you don't bother to take the root:
      if (root >= 0.D0 .and. so4(i) > 0.D0) then
         hso4(i) = 0.5D0 * ( - b - sqrt(root))
      else
         hso4(i) = 0.D0
      end if
      h(i) = max(h(i) - hso4(i), 0.D0)
      so4(i) = max(so4(i) - hso4(i), 0.D0)
   end do

!  Convert values back to moles/kg air from moles/kg water by
!  multiplying by the liquid water content:
   do i = 1, ne12
      nh3(i)   = nh3(i)   * lwn(i)
      hno3(i)  = hno3(i)  * lwn(i)
      so4(i)   = so4(i)   * lwn(i)
      h(i)     = h(i)     * lwn(i)
      no3(i)   = no3(i)   * lwn(i)
      nh4(i)   = nh4(i)   * lwn(i)
      hso4(i)  = hso4(i)  * lwn(i)
   end do

   so4_i  = unpack(so4,  ncas == 12, so4_i)
   nh4_i  = unpack(nh4,  ncas == 12, nh4_i)
   no3_i  = unpack(no3,  ncas == 12, no3_i)
   hso4_i = unpack(hso4, ncas == 12, hso4_i)
   hno3_i = unpack(hno3, ncas == 12, hno3_i)
   nh3_i  = unpack(nh3,  ncas == 12, nh3_i)
   h_i    = unpack(h,    ncas == 12, h_i)
   lwn_i  = unpack(lwn,  ncas == 12, lwn_i)

   return
end subroutine mach_hetv_case12

