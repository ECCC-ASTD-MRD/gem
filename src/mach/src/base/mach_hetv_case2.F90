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
! Fichier/File   : mach_hetv_case2.ftn90
! Creation       : P. Makar, V. Bouchet, A. Nenes, S. Gravel, B. Pabla, S. Menard
! Description    : Heterogeneous chemistry solver for case2.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne2 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case 2's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne2.
!
!                  Units on input are moles/kg air.
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                  tats <  1.0, rh >= drh_ambis
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!
!                     HSO4 <=> H + SO4 , kHSO4
!                     NH4 = TA
!                     H + NH4 = 2 SO4 + HSO4,
!                     TS = HSO4 + SO4
!                  The system of equations to be solved in this case are as follows:
!
!                     kHSO4 (HSO4) - (H)(SO4) = 0,
!                     NH4 = TA
!                     H + NH4 - 2 SO4 - HSO4,
!                     TS - HSO4 - SO4 = 0
!                  The solution to the system of equations:
!
!                     Let
!                     b = (kHSO4 - TA + TS),
!                     c  = - kHSO4 TS.
!                     SO4 = 0.5*(-b + sqrt(b**2 -4c)),
!                     H = - TA + SO4 + TS
!                     HSO4 = TS - SO4
!                     NH4 = TA
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case2(npts, nr, ne2, so4_i, no3_i, nh4_i, hso4_i, hno3_i, &
                           h_i, nh3_i, lwn_i, t_i, rh_i, ts_i, ta_i, tn_i,     &
                           k0, p1, p2, ncas)
!!if_off
   use mach_hetv_headers_mod, only: mach_hetv_water, mach_hetv_activity
   use mach_hetv_mod,         only: tstd, smrt, itero, eps2, small
   use chm_utils_mod,         only: chm_error_l
   implicit none   
!!if_on
   integer(kind=4), intent   (in) :: nr
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne2
   integer(kind=4), intent   (in) :: ncas (npts)
   real(kind=8),    intent   (in) :: k0   (nr)
   real(kind=8),    intent   (in) :: p1   (nr)
   real(kind=8),    intent   (in) :: p2   (nr)
   real(kind=8),    intent(inout) :: so4_i   (npts)
   real(kind=8),    intent(inout) :: no3_i   (npts)
   real(kind=8),    intent(inout) :: nh4_i   (npts)
   real(kind=8),    intent(inout) :: hso4_i  (npts)
   real(kind=8),    intent(inout) :: hno3_i  (npts)
   real(kind=8),    intent(inout) :: h_i     (npts)
   real(kind=8),    intent(inout) :: nh3_i   (npts)
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
   real(kind=8), dimension(ne2) :: so4, no3, nh4, hso4, hno3, h, nh3, lwn
   real(kind=8), dimension(ne2) :: t, aw
   real(kind=8), dimension(ne2) :: ta, ts, tn
!
   integer(kind=4)               :: i, jo
   real(kind=8)                  :: del, b, c, v, khso4_m, d, lwo_lwn
!  Workspace arrays:
   real(kind=8), dimension(ne2) :: g_h_hso4, g_h_hso4_o
   real(kind=8), dimension(ne2) :: g_h2_so4, g_h2_so4_o
   real(kind=8), dimension(ne2) :: lwo, awu, law
!  arrays for reaction rates:
   real(kind=8), dimension(ne2) :: khso4

!   so4   = pack(so4_i,  ncas == 2)
!   nh4   = pack(nh4_i,  ncas == 2)
   no3   = pack(no3_i,  ncas == 2)
!   hso4  = pack(hso4_i, ncas == 2)
   hno3  = pack(hno3_i, ncas == 2)
!   h     = pack(h_i,    ncas == 2)
   ta    = pack(ta_i,   ncas == 2)
   ts    = pack(ts_i,   ncas == 2)
   tn    = pack(tn_i,   ncas == 2)
   aw    = pack(rh_i,   ncas == 2)
   t     = pack(t_i,    ncas == 2)

   do i = 1, ne2
!  Initial guess for ions:  TS split evenly between SO4 and HSO4,
!  NH4 = TA, and H = max(0, 2 SO4 + HSO4 - NH4):
      so4(i)  = 0.5d0 * ts(i)
      hso4(i) = 0.5d0 * ts(i)
      nh4(i)  = ta(i)
      h(i)    = max(0.D0, 2.D0 * so4(i) + hso4(i) - nh4(i))

!   Initial guess for liquid water content (will be replaced in calculations to follow:
      lwo(i) = 1.D0

!  Zero the 'old value' activity coefficient arrays
      g_h_hso4_o(i)     = 0.D0
      g_h2_so4_o(i)     = 0.D0
!
!  Calculate rates for Equilibrium constants
!  5:  HSO4- <=> SO4= + H+
      khso4(i) = k0(5) * exp( p1(5) * (tstd / t(i) - 1.D+00)             &
               + p2(5) * (1.d+00 + log( tstd / t(i) ) - tstd / t(i) )  )
!
!  polynomials for electrolytes only good between aw=0.002 and aw=0.98
      awu(i) = min(max(aw(i), 2.0d-03), 0.98d0)
      law(i) = dble(log(real(awu(i))))    !double precision log unnec.
   end do

!  Looping point for outer iteration:
   jo = 0

   del = 1.0d0
!  Convergence check:  Have the activity coefficients
!  changed by more than eps2 in the last iteration?  If
!  so, repeat the calculation with the new activity
!  coefficients.  Keep track of the number of iterations.
   do while (del > eps2 .and. jo < itero)
      jo = jo + 1

!  Calculate liquid water, based on current ion concentration.
      call mach_hetv_water(ne2, so4, h, no3, nh4, hso4, awu, law, lwo, lwn)
!
!
!  (2)  Calculate new concentrations of ions and other phases; multiply old concentration by
!  ratio of old to new liquid water contents
      do i = 1, ne2
         lwo_lwn = lwo(i) / lwn(i)
         ta(i)   = ta(i)   * lwo_lwn
         ts(i)   = ts(i)   * lwo_lwn
         tn(i)   = tn(i)   * lwo_lwn
         hno3(i) = hno3(i) * lwo_lwn
         so4(i)  = so4(i)  * lwo_lwn
         h(i)    = h(i)    * lwo_lwn
         no3(i)  = no3(i)  * lwo_lwn
         nh4(i)  = nh4(i)  * lwo_lwn
         hso4(i) = hso4(i) * lwo_lwn
!
!  Update liquid water content
         lwo(i) = lwn(i)
      end do

!  Calculate activity coefficients
      call mach_hetv_activity(ne2, t, h, hso4, nh4, no3, so4, g_h_hso4, &
                              g_h2_so4)

!  Calculate maximum change in activity coefficients since the previous call in this subroutine:
      del = 0.D0
      do i = 1, ne2
         del = max(del, abs((g_h_hso4(i) - g_h_hso4_o(i)) / g_h_hso4(i)))
         g_h_hso4_o(i) = g_h_hso4(i)
         del = max(del, abs((g_h2_so4(i) - g_h2_so4_o(i)) / g_h2_so4(i)))
         g_h2_so4_o(i) = g_h2_so4(i)
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

! 5:  HSO4- <=> SO4= + H+
!  Multiply by the square of the H-HSO4 activity coefficient,
!  and divide by the cube of the H2SO4 activity coefficient:
         khso4_m = khso4(i) * g_h_hso4(i)**2 / (g_h2_so4(i) * g_h2_so4(i)**2)

!   Solve system of equations.
         b = khso4_m - ta(i) + ts(i)
         c = - khso4_m * ts(i)
         if (b /= 0.D0) then
            d = c / (b * b)
            v = 4.d0 * d
         else
            v = 1.D+03
         end if
         if (abs(v) <= smrt .and. b > 0.D0) then
            so4(i) = - ((((14.D0 * d + 5.D0) * d + 2.D0) * &
                           d + 1.D0) * d + 1.D0) * c / b
         else
             so4(i) = 0.5d0 * (-b + sqrt(b * b - 4.d0 * c))
         end if
         so4(i) = max(so4(i), 0.D0)
         h(i) = max((- ta(i) + so4(i) + ts(i)), 0.D0)
         hso4(i) = max((ts(i) - so4(i)), 0.D0)
         nh4(i) = ta(i)
      end do
      
   end do
!
!   Activity coefficients changed by less than eps2 during
!   the previous iteration; the problem has converged.  Zero
!   the particle nitrate, the gas-phase ammonia, and the salts,
!   and return to the calling code.

   do i = 1, ne2
      no3(i)   = 0.D0
      nh3(i)   = 0.D0
      hno3(i)  = tn(i)

!  Convert values back to moles/kg air from moles/kg water by
!  multiplying by the liquid water content:
      hno3(i)  = hno3(i)  * lwn(i)
      so4(i)   = so4(i)   * lwn(i)
      h(i)     = h(i)     * lwn(i)
      nh4(i)   = nh4(i)   * lwn(i)
      hso4(i)  = hso4(i)  * lwn(i)
   end do

!  CHECK
   if (any(h < 0.0D0)) then
      write(0, *) '### Error in mach_hetv_case2 ###'
      i = minloc(h, 1)
      write(0, *) 'Negative H+ value detected in case2 for i = ',i
      write(0, *) 'H+ : '   ,h(i)
      write(0, *) 'HSO4-: ' ,hso4(i)
      write(0, *) 'SO4=: '  ,so4(i)
      write(0, *) 'NO3-: '  ,no3(i)
      write(0, *) 'NH4+: '  ,nh4(i)
      write(0, *) 'HNO3: '  ,hno3(i)
      write(0, *) '###         ABORT         ###'
      chm_error_l = .true.
      return
   end if

   so4_i  = unpack(so4,  ncas == 2, so4_i)
   nh4_i  = unpack(nh4,  ncas == 2, nh4_i)
   no3_i  = unpack(no3,  ncas == 2, no3_i)
   hso4_i = unpack(hso4, ncas == 2, hso4_i)
   hno3_i = unpack(hno3, ncas == 2, hno3_i)
   h_i    = unpack(h,    ncas == 2, h_i)
   nh3_i  = unpack(nh3,  ncas == 2, nh3_i)
   lwn_i  = unpack(lwn,  ncas == 2, lwn_i)

   return
end subroutine mach_hetv_case2
