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
! Fichier/File   : mach_hetv_case9.ftn90
! Creation       : P. Makar,  V. Bouchet,  A. Nenes,  S. Gravel,  B. Pabla,  S. Menard
! Description    : Heterogeneous chemistry solver for case9.  Based on algorithm's
!                  in Nenes et al's ISORROPIA solver, recoded to vectorize over the
!                  gridpoint dimension.  All input arrays are of length npts (total
!                  number of gridpoints submitted from the calling code.  The
!                  subsection of each 1-D array between index 1 and ne9 (inclusive)
!                  has been pre-sorted, and contains the gridpoint data that must
!                  be solved using case 9's algorithm.  Operations within this
!                  subroutine therefore take place over array bounds from 1 to ne9.
!
!                  Units on input are moles/kg air!
!
! Extra info     : Athanasios Nenes, Mail Code 210-41, Dept. of Chemical Engineering,
!                  California Institute of Technology, Pasadena, Ca., 91125, USA
!                  nenes@its.caltech.edu
!
!                  The conditions for which this case is called are as follows:
!                        tats >= 2.0, rh <  mdrh_amnit_amsul
!                  The reactions representing this case are as follows(bracketed numbers
!                  indicate that the solution is done in stages, as numbered):
!                     NH4NO3 (s) <=> NH3 + HNO3, kNH4NO3
!                     TA = 2 (NH4)2SO4 + NH4NO3 + NH3
!                     TS = (NH4)2SO4
!                     TN = HNO3 + NH4NO3
!                  The system of equations to be solved in this case are as follows:
!                     kNH4NO3 - (NH3)(HNO3) = 0,
!                     TA - 2 (NH4)2SO4 - NH4NO3 - NH3 = 0,
!                     TN - HNO3 = 0,
!                     TS = (NH4)2SO4
!                  The solution to the system of equations:
!                     b = - TA + 2 (NH4)2SO4 + TN,
!                     c = - kNH4NO3,
!                     NH3 = 0.5*(-b + sqrt(b**2 -4c)),
!                     HNO3 = - TA + 2 (NH4)2SO4 + TN + NH3
!                     NH4NO3 = TA - 2 (NH4)2SO4 - NH3
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_case9(npts, ne9, hno3, nh3, amsul, amnit, &
                           kamnit, ta, ts, tn)
!!if_off
   use mach_hetv_mod,    only: smrt
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: npts
   integer(kind=4), intent   (in) :: ne9
   real(kind=8),    intent(inout) :: hno3  (npts)
   real(kind=8),    intent(inout) :: nh3   (npts)
   real(kind=8),    intent(inout) :: amsul (npts)
   real(kind=8),    intent(inout) :: amnit (npts)
   real(kind=8),    intent   (in) :: ta    (npts)
   real(kind=8),    intent   (in) :: ts    (npts)
   real(kind=8),    intent   (in) :: tn    (npts)
   real(kind=8),    intent   (in) :: kamnit(npts)
!!if_off
!
!  Local variables:
!
   integer(kind=4)              :: i
!  arrays for reaction rates:
   real(kind=8)                 :: b, c, d, v

!  Dry particle, so no water calculation and no iteration required
   do i = 1, ne9
!   Solve system of equations.
      amsul(i) = ts(i)
      b = -ta(i) + 2.0d0 * amsul(i) + tn(i)
      c = -kamnit(i)
      if (b /= 0.0d0) then
         d = c / (b * b)
         v = 4.0d0 * d
      else
         v = 1.0d+03
      end if
      if (abs(v) <= smrt .and. b > 0.0d0) then
         nh3(i) = -((((14.0d0 * d + 5.0d0) * d + 2.0d0) * &
                       d + 1.0d0) * d + 1.0d0) * c / b
      else
         nh3(i) = 0.5d0 * (-b + sqrt(b * b - 4.0d0 * c))
      end if
   end do

!  Note:  system of equations ensures mass conservation, but not
!  monotonicity.  Extra boundary conditions are therefore required on
!  the value of NH3.  Gas phase NH3 can be no larger than the value which
!  would zero ammonium nitrate (otherwise, the ammonium nitrate concentration
!  would be negative).  Values of less than zero for NH3 are possible if
!  round-off error affects the calculation; this is also corrected.
   do i = 1, ne9
      nh3(i)   = min(nh3(i), ta(i) - 2.0d0 * amsul(i))
      nh3(i)   = max(nh3(i), 0.0d0)
      amnit(i) = max(0.0d0, (ta(i) - 2.0d0 * amsul(i) - nh3(i)))
!  Correction:  P.A. Makar, Nov 2008
      if (amnit(i) > 0.0d0) then
         hno3(i) = max(-ta(i) + 2.0d0 * amsul(i) + tn(i) + nh3(i), 0.0d0)
      else
         hno3(i) = tn(i)
      end if
   end do

   return
end subroutine mach_hetv_case9
