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
! Fichier/File   : mach_hetv_poly3v.ftn90
! Creation       : V. Bouchet, P. Makar, S. Menard, S. Gravel, B. Pabla
! Description    :
! Extra info     :
!
! Arguments  IN
!
!            OUT
!
!=============================================================================
!
!!if_on
subroutine mach_hetv_poly3v (a1, a2, a3, root, islv, ne)
!!if_off
   use chm_consphychm_mod, only: pi
   use chm_utils_mod,      only: chm_error_l
   implicit none
!!if_on
   integer(kind=4), intent (in) :: ne
   real(kind=8),    intent (in) :: a1  (ne)
   real(kind=8),    intent (in) :: a2  (ne)
   real(kind=8),    intent (in) :: a3  (ne)
   real(kind=8),    intent(out) :: root(ne)
   integer(kind=4), intent(out) :: islv(ne)
!!if_off
!
! Local variables
!
   logical(kind=4), parameter :: bisect=.false.
   integer(kind=4), parameter :: jmax = 40
   integer(kind=4)            :: i, j
   integer(kind=4)            :: ix(ne)
   real(kind=8)               :: thet, coef, s, sqd, ssig, tsig, t, fu
   real(kind=8)               :: x(ne, 3)
   real(kind=8)               :: d(ne), q(ne), r(ne)
   real(kind=8), parameter    :: expon = 1.0d0 / 3.0d0, zero = 0.0d0,                  &
                                 thet1 = 120.0d0 / 180.0d0, thet2 = 240.0d0 / 180.0d0, &
                                 eps = 1.d-50

   do i = 1, ne
      islv(i) = 1
      root(i) = 0.0d0
   end do

   do j = 1, 3
      do i = 1, ne
         x(i, j) = 0.0d0
      end do
   end do

!  special case : quadratic*x equation
   do i = 1, ne
      d(i) = a1(i) * a1(i) - 4.0d0 * a2(i)
   end do

   do i = 1, ne
      if (abs(a3(i)) <= eps) then
         ix(i) = 1
         x(i, 1) = zero
      end if
   end do

   do i = 1, ne
      if (abs(a3(i)) <= eps .and. d(i) >= zero) then
         ix(i) = 3
         x(i, 2) = 0.5d0 * (-a1(i) + sqrt(d(i)))
         x(i, 3) = 0.5d0 * (-a1(i) - sqrt(d(i)))
      end if
   end do

!  normal case : cubic equation
!  define parameters q, r, s, t, d
   do i = 1, ne
      q(i) = (3.0d0 * a2(i) - a1(i) * a1(i)) / 9.0d0
      r(i) = (9.0d0 * a1(i) * a2(i) - 27.0d0 * a3(i) &
                           - 2.0d0 * a1(i) * a1(i) * a1(i)) / 54.0d0
      d(i) = q(i) * q(i) * q(i) + r(i) * r(i)
   end do

!  calculate roots *
!  d < 0, three real roots
   do i = 1, ne
      if (abs(a3(i)) > eps .and. d(i) < -eps) then
         ix(i) = 3
         thet = expon * acos(r(i) / sqrt(-q(i) * q(i) * q(i)))
         coef = 2.0d0 * sqrt(-q(i))
         x(i, 1) = coef * cos(thet) - expon * a1(i)
         x(i, 2) = coef * cos(thet + thet1 * pi) - expon * a1(i)
         x(i, 3) = coef * cos(thet + thet2 * pi) - expon * a1(i)
      end if
   end do

!  d = 0, three real (one double) roots
   do i = 1, ne
      if (abs(a3(i)) > eps .and. d(i) >= -eps .and. d(i) <= eps) then
         ix(i) = 2
         s = (abs(r(i)) ** expon) * sign(1.0e0, real(r(i)))
         x(i, 1) = 2.0d0 * s - expon * a1(i)
         x(i, 2) = -s - expon * a1(i)
      end if
   end do

!  d > 0, one real root
   do i = 1, ne
      if (abs(a3(i)) > eps .and. d(i) > eps) then
         ix(i) = 1
         sqd = sqrt(d(i))
         ssig = sign(1.0e0, real(r(i) + sqd))
         tsig = sign(1.0e0, real(r(i) - sqd))
         s = ssig * (abs(r(i) + sqd)) ** expon
         t = tsig * (abs(r(i) - sqd)) ** expon
         x(i, 1) = s + t - expon * a1(i)
      end if
   end do

!  select appropriate root
!  the following section of code sets islv to 1 if no positive roots
!  with magnitudes greater than eps are found:
      do i = 1, ne
         islv(i) = 1
      end do
      do j = 1, 3
         do i = 1, ne
            islv(i) =  min(islv(i), int(0.5d0 * (1.0d0 + sign(1.0d0, (eps - x(i, j))))))
         end do
      end do

!  next, set any calculated roots that are less than or equal to eps
!  to 1e30, then search for the smallest root of the system of equations
! (all expected values of the positive roots should be << 1e30):
   do i = 1, ne
      root(i) = 1.0d30
   end do
   do j = 1, 3
      do i = 1, ne
         x(i, j) = max(x(i, j), 1d30 * (0.5d0 * (1.0d0 - sign(1.0d0, x(i, j) - eps))))
      end do
   end do
   do j = 1, 3
      do i = 1, ne
         root(i) = min(root(i), x(i, j))
      end do
   end do

!  The given root may not be very accurate, due to the potentially
!  large dynamic range of the constants a1, a2, a3 (e.g. 10 orders of
!  magnitude between each).  A bisection search is employed to refine
!  the value of the root.  Algorithm for the boundaries of the search:
!    (a)  If the root was the largest one of those found, the upper
!  limit of the search is set to 10x that value.
!    (b)  If a larger root was available, the upper limit of the search
!  was set to the midpoint between the chosen root and that larger
!  root.
!    (c) The lower limit for the search is always set to zero, since
!  the desire is to find the smallest positive real root.
!  Note: q is now the upper x limit, r is the lower x limit,
!  d is the value of the cubic polynomial, x(i, 1) is the value of
!  dx.
!    The bisection search algorithm itself is from Numerical Recipes,
!  page 347, subroutine "rtbis"
   if (bisect) then
      do i = 1, ne
         do j = 1, ix(i)
            if (x(i, j) > root(i)) then
               q(i) = (root(i) + x(i, j)) * 0.5d0
            else
               q(i) = root(i) * 10.0d0
            end if
         end do
         r(i) = 0.0d0
      end do

!  Commence bisection search for the polynomial x^3 + a1 x^2 + a2 x + a3 = 0:
      do i = 1, ne
         fu = q(i) ** 3 + a1(i) * q(i) ** 2 + a2(i) * q(i) + a3(i)
         d(i) = r(i) ** 3 + a1(i) * r(i) ** 2 + a2(i) * r(i) + a3(i)
         if (fu * d(i) >= 0.0d0) then
            write(0, *) '### Error in mach_hetv_poly3v ###'
            write(0, *) '# root not bracketed in poly3v'
            write(0, *) '# gridpoint number: ', i
            write(0, *) '###         ABORT         ###'
            chm_error_l = .true.
            return
         end if
      end do

!  Orient direction of search:
      do i = 1, ne
         if (d(i) < 0.0d0) then
            root(i) = r(i)
            x(i, 1) = q(i) - r(i)
         else
            root(i) = q(i)
            x(i, 1) = r(i) - q(i)
         end if
      end do

!   Bisect jmax times:
      do j = 1, jmax
         do i = 1, ne
            x(i, 1) = x(i, 1) * 0.5d0
            q(i) = root(i) + x(i, 1)  !q now the midpoint value
            d(i) = q(i) ** 3 + a1(i) * q(i) ** 2 + a2(i) * q(i) + a3(i)
            if (d(i) <= 0.0d0) root(i) = q(i)
         end do
      end do
   end if

return
end
