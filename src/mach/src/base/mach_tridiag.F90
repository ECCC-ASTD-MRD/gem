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
! Fichier/File   : mach_tridiag.ftn90
! Creation       : Paul Makar, nov 2007
! Description    : Solution of the tridiagonal system of equations
!
! Extra info     : Based on Janusz Pudykiewicz, March 1995
!                  Follows algorithm found in "Numerical Recipes".
!                  Note that the potential exists for a division by zero error in the
!                  following, in that if array "bet" happens to = zero, then a division
!                  by zero error could occur.  However, the same algorithm has been
!                  used in AURAMS and CHRONOS for many years without encountering
!                  this error condition, and to add a check here would slow down the
!                  processing.
! Arguments:
!           IN
!              a, b, c -> matrix terms
!              r       -> field to be diffused   
!
!           OUT
!              u       -> updated diffused field 
!
!=============================================================================
!
!!if_on
subroutine mach_tridiag (a, b, c, r, u, dni, dnk)
!!if_off
   implicit none
!!if_on
   integer(kind=4), intent (in) :: dni, dnk
   real(kind=4),    intent (in) :: a(dni, dnk)
   real(kind=4),    intent (in) :: b(dni, dnk)
   real(kind=4),    intent (in) :: c(dni, dnk)
   real(kind=4),    intent (in) :: r(dni, dnk)
   real(kind=4),    intent(out) :: u(dni, dnk)
!!if_off
!
! Local Variables
!
   integer(kind=4) :: i, j
   real(kind=4)    :: bet(dni)
   real(kind=4)    :: gam(dni, dnk)
!
!--------------------------------------------------------------------------
!
!  solves for a vector u of length n the tridiagonal linear set:
!
!       b1  c1                        u1      r1
!       a2  b2  c2                 =  u2      r2
!
!
!                           aN  bN    uN      rN
!
!  a, b, c and r are input vectors and are not modified
!
!--------------------------------------------------------------------------
!
   do i = 1, dni
      bet(i) = b(i, 1)
      u(i, 1) = r(i, 1) /  bet(i)
   end do

      do j = 2, dnk
         do i = 1, dni
            gam(i, j) = c(i, j - 1) / bet(i)
            bet(i)    = b(i, j) - a(i, j) * gam(i, j)
            u  (i, j) = (r(i, j) - a(i, j) * u(i, j - 1)) / bet(i)
         end do
      end do
!
!  Backsubstitution
!
   do j = dnk - 1, 1, -1
      do i = 1, dni
         u(i, j) = u(i, j) - gam(i, j + 1) * u(i, j + 1)
      end do
   end do
   return
end subroutine mach_tridiag
