!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r wil_gauleg - Computation of the abscisses and weights of Gauss-Legendre

      subroutine wil_gauleg (x1,x2,x,w,n)

      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: n
      real, intent(in) :: x1,x2
      real, intent(out) :: x(n),w(n)

      !object
      !==================================================================================================
      !     Computation of the abscisses and weights of Gauss-Legendre
      !     Based on GAULEG: Extracted from Numerical Recipes in FORTRAN; The Art of Scientific Computing
      !                      Second Edition, Cambridge University Press New York, NY, USA, 1993
      !==================================================================================================

      !---------------------------------------------------------------

      real(kind=REAL64), parameter :: EPS_8 = 3.d-14
      integer :: i,j,m
      real(kind=REAL64) :: p1_8, p2_8, p3_8, pp_8, xl_8, xm_8, z_8, z1_8

      !---------------------------------------------------------------

      m = (n+1)/2

      xm_8 = 0.5d0*(x2+x1)
      xl_8 = 0.5d0*(x2-x1)

      do i=1,m

         z_8  = cos(pi_8 * (i-0.25d0) / (n+0.5d0))
         z1_8 = z_8 + 1.d14 * EPS_8 ! Anything bigger than EPS_8 will work

         do while (abs(z_8-z1_8) > EPS_8)

            p1_8 = 1.d0
            p2_8 = 0.d0

            do j=1,n

               p3_8 = p2_8
               p2_8 = p1_8
               p1_8 = ((2.d0*j-1.d0)*z_8*p2_8-(j-1.d0)*p3_8)/j

            enddo

            pp_8 = n*(z_8*p1_8-p2_8)/(z_8*z_8-1.d0)

            z1_8 = z_8
            z_8  = z1_8-p1_8/pp_8

         end do

         x(i)     = xm_8-xl_8*z_8
         x(n+1-i) = xm_8+xl_8*z_8
         w(i)     = 2.d0*xl_8/((1.d0-z_8*z_8)*pp_8*pp_8)
         w(n+1-i) = w(i)

      enddo

      !---------------------------------------------------------------

      return
      end
