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

!**s/r vterp - Finds the x point in f
!
   subroutine vterp (F_x, F_fx, F_f, F_g, F_y, F_acc, np, nn)
   implicit none
#include <arch_specific.hf>
   !
   integer, intent(in) :: nn, np
   real, dimension(np), intent(out) :: F_x
   real, dimension(np), intent(in) :: F_fx
   real, dimension(nn), intent(in) :: F_y
   real, dimension(np,nn), intent(in) :: F_f, F_g
   real, intent(in) :: F_acc
   !
   !OBJECT
   !      Given the values of a monotonic function F and the values of its
   !      derivative G at NN points Y(1) TO Y(NN), this routine finds X
   !      point at which F assumes the specified value FX.
   !      AT INPUT A FIRST GUESS SHOULD BE PROVIDED FOR X.
   !      WE ASSUME FX LE F(1)
   !
   !ARGUMENTS
   !  Name                      Description
   !----------------------------------------------------------
   ! F_X               X point to find
   ! F_F               function
   ! F_FX              value at one point of the function
   ! F_G               slopes of points
   ! F_Y               value of points
   ! F_ACC             accuracy requested for the interpolation
   !----------------------------------------------------------------
   !
   !*
   integer :: n, i
   real :: f0, g0, y0, dy, a, c, p, er, der, root, &
           f1, g1, y1, cd, b, r, q

   do i=1,np

      if (F_fx(i) == F_f(i,nn)) then
         F_x(i) = F_y(nn)
      elseif (F_fx(i) < F_f(i,nn)) then
         ! extrapolation
         f1 = F_f(i,nn-1)
         f0 = F_f(i,nn)
         g1 = F_g(i,nn-1)
         g0 = F_g(i,nn)
         y1 = F_y(nn-1)
         y0 = F_y(nn)
         root = g0**2 - 2.*(g0-g1)*(f0-F_fx(i))/(y0-y1)
         if (root >= 0.0 .and. abs((g0-g1)/g0) > 0.01) then
            ! use quadratic extrapolation
            root = sqrt(root)
            F_x(i)    = y0-(y0-y1)*(g0+root)/(g0-g1)
         else
            ! use linear formula
            F_x(i)    = y0+(F_fx(i)-f0)/g0
         end if
      elseif (F_fx(i) > F_f(i,1)) then
         stop 'vterp : valeur a interpoler trop grande'
      else

         if (F_fx(i) > F_f(i,n)) then
            f1 = F_f(i,n-1)
            f0 = F_f(i,n)
            g1 = F_g(i,n-1)
            g0 = F_g(i,n)
            y1 = F_y(n-1)
            y0 = F_y(n)
            dy = y1-y0
            a  = +f1/dy
            b  = -f0/dy
            c  = (g0+g1)/dy**2 - 2.*(f1-f0)/dy**3
            cd = (y1*g0+y0*g1-(a+b)*(y1+y0))/dy**2
            ! Newton formula iteration loop
            do
               p  = F_x(i)-y0
               q  = F_x(i)-y1
               r  = c*F_x(i)-cd
               er = a*p+b*q+p*q*r-F_fx(i)
               if (abs(er) < F_acc) exit
               der = a+b+p*r+q*r+c*p*q
               F_x(i)  = F_x(i) -er/der
            end do
         end if

      end if

   end do

   return
   end
