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
!
   function linz_cubyx(src, x, y, z, ni, ninj)  result (xyz)
   implicit none
   real, dimension(0:*), intent(in) :: src
   real*8, dimension(0:3), intent(in) :: x, y
   real*8, intent(in) :: z
   real :: xyz
   integer :: i,ni, ninj, ni2, ni3
   real*8, dimension(0:3) :: va4, vb4, vc4, vd4
   real , dimension(0:3) ::  dst
   real*8 :: z0, z1

   ni2 = ni + ni
   ni3 = ni2 + ni

   z1=z
   z0=1. - z

   do i = 0,3
      va4(i) =  src(i    )*z0 +  src(i    +ninj)*z1
      vb4(i) =  src(i+ni )*z0 +  src(i+ni +ninj)*z1
      vc4(i) =  src(i+ni2)*z0 +  src(i+ni2+ninj)*z1
      vd4(i) =  src(i+ni3)*z0 +  src(i+ni3+ninj)*z1
      dst(i) = va4(i)*y(0) + vb4(i)*y(1) + vc4(i)*y(2) + vd4(i)*y(3)
   end do

   xyz  = dst(0)*x(0) + dst(1)*x(1) + dst(2)*x(2) + dst(3)*x(3)

   return
   end
