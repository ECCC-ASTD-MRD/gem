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

   subroutine tricub_zyx3 (src, x, y, abcd, ni, ninj, xyz)
   implicit none
   integer, intent(in) :: ni, ninj
   real, dimension(3,0:*), intent(in) :: src
   real*8, dimension(0:3), intent(in) :: abcd,x,y
   real, dimension(3), intent(out) :: xyz

   integer :: i, j, ni2, ni3, ninj2, ninj3
   real*8, dimension(0:3) :: va4, vb4, vc4, vd4, dst

   ninj2 = ninj + ninj
   ninj3 = ninj2 + ninj
   ni2 = ni + ni
   ni3 = ni2 + ni

   do j = 1,3

   do i = 0,3
      va4(i) = src(j,i    )*abcd(0) + src(j,i    +ninj)*abcd(1) +  src(j,i    +ninj2)*abcd(2) + src(j,i    +ninj3)*abcd(3)
      vb4(i) = src(j,i+ni )*abcd(0) + src(j,i+ni +ninj)*abcd(1) +  src(j,i+ni +ninj2)*abcd(2) + src(j,i+ni +ninj3)*abcd(3)
      vc4(i) = src(j,i+ni2)*abcd(0) + src(j,i+ni2+ninj)*abcd(1) +  src(j,i+ni2+ninj2)*abcd(2) + src(j,i+ni2+ninj3)*abcd(3)
      vd4(i) = src(j,i+ni3)*abcd(0) + src(j,i+ni3+ninj)*abcd(1) +  src(j,i+ni3+ninj2)*abcd(2) + src(j,i+ni3+ninj3)*abcd(3)
      dst(i) = va4(i)*y(0) + vb4(i)*y(1) + vc4(i)*y(2) + vd4(i)*y(3)
   end do
   xyz(j) = dst(0)*x(0) + dst(1)*x(1) + dst(2)*x(2) + dst(3)*x(3)

   end do

  return
end
