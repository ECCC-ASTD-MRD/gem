!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

!/@*
subroutine fill_halo (f,minx,maxx,miny,maxy,ni,nj,nk)
   implicit none
!!!#include <arch_specific.hf>
   !@Object 
   !@Arguments
   integer, intent(in) :: minx,maxx,miny,maxy,ni,nj,nk
   real, intent(inout) :: f(minx:maxx,miny:maxy,nk)
   !@Author 2018, Andre Plante
   !*@/
   integer :: i0,j0,i,j,k
   !----------------------------------------------------------------------
   i0 = 1
   j0 = 1
   do k=1,nk
      do i= minx , i0-1
         f(i,j0:nj,k)= f(i0  ,j0:nj,k)
      end do
      do i= ni+1, maxx
         f(i,j0:nj,k)= f(ni ,j0:nj,k)
      end do
   end do

   do k=1,nk
      do j= miny, j0-1
         f(minx:maxx,j,k)= f(minx:maxx,j0,k)
      end do
      do j= nj+1, maxy
         f(minx:maxx,j,k)= f(minx:maxx,nj,k)
      end do
   end do
   !----------------------------------------------------------------------
   return
end subroutine fill_halo
