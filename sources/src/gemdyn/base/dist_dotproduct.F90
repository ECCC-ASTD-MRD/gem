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
real*8 function dist_dotproduct(x, y, minx, maxx, miny, maxy, i0, il, j0, jl, nk) result(dotprod)
   ! Distributed dot product
   !
   !author
   !     Stéphane Gaudreault -- June 2014
   !
   !revision
   !     v4-80 - Gaudreault S.      - code refactoring
   !
   implicit none
#include <arch_specific.hf>

   integer, intent(in) :: minx, maxx, miny, maxy, i0, il, j0, jl, nk
   real*8, dimension(minx:maxx, miny:maxy, nk), intent(in) :: x, y

   real*8 :: local_dot
   integer :: i, j, k, ierr

   local_dot = 0.0d0
   do k=1,nk
      do j=j0,jl
         do i=i0,il
            local_dot = local_dot + (x(i, j, k) * y(i, j, k))
         end do
      end do
   end do

   call RPN_COMM_allreduce(local_dot, dotprod, 1, "MPI_double_precision", "MPI_sum", "grid", ierr)
end function dist_dotproduct
