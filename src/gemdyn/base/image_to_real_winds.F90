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

!**s/r image_to_real_winds - Convert image winds to real winds

   subroutine image_to_real_winds( F_u,F_v,Minx,Maxx,Miny,Maxy,F_nk )
      use dcst
      use geomh
      use tdpack
      use glb_ld
   implicit none
#include <arch_specific.hf>

   integer                                  , intent(IN)     :: Minx,Maxx,Miny,Maxy,F_nk
   real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(IN OUT) :: F_u, F_v


   ! Local variables
   integer :: i,j,k
   real :: c1

!$omp parallel private(c1)
!$omp do
   do k = 1, F_nk
      do j = 1, l_nj
         c1 = Dcst_rayt_8 / geomh_cy_8(j)
         do i = 1, l_niu
            F_u(i,j,k) = c1 * F_u(i,j,k)
         end do
      end do
   end do
!$omp enddo
!$omp do
   do k = 1, F_nk
      do j = 1, l_njv
         c1 = Dcst_rayt_8 / geomh_cyv_8(j)
         do i = 1, l_ni
            F_v(i,j,k) = c1 * F_v(i,j,k)
         end do
      end do
   end do
!$omp enddo
!$omp end parallel

   return

end subroutine image_to_real_winds
