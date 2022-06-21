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

!**s/r sumhydro - Sum over Hydrometeors
!
      subroutine sumhydro_hlt (F_qh,minx,maxx,miny,maxy,nk,ntr,F_tracers)
      use dyn_fisl_options
      use gem_options
      use glb_ld
      use tr3d
      use mem_tracers
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: minx,maxx,miny,maxy,nk,ntr
      real, intent(in ) :: F_tracers(minx:maxx,miny:maxy,nk,ntr)
      real, intent(out) :: F_qh(minx:maxx,miny:maxy,nk)

      integer i, j, k, n
!     ________________________________________________________________
!
!$omp do collapse(2)
      do k = 1, l_nk
        do j= 1-G_haloy, l_nj+G_haloy
          do i= 1-G_halox, l_ni+G_halox
             F_qh(i,j,k)=0.
          end do
        end do
      end do
!$omp end do

      if (.not.Schm_wload_L) return

!     Sum over Hydrometeors
      do n = 1, Tr3d_ntr
         if (Tr3d_wload (n)) then
!$omp do collapse(2)
            do k = 1, l_nk
               do j= 1-G_haloy, l_nj+G_haloy
                  do i= 1-G_halox, l_ni+G_halox
                     F_qh(i,j,k)=F_qh(i,j,k)+F_tracers(i,j,k,n)
                  end do
               end do
            end do
!$omp end do
         end if
      end do
!     ________________________________________________________________
!
      return
      end subroutine sumhydro_hlt
