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

!**s/r slabsym - symmetrical boundary conditions for theoretical cases
!
      subroutine slabsym ()
      use gmm_vt0
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer i,j,k,jin,jj
!
!----------------------------------------------------------------------
!
      if (l_north) then
!$omp do
         do k=1,G_nk
            jin = l_nj-pil_n-1
            jj  = l_nj-pil_n
            do i=1,l_ni
               vt0  (i,jj,k) = vt0  (i,jin,k)
            end do
            jin = l_nj-pil_n
            do j=1,pil_n
               jj  = l_nj-pil_n+j
               do i=1,l_ni
                  vt0 (i,jj,k) = vt0 (i,jin,k)
                  tt0 (i,jj,k) = tt0 (i,jin,k)
                  wt0 (i,jj,k) = wt0 (i,jin,k)
                  zdt0(i,jj,k) = zdt0(i,jin,k)
                  qt0 (i,jj,k) = qt0 (i,jin,k)
               end do
               do i=1,l_niu
                  ut0 (i,jj,k) = ut0 (i,jin,k)
               end do
            end do
         end do
!$omp end do
         jin = l_nj-pil_n
!$omp do
         do j=1,pil_n
            jj  = l_nj-pil_n+j
            do i=1,l_ni
               st0(i,jj)        = st0(i,jin)
               qt0(i,jj,G_nk+1) = qt0(i,jin,G_nk+1)
            end do
         end do
!$omp end do
      end if

      if (l_south) then
!$omp do
         do k=1,G_nk
            jin = pil_s+1
            do j=1,pil_s
               jj  = pil_s-j+1
               do i=1,l_ni
                  vt0 (i,jj,k) = vt0 (i,jin,k)
                  wt0 (i,jj,k) = wt0 (i,jin,k)
                  tt0 (i,jj,k) = tt0 (i,jin,k)
                  zdt0(i,jj,k) = zdt0(i,jin,k)
                  qt0 (i,jj,k) = qt0 (i,jin,k)
               end do
               do i=1,l_niu
                  ut0 (i,jj,k) = ut0 (i,jin,k)
               end do
            end do
         end do
!$omp end do
         jin = pil_s+1
!$omp do
         do j=1,pil_s
            jj  = pil_s-j+1
            do i=1,l_ni
               st0(i,jj)        = st0(i,jin)
               qt0(i,jj,G_nk+1) = qt0(i,jin,G_nk+1)
            end do
         end do
!$omp end do
      end if
!
!----------------------------------------------------------------------
!
      return
      end
