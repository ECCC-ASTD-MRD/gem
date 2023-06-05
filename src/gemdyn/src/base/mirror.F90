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

!**s/r mirror - mirror boundary condition for BUBBLE theoretical case
!
      subroutine mirror ()
      use glb_ld
      use gmm_vt0
      implicit none
#include <arch_specific.hf>

      integer i, j, k, ii, iin, jj, jin
!
!----------------------------------------------------------------------
!
      if (l_north) then
!$omp do
         do k=1,G_nk
            do i=1,l_ni
               vt0  (i,l_nj-pil_n,k) = 0.
            end do
            do j=1,pil_n-1
               jin = l_nj-pil_n-j
               jj = l_nj-pil_n+j
               do i=1,l_ni
                  vt0  (i,jj,k) = - vt0  (i,jin,k)
               end do
            end do
            do j=1,pil_n
               jin = l_nj-pil_n-j+1
               jj  = l_nj-pil_n+j
               do i=1,l_ni
                  tt0 (i,jj,k) = tt0 (i,jin,k)
                  zdt0(i,jj,k) = zdt0(i,jin,k)
                  wt0 (i,jj,k) = wt0 (i,jin,k)
                  qt0 (i,jj,k) = qt0 (i,jin,k)
               end do
               do i=1,l_niu
                  ut0   (i,jj,k) = ut0   (i,jin,k)
               end do
            end do
         end do
!$omp end do
!$omp do
         do j=1,pil_n
            jin = l_nj-pil_n-j+1
            jj  = l_nj-pil_n+j
            do i=1,l_ni
               st0(i,jj)        = st0(i,jin)
               qt0(i,jj,G_nk+1) = qt0(i,jin,G_nk+1)
            end do
         end do
!$omp end do
      end if

      if (l_east) then
!$omp do
         do k=1,G_nk
            do j=1,l_nj
               ut0  (l_ni-pil_e,j,k) = 0.
            end do
            do j=1,l_nj
               do i=1,pil_e-1
                  iin = l_ni-pil_e-i
                  ii  = l_ni-pil_e+i
                  ut0  (ii,j,k) = - ut0  (iin,j,k)
               end do
            end do
            do j=1,l_nj
               do i=1,pil_e
                  iin = l_ni-pil_e-i+1
                  ii  = l_ni-pil_e+i
                  tt0 (ii,j,k) = tt0 (iin,j,k)
                  zdt0(ii,j,k) = zdt0(iin,j,k)
                  wt0 (ii,j,k) = wt0 (iin,j,k)
                  qt0 (ii,j,k) = qt0 (iin,j,k)
               end do
            end do
            do j=1,l_njv
               do i=1,pil_e
                  iin = l_ni-pil_e-i+1
                  ii  = l_ni-pil_e+i
                  vt0  (ii,j,k) = vt0  (iin,j,k)
               end do
            end do
         end do
!$omp end do
!$omp do
         do j=1,l_nj
            do i=1,pil_e
               iin = l_ni-pil_e-i+1
               ii  = l_ni-pil_e+i
               st0(ii,j)        = st0(iin,j)
               qt0(ii,j,G_nk+1) = qt0(iin,j,G_nk+1)
            end do
         end do
!$omp end do
      end if

      if (l_south) then
!$omp do
         do k=1,G_nk
            do i=1,l_ni
               vt0  (i,pil_s,k) = 0.
            end do
            do j=1,pil_s-1
               jin = pil_s+j
               jj  = pil_s-j
               do i=1,l_ni
                  vt0  (i,jj,k) = - vt0 (i,jin,k)
               end do
            end do
            do j=1,pil_s
               jin = pil_s+j
               jj  = pil_s-j+1
               do i=1,l_ni
                  tt0 (i,jj,k) = tt0 (i,jin,k)
                  zdt0(i,jj,k) = zdt0(i,jin,k)
                  wt0 (i,jj,k) = wt0 (i,jin,k)
                  qt0 (i,jj,k) = qt0 (i,jin,k)
               end do
               do i=1,l_niu
                  ut0  (i,jj,k) = ut0  (i,jin,k)
               end do
            end do
         end do
!$omp end do
!$omp do
         do j=1,pil_s
            jin = pil_s+j
            jj  = pil_s-j+1
            do i=1,l_ni
               st0(i,jj)        = st0(i,jin)
               qt0(i,jj,G_nk+1) = qt0(i,jin,G_nk+1)
            end do
         end do
!$omp end do
      end if

      if (l_west) then
!$omp do
         do k=1,G_nk
            do j=1,l_nj
               ut0  (pil_w,j,k) = 0.
            end do
            do j=1,l_nj
               do i=1,pil_w-1
                  iin = pil_w+i
                  ii  = pil_w-i
                  ut0  (ii,j,k) = - ut0  (iin,j,k)
               end do
            end do
            do j=1,l_nj
               do i=1,pil_w
                  iin = pil_w+i
                  ii  = pil_w-i+1
                  tt0 (ii,j,k) = tt0 (iin,j,k)
                  zdt0(ii,j,k) = zdt0(iin,j,k)
                  wt0 (ii,j,k) = wt0 (iin,j,k)
                  qt0 (ii,j,k) = qt0 (iin,j,k)
               end do
            end do
            do j=1,l_njv
               do i=1,pil_w
                  iin = pil_w+i
                  ii  = pil_w-i+1
                  vt0  (ii,j,k) = vt0  (iin,j,k)
               end do
            end do
         end do
!$omp end do
!$omp do
         do j=1,l_nj
            do i=1,pil_w
               iin = pil_w+i
               ii  = pil_w-i+1
               st0(ii,j)        = st0(iin,j)
               qt0(ii,j,G_nk+1) = qt0(iin,j,G_nk+1)
            end do
         end do
!$omp end do
      end if
!
!----------------------------------------------------------------------
!
      return
      end
