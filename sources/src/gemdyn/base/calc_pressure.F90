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

!**s/r calc_pressure - Compute pressure on momentum and thermodynamic
!                      levels and also surface pressure.
!
      subroutine calc_pressure (F_pm, F_pt, F_p0, F_s, Minx,Maxx,Miny,Maxy,Nk)
      use gmm_vt1
      use gmm_geof
      use glb_ld
      use ver
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

      integer Minx,Maxx,Miny,Maxy,Nk
      real F_pm(Minx:Maxx,Miny:Maxy,Nk), F_pt(Minx:Maxx,Miny:Maxy,Nk), &
           F_p0(Minx:Maxx,Miny:Maxy)   , F_s (Minx:Maxx,Miny:Maxy)

      integer i, j, k, istat
      real wk1(l_ni,l_nj,2),wk2(l_ni,l_nj,2)
!     ________________________________________________________________
!
      istat = gmm_get (gmmk_sls_s ,sls )

!$omp parallel private(wk1,wk2,i,j)
!$omp do
      do k=1,l_nk+1
         do j=1,l_nj
         do i=1,l_ni
            wk1(i,j,1) = Ver_a_8%m(k) + Ver_b_8%m(k) * F_s(i,j) + Ver_c_8%m(k) * sls(i,j)
            wk1(i,j,2) = Ver_a_8%t(k) + Ver_b_8%t(k) * F_s(i,j) + Ver_c_8%t(k) * sls(i,j)
         end do
         end do
         call vsexp(wk2(1,1,1), wk1(1,1,1), 2*l_ni*l_nj)
         if (k < l_nk+1) then
            F_pm(1:l_ni,1:l_nj,k) = wk2(1:l_ni,1:l_nj,1)
            F_pt(1:l_ni,1:l_nj,k) = wk2(1:l_ni,1:l_nj,2)
         else
            F_p0(1:l_ni,1:l_nj)= wk2(1:l_ni,1:l_nj,1)
         end if
      end do
!$omp enddo
!$omp end parallel
!     ________________________________________________________________
!
      return
      end
