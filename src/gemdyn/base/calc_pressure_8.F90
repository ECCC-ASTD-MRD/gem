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

!**s/r calc_pressure_8 - see calc_pressure
!
      subroutine calc_pressure_8 ( F_pm, F_pt, F_p0, F_s, &
                                   Minx,Maxx,Miny,Maxy, Nk )
      use gmm_vt1
      use gmm_geof
      use glb_ld
      use ver
      use gmm_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,Nk
      real, dimension(Minx:Maxx,Miny:Maxy   ), intent(IN ) :: F_s
      real (kind=REAL64),dimension(Minx:Maxx,Miny:Maxy   ),intent(OUT)::&
                                                                    F_p0
      real (kind=REAL64),dimension(Minx:Maxx,Miny:Maxy,Nk),intent(OUT)::&
                                                              F_pm, F_pt
      integer i, j, k, istat
!     ________________________________________________________________
!
      istat = gmm_get (gmmk_sls_s ,sls )

!$omp parallel
!$omp do
      do k=1,l_nk
         do j=1,l_nj
         do i=1,l_ni
            F_pm(i,j,k) = exp( Ver_a_8%m(k) + Ver_b_8%m(k) * F_s(i,j)&
                              +Ver_c_8%m(k) * sls(i,j))
            F_pt(i,j,k) = exp( Ver_a_8%t(k) + Ver_b_8%t(k) * F_s(i,j)&
                              +Ver_c_8%t(k) * sls(i,j))
         end do
         end do
      end do
!$omp enddo

      k= l_nk+1
!$omp do
      do j=1,l_nj
         do i=1,l_ni
            F_p0(i,j)= exp( Ver_a_8%m(k) + Ver_b_8%m(k) * F_s(i,j)&
                           +Ver_c_8%m(k) * sls(i,j))
        end do
      end do
!$omp enddo
!$omp end parallel
!     ________________________________________________________________
!
      return
      end subroutine calc_pressure_8
