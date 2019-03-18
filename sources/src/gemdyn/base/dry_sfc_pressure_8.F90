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

!**s/r dry_sfc_pressure_8 - Compute dry air surface pressure REAL*8
!
      subroutine dry_sfc_pressure_8 (F_drysfcp0_8, presT_8, p0T_8, &
                                     Minx,Maxx,Miny,Maxy,Nk,F_timelevel_S)
      use glb_ld
      use cstv
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

      character(len=1) :: F_timelevel_S
      integer Minx,Maxx,Miny,Maxy,Nk
      real*8 F_drysfcp0_8(Minx:Maxx,Miny:Maxy),presT_8(Minx:Maxx,Miny:Maxy,Nk),&
             p0T_8(Minx:Maxx,Miny:Maxy)

!author
!     Michel Desgagne --  fall 2014
!
!revision
! v4_70 - M. Desgagne      - Initial version
! v4_XX - M. Tanguay       - REAL*8


      integer i,j,k,istat
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: sumq
      real, pointer, dimension(:,:,:)         :: tr
!     ________________________________________________________________
!
      call sumhydro (sumq,Minx,Maxx,Miny,Maxy,Nk,F_timelevel_S)

      istat = gmm_get('TR/HU:'//F_timelevel_S,tr)

!$omp parallel private(i,j,k) &
!$omp shared(sumq,tr,F_drysfcp0_8,presT_8,p0T_8)

!$omp do
      do k=1,Nk
         sumq(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,k)= &
         sumq(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,k)+ &
         tr  (1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,k)
      end do
!$omp enddo
!$omp do
      do j=1+pil_s,l_nj-pil_n
         F_drysfcp0_8(:,j) = 0.0d0
         do k=1,Nk-1
         do i=1+pil_w,l_ni-pil_e
            F_drysfcp0_8(i,j)= F_drysfcp0_8(i,j) + &
                 (1.-sumq(i,j,k))*(presT_8(i,j,k+1) - presT_8(i,j,k))
         end do
         end do
         do i=1+pil_w,l_ni-pil_e
            F_drysfcp0_8(i,j)= F_drysfcp0_8(i,j) + &
                 (1.-sumq(i,j,Nk))*(p0T_8(i,j) - presT_8(i,j,Nk)) - Cstv_pref_8
         end do
      end do
!$omp enddo

!$omp end parallel

!     ________________________________________________________________
!
      return
      end
