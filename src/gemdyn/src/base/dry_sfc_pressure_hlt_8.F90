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

!**s/r dry_sfc_pressure_8 - Compute dry air surface pressure REAL64
!
      subroutine dry_sfc_pressure_hlt_8 (F_drysfcp0_8,presT_8,p0T_8,F_k0,F_timelevel_S)
      use dyn_fisl_options
      use glb_ld
      use tr3d
      use mem_tracers
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=1) :: F_timelevel_S
      real(kind=REAL64) F_drysfcp0_8(l_minx:l_maxx,l_miny:l_maxy),&
                             presT_8(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                               p0T_8(l_minx:l_maxx,l_miny:l_maxy)
      integer i,j,k,F_k0
      real, pointer, dimension(:    ) :: tr
      real, pointer, dimension(:,:,:) :: hu
!     ________________________________________________________________
!
      if (F_timelevel_S == 'P') then
         tr   => trt1
         hu   => tracers_P(Tr3d_hu)%pntr
      endif
      if (F_timelevel_S == 'M') then
         tr   => trt0
         hu   => tracers_M(Tr3d_hu)%pntr
      endif
      call sumhydro_hlt (sumq_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,Tr3d_ntr, tr, Schm_wload_L)

!$omp do collapse(2)
      do k=F_k0,l_nk
        do j=1+pil_s,l_nj-pil_n
          do i=1+pil_w,l_ni-pil_e
          sumq_8(i,j,k)= sumq_8(i,j,k) + hu(i,j,k)
          end do
        end do
      end do
!$omp end do

!$omp do
      do j=1+pil_s,l_nj-pil_n
         F_drysfcp0_8(:,j) = 0.0d0
         do k=F_k0,l_nk-1
            do i=1+pil_w,l_ni-pil_e
               F_drysfcp0_8(i,j)= F_drysfcp0_8(i,j) + &
                    (1.-sumq_8(i,j,k))*(presT_8(i,j,k+1) - presT_8(i,j,k))
            end do
         end do
         do i=1+pil_w,l_ni-pil_e
            F_drysfcp0_8(i,j)= F_drysfcp0_8(i,j) + &
                 (1.-sumq_8(i,j,l_nk))*(p0T_8(i,j) - presT_8(i,j,l_nk))
         end do
      end do
!$omp end do

!     ________________________________________________________________
!
      return
      end
