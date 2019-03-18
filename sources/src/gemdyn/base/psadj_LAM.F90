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

!**s/r psadj_LAM - Adjust surface pressure for conservation using Flux calculations based on Aranami et al. (2015)

      subroutine psadj_LAM (F_cub_o,F_cub_i,Minx,Maxx,Miny,Maxy,F_nk,F_k0)

      use cstv
      use dyn_fisl_options
      use geomh
      use glb_ld
      use gmm_vt1
      use gmm_vt0
      use gmm_itf_mod
      use lun

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: Minx,Maxx,Miny,Maxy                         !Dimension H
      integer, intent(in) :: F_k0                                        !Scope of operator
      integer, intent(in) :: F_nk                                        !Number of vertical levels
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(in)   :: F_cub_o !High-order SL solution FLUX_out
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(in)   :: F_cub_i !High-order SL solution FLUX_in

      !object
      !=====================================================================
      !     Adjust surface pressure for conservation using Flux calculations
      !     based on Aranami et al.,QJRMS,141,1795-1803 (2015)
      !=====================================================================

      !---------------------------------------------------------------------

      real, pointer, contiguous, dimension(:,:,:) :: tr
      real*8,dimension(Minx:Maxx,Miny:Maxy,F_nk+1):: pr_m_8,pr_t_8
      real*8,dimension(Minx:Maxx,Miny:Maxy)       :: pr_p0_1_8,pr_p0_0_8,pr_fl_w_8,pr_p0_w_8
      real,  dimension(Minx:Maxx,Miny:Maxy,F_nk)  :: sumq
      integer :: err,i,j,k,n,istat,MAX_iteration
      real*8  :: l_avg_8(2),g_avg_8(2),g_avg_ps_w_1_8,g_avg_ps_w_0_8,g_avg_fl_w_0_8
      logical, save :: done_LAM_scale_L = .false.
      real*8,  save :: PSADJ_LAM_scale_8

      !---------------------------------------------------------------------

      if (Lun_out>0) then
         write(Lun_out,*) ''
         write(Lun_out,*) '--------------------------------------'
         if (Schm_psadj==2) write(Lun_out,*) 'PSADJ LAM is done for DRY AIR (REAL*8)'
         if (Schm_psadj==1) write(Lun_out,*) 'PSADJ LAM is done for WET AIR (REAL*8)'
         write(Lun_out,*) '--------------------------------------'
      end if

      !Estimate area
      !-------------
      if (.not.done_LAM_scale_L) then

         l_avg_8(1) = 0.0d0
         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            l_avg_8(1) = l_avg_8(1) + geomh_area_8(i,j) * geomh_mask_8(i,j)
         end do
         end do

         call RPN_COMM_allreduce (l_avg_8(1), PSADJ_LAM_scale_8, 1, &
           "MPI_DOUBLE_PRECISION","MPI_SUM","GRID",err)

         PSADJ_LAM_scale_8 = 1.0d0/PSADJ_LAM_scale_8

         done_LAM_scale_L = .true.

      end if

      !-------------
      !Treat TIME T1
      !-------------

      istat = gmm_get(gmmk_st1_s,st1)

      !Obtain pressure levels
      !----------------------
      call calc_pressure_8 (pr_m_8,pr_t_8,pr_p0_1_8,st1,Minx,Maxx,Miny,Maxy,F_nk)

      pr_m_8(:,:,F_nk+1) = pr_p0_1_8(:,:)

      !Compute wet (psadj=1)/ dry (psadj=2) surface pressure on CORE
      !-------------------------------------------------------------
      if (Schm_psadj==2) then

         call sumhydro (sumq,Minx,Maxx,Miny,Maxy,F_nk,'P')

         istat = gmm_get('TR/HU:P',tr)

!$omp parallel do private(k) shared(sumq,tr)
         do k=F_k0,F_nk
            sumq(1:l_ni,1:l_nj,k) = sumq(1:l_ni,1:l_nj,k) + tr(1:l_ni,1:l_nj,k)
         end do
!$omp end parallel do

      else

         sumq = 0.

      end if

      if (Schm_psadj==1) pr_p0_w_8(:,:) = pr_m_8(:,:,1)
      if (Schm_psadj==2) pr_p0_w_8(:,:) = 0.0d0

!$omp parallel private(i,j,k) &
!$omp shared(pr_m_8,pr_p0_w_8,sumq)
!$omp do
      do j=1+pil_s,l_nj-pil_n
         do k=F_k0,F_nk
            do i=1+pil_w,l_ni-pil_e
               pr_p0_w_8(i,j) = pr_p0_w_8(i,j) + (1.-sumq(i,j,k)) * (pr_m_8(i,j,k+1)-pr_m_8(i,j,k))
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      !Estimate air mass on CORE
      !-------------------------
      l_avg_8(1) = 0.0d0

      do j=1+pil_s,l_nj-pil_n
      do i=1+pil_w,l_ni-pil_e

         l_avg_8(1) = l_avg_8(1) + pr_p0_w_8(i,j) * geomh_area_8(i,j) * geomh_mask_8(i,j)

      end do
      end do

      call RPN_COMM_allreduce (l_avg_8,g_avg_8,1,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",err)

      g_avg_ps_w_1_8 = g_avg_8(1) * PSADJ_LAM_scale_8

      MAX_iteration = 3
      if (Schm_psadj==1) MAX_iteration = 1

      do n = 1,MAX_iteration

         !-------------
         !Treat TIME T0
         !-------------

         istat = gmm_get(gmmk_st0_s,st0)

         !Obtain pressure levels
         !----------------------
         call calc_pressure_8 (pr_m_8,pr_t_8,pr_p0_0_8,st0,Minx,Maxx,Miny,Maxy,F_nk)

         pr_m_8(:,:,F_nk+1) = pr_p0_0_8(:,:)

         !Compute wet (psadj=1)/ dry (psadj=2) surface pressure on CORE
         !-------------------------------------------------------------
         if (Schm_psadj==2) then

            call sumhydro (sumq,Minx,Maxx,Miny,Maxy,F_nk,'M')

            istat = gmm_get('TR/HU:M',tr)

!$omp parallel do private(k) shared(sumq,tr)
            do k=F_k0,F_nk
               sumq(1:l_ni,1:l_nj,k) = sumq(1:l_ni,1:l_nj,k) + tr(1:l_ni,1:l_nj,k)
            end do
!$omp end parallel do

         else

            sumq = 0.

         end if

         if (Schm_psadj==1) pr_p0_w_8(:,:) = pr_m_8(:,:,1)
         if (Schm_psadj==2) pr_p0_w_8(:,:) = 0.0d0

!$omp parallel private(i,j,k) &
!$omp shared(pr_m_8,pr_p0_w_8,sumq)
!$omp do
         do j=1+pil_s,l_nj-pil_n
            do k=F_k0,F_nk
               do i=1+pil_w,l_ni-pil_e
                  pr_p0_w_8(i,j) = pr_p0_w_8(i,j) + (1.-sumq(i,j,k)) * (pr_m_8(i,j,k+1)-pr_m_8(i,j,k))
               end do
            end do
         end do
!$omp end do
!$omp end parallel

         !Compute FLUX on NEST+CORE
         !-------------------------
!$omp parallel private(i,j,k) &
!$omp shared(pr_m_8,pr_fl_w_8,sumq,F_cub_i,F_cub_o)
!$omp do
         do j=1,l_nj
            pr_fl_w_8(:,j) = 0.0d0
            do k=F_k0,F_nk
            do i=1,l_ni
               pr_fl_w_8(i,j) = pr_fl_w_8(i,j) + (1.-sumq(i,j,k)) * (pr_m_8(i,j,k+1) - pr_m_8(i,j,k)) &
                                * (F_cub_i(i,j,k) - F_cub_o(i,j,k))
            end do
            end do
         end do
!$omp end do
!$omp end parallel

         !Estimate air mass on CORE
         !-------------------------
         l_avg_8(1) = 0.0d0

         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e

            l_avg_8(1) = l_avg_8(1) + pr_p0_w_8(i,j) * geomh_area_8(i,j) * geomh_mask_8(i,j)

         end do
         end do

         !Estimate FLUX mass on NEST+CORE
         !-------------------------------
         l_avg_8(2) = 0.0d0

         do j=1,l_nj
         do i=1,l_ni

            l_avg_8(2) = l_avg_8(2) + pr_fl_w_8(i,j) * geomh_area_8(i,j) * geomh_mask_8(i,j)

         end do
         end do

         call RPN_COMM_allreduce (l_avg_8,g_avg_8,2,"MPI_DOUBLE_PRECISION","MPI_SUM","GRID",err)

         g_avg_ps_w_0_8 = g_avg_8(1) * PSADJ_LAM_scale_8
         g_avg_fl_w_0_8 = g_avg_8(2) * PSADJ_LAM_scale_8

         !Correct surface pressure in order to preserve air mass on CORE when taking FLUX mass into account
         !-------------------------------------------------------------------------------------------------
         pr_p0_0_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n) = pr_p0_0_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n) + &
                                                            (g_avg_ps_w_1_8 - g_avg_ps_w_0_8) + g_avg_fl_w_0_8

         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            st0(i,j)= log(pr_p0_0_8(i,j)/Cstv_pref_8)
         end do
         end do

      end do

      !-----------
      !Diagnostics
      !-----------

      if (Lun_out>0) then
         write(Lun_out,*)    ''
         write(Lun_out,*)    '---------------------------------------------------------'
         write(Lun_out,1004) 'PSADJ_LAM: C=',g_avg_ps_w_0_8,' FLUX=',g_avg_fl_w_0_8
         write(Lun_out,*)    '---------------------------------------------------------'
         write(Lun_out,*)    ''
      end if

      !---------------------------------------------------------------------

      return

1004 format(1X,A13,E19.12,A6,E19.12)

      end subroutine psadj_LAM
