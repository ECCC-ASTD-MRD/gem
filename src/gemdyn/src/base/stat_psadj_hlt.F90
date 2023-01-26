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

!**s/r stat_psadj - Print Wet/Dry air mass at appropriate time

      subroutine stat_psadj_hlt(F_time,F_comment_S)

      use adz_mem
      use adz_options 
      use geomh
      use glb_ld
      use gmm_pw
      use HORgrid_options
      use lun
      use psadjust
      use masshlt
      use ptopo
      use omp_timing
      use omp_lib


      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,          intent(in) :: F_time       !Time 0 or Time 1
      character(len=*), intent(in) :: F_comment_S  !Comment

      !object
      !===============================================
      !     Print Wet/Dry air mass at appropriate time
      !===============================================

      include 'mpif.h'
      include 'rpn_comm.inc'
      integer :: err,i,j,i0_sb,in_sb,j0_sb,jn_sb,k0_sb,comm
      logical :: do_subset_GY_L
      real(kind=REAL64), dimension(:,:), pointer :: p0_wet_8,p0_dry_8
      real(kind=REAL64) :: l_avg_8(2),g_avg_ps_wet_8,g_avg_ps_dry_8,sum_lcl(2)
      real(kind=REAL64) :: gathV(2,Ptopo_numproc*Ptopo_ncolors)
      character(len= 7) :: time_S
      thread_sum=0.

!---------------------------------------------------------------------
!
      comm = RPN_COMM_comm ('MULTIGRID')

      p0_dry_8 => p0_0_8
      p0_wet_8 => p0_1_8

      if (F_time==1) time_S = "TIME T1"
      if (F_time==0) time_S = "TIME T0"

!$omp do
        do j=l_miny,l_maxy
          do i=l_minx,l_maxx
              p0_wet_8(i,j) = 0.
              p0_dry_8(i,j) = 0.
          end do
        end do
!$omp end do

      !Compute Dry Surface pressure at appropriate time
      !------------------------------------------------
      if (F_time==1) then
         call dry_sfc_pressure_hlt_8(p0_dry_8,pw_pm_plus_8,pw_p0_plus_8,Adz_k0t,'P')
      else if (F_time==0) then
         call dry_sfc_pressure_hlt_8(p0_dry_8,pw_pm_moins_8,pw_p0_moins_8,Adz_k0t,'M')
      endif  

      !Obtain Wet surface pressure/Pressure Momentum at appropriate time
      !-----------------------------------------------------------------
      if (F_time==1) then
       if (Adz_k0t==1) then
!$omp do
        do j=1+pil_s,l_nj-pil_n
          do i=1+pil_w,l_ni-pil_e
              p0_wet_8(i,j) = pw_p0_plus_8(i,j)
          end do
        end do
!$omp end do
       else if (Adz_k0t/=1) then
!$omp do
        do j=1+pil_s,l_nj-pil_n
          do i=1+pil_w,l_ni-pil_e
              p0_wet_8(i,j) = pw_pm_plus_8(i,j,l_nk+1) - pw_pm_plus_8(i,j,Adz_k0t)
          end do
        end do
!$omp end do
       endif
      else if (F_time==0) then
       if (Adz_k0t==1) then
!$omp do
        do j=1+pil_s,l_nj-pil_n
          do i=1+pil_w,l_ni-pil_e
              p0_wet_8(i,j) = pw_p0_moins_8(i,j)
          end do
        end do
!$omp end do
       else if (Adz_k0t/=1) then
!$omp do
        do j=1+pil_s,l_nj-pil_n
          do i=1+pil_w,l_ni-pil_e
              p0_wet_8(i,j) = pw_pm_moins_8(i,j,l_nk+1) - pw_pm_moins_8(i,j,Adz_k0t)
          end do
        end do
!$omp end do
       endif
      endif !F_TIME

      !Estimate Wet/Dry air mass on CORE
      !---------------------------------
      sum_lcl(1) = 0.0d0
      sum_lcl(2) = 0.0d0
!$omp do
      do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            sum_lcl(1) = sum_lcl(1) + p0_wet_8(i,j) * geomh_area_mask_8(i,j)
            sum_lcl(2) = sum_lcl(2) + p0_dry_8(i,j) * geomh_area_mask_8(i,j)
         end do
      end do
!$omp end do nowait
      thread_sum(1,OMP_get_thread_num()) = sum_lcl(1)
      thread_sum(2,OMP_get_thread_num()) = sum_lcl(2)
!$OMP BARRIER

!$omp single
      l_avg_8(1) = sum(thread_sum(1,:))
      l_avg_8(2) = sum(thread_sum(2,:))
      call MPI_Allgather(l_avg_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
      do i=1,2
         g_avg_8(i) = sum(gathV(i,:))
      end do
!$omp end single

      g_avg_ps_wet_8 = g_avg_8(1) * PSADJ_scale_8
      g_avg_ps_dry_8 = g_avg_8(2) * PSADJ_scale_8

!$omp single
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ""
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "C:WET air mass",time_S,g_avg_ps_wet_8,F_comment_S
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "C:DRY air mass",time_S,g_avg_ps_dry_8,F_comment_S
!$omp end single

      do_subset_GY_L = Grd_yinyang_L.and.adz_pil_sub_s_g /= -1

      !Estimate Wet/Dry air mass on SUBSET of YIN
      !------------------------------------------
      if (do_subset_GY_L) then

         i0_sb = 1+Adz_pil_sub_w ; j0_sb = 1+Adz_pil_sub_s ; in_sb = l_ni-Adz_pil_sub_e ; jn_sb = l_nj-Adz_pil_sub_n ; k0_sb = Adz_k0t_sub

         !Compute Dry Surface pressure at appropriate time
         !------------------------------------------------
         if (F_time==1) then
           call dry_sfc_pressure_hlt_8(p0_dry_8,pw_pm_plus_8,pw_p0_plus_8,k0_sb,'P')
         else if (F_time==0) then
           call dry_sfc_pressure_hlt_8(p0_dry_8,pw_pm_moins_8,pw_p0_moins_8,k0_sb,'M')
         endif  

         !Obtain Wet surface pressure/Pressure Momentum at appropriate time
         !-----------------------------------------------------------------
         if (k0_sb/=1) then
           if (F_time==1) then
!$omp do
            do j=j0_sb,jn_sb
             do i=i0_sb,in_sb
               p0_wet_8(i,j) = pw_pm_plus_8(i,j,l_nk+1) - pw_pm_plus_8(i,j,k0_sb)
             end do
            end do
!$omp end do
           else if (F_time==0) then
!$omp do
            do j=j0_sb,jn_sb
             do i=i0_sb,in_sb
               p0_wet_8(i,j) = pw_pm_moins_8(i,j,l_nk+1) - pw_pm_moins_8(i,j,k0_sb)
             end do
            end do
!$omp end do
           endif
         endif !F_TIME

         sum_lcl(1) = 0.0d0
         sum_lcl(2) = 0.0d0
!$omp do
         do j=j0_sb,jn_sb
           do i=i0_sb,in_sb
            sum_lcl(1) = sum_lcl(1) + p0_wet_8(i,j) * geomh_area_mask_8(i,j)
            sum_lcl(2) = sum_lcl(2) + p0_dry_8(i,j) * geomh_area_mask_8(i,j)
           end do
         end do
!$omp end do nowait
         thread_sum(1,OMP_get_thread_num()) = sum_lcl(1)
         thread_sum(2,OMP_get_thread_num()) = sum_lcl(2)
!$OMP BARRIER

!$omp single
         l_avg_8(1) = sum(thread_sum(1,:))
         l_avg_8(2) = sum(thread_sum(2,:))
         if (Ptopo_couleur/=0) l_avg_8(1:2) = 0.0d0
         call MPI_Allgather(l_avg_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
         do i=1,2
           g_avg_8(i) = sum(gathV(i,:))
         end do
!$omp end single

         g_avg_ps_wet_8 = g_avg_8(1) / Adz_gs_area_8
         g_avg_ps_dry_8 = g_avg_8(2) / Adz_gs_area_8

!$omp single
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ""
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1001) "SB:WET air mass",time_S,g_avg_ps_wet_8,F_comment_S
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1001) "SB:DRY air mass",time_S,g_avg_ps_dry_8,F_comment_S
!$omp end single

      end if

!$omp do
        do j=l_miny,l_maxy
          do i=l_minx,l_maxx
              p0_wet_8(i,j) = 0.
              p0_dry_8(i,j) = 0.
          end do
        end do
!$omp end do

      nullify(p0_dry_8,p0_wet_8)
!---------------------------------------------------------------------
!
      return

 1000 format(1X,A14,1X,A7,1X,E19.12,1X,A16)
 1001 format(1X,A15,1X,A7,1X,E19.12,1X,A16)

      end subroutine stat_psadj_hlt
