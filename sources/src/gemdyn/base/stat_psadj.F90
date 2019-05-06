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

!**s/r stat_psadj - Print dry/wet air mass

      subroutine stat_psadj (F_time,F_comment_S)

      use cstv
      use dynkernel_options
      use geomh
      use glb_ld
      use glb_pil
      use gmm_itf_mod
      use gmm_vt1
      use gmm_vt0
      use gem_options
      use HORgrid_options
      use lun
      use metric
      use ptopo
      use tdpack, only: rgasd_8
      use ver

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,          intent(in) :: F_time       !Time 0 or Time 1
      character(len=*), intent(in) :: F_comment_S  !Comment

      !object
      !===========================
      !     Print dry/wet air mass
      !===========================

      !---------------------------------------------------------------

      real, pointer, dimension (:,:)   :: w2d
      real, pointer, dimension (:,:,:) :: wqt
      integer err,i,j,istat
      real*8,dimension(l_minx:l_maxx,l_miny:l_maxy,1:l_nk):: pr_m_8,pr_t_8
      real*8,dimension(l_minx:l_maxx,l_miny:l_maxy)       :: pr_p0_8,pr_p0_dry_8,log_p0_8
      real*8 l_avg_8(2),g_avg_ps_8(2),scale_8
      real*8 ll_avg_8(l_ni,l_nj,2),lx_avg_8(l_ni,l_nj)
      character(len= 9) communicate_S
      character(len= 1) in_S
      character(len= 7) time_S

      !---------------------------------------------------------------

      !Estimate area
      !-------------
      l_avg_8(1) = 0.0d0
      lx_avg_8 = 0.0d0
      do j=1+pil_s,l_nj-pil_n
      do i=1+pil_w,l_ni-pil_e
         l_avg_8(1) = l_avg_8(1) + geomh_area_8(i,j) * geomh_mask_8(i,j)
         lx_avg_8(i,j) =           geomh_area_8(i,j) * geomh_mask_8(i,j)
      end do
      end do

      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"
      if ( Lctl_rxstat_S == 'GLB_8') then
         call glbsum8 (l_avg_8(1),lx_avg_8,1,l_ni,1,l_nj,1,          &
                   1+Lam_pil_w, G_ni-Lam_pil_e, 1+Lam_pil_s, G_nj-Lam_pil_n)
         if (Grd_yinyang_L) communicate_S = "GRIDPEERS"
      end if

      call RPN_COMM_allreduce (l_avg_8, scale_8, 1, &
        "MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      scale_8 = 1.0d0/scale_8

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then
         if (F_time==1) istat = gmm_get(gmmk_st1_s,w2d)
         if (F_time==0) istat = gmm_get(gmmk_st0_s,w2d)
      elseif (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then
         if (F_time==1) istat = gmm_get(gmmk_qt1_s,wqt)
         if (F_time==0) istat = gmm_get(gmmk_qt0_s,wqt)
      end if

      if (F_time==1) in_S = 'P'
      if (F_time==0) in_S = 'M'

      if (F_time==1) time_S = "TIME T1"
      if (F_time==0) time_S = "TIME T0"

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

         !Obtain pressure levels
         !----------------------
         call calc_pressure_8 ( pr_m_8, pr_t_8, pr_p0_8, w2d, &
                                l_minx,l_maxx,l_miny,l_maxy,l_nk )

         !Compute dry surface pressure (- Cstv_pref_8)
         !--------------------------------------------
         call dry_sfc_pressure_8 ( pr_p0_dry_8, pr_m_8, pr_p0_8, &
                                   l_minx,l_maxx,l_miny,l_maxy,l_nk,in_S)

         !Add Cstv_pref_8 to dry surface pressure
         !---------------------------------------
         pr_p0_dry_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n) = pr_p0_dry_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n) + Cstv_pref_8

      elseif (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

         !Obtain surface pressure
         !-----------------------
         log_p0_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n) = (dble(wqt(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,l_nk+1))/&
                                                           (rgasd_8*Ver_Tstar_8%m(l_nk+1))  &
                                                          + dble(lg_pstar(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n,l_nk+1)))
         pr_p0_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n) = exp(log_p0_8(1+pil_w:l_ni-pil_e,1+pil_s:l_nj-pil_n))

         pr_p0_dry_8 = 0.0 !NOT COMPLETED GEM-H

      end if

      !Compute dry air mass
      !--------------------
      l_avg_8(1) = 0.0d0
      ll_avg_8(:,:,1) = 0.0d0
      do j=1+pil_s,l_nj-pil_n
      do i=1+pil_w,l_ni-pil_e
         l_avg_8(1) = l_avg_8(1) + pr_p0_dry_8(i,j) * geomh_area_8(i,j) &
                                                    * geomh_mask_8(i,j)
         ll_avg_8(i,j,1) =         pr_p0_dry_8(i,j) * geomh_area_8(i,j) &
                                                    * geomh_mask_8(i,j)
      end do
      end do

      !Compute wet air mass
      !--------------------
      l_avg_8(2) = 0.0d0
      ll_avg_8(:,:,2) = 0.0d0
      do j=1+pil_s,l_nj-pil_n
      do i=1+pil_w,l_ni-pil_e
         l_avg_8(2) = l_avg_8(2) + pr_p0_8(i,j) * geomh_area_8(i,j) &
                                                * geomh_mask_8(i,j)
         ll_avg_8(i,j,2) =         pr_p0_8(i,j) * geomh_area_8(i,j) &
                                                * geomh_mask_8(i,j)
      end do
      end do

      if ( Lctl_rxstat_S == 'GLB_8') then
         call glbsum8 (l_avg_8,ll_avg_8,1,l_ni,1,l_nj,2,          &
                   1+Lam_pil_w, G_ni-Lam_pil_e, 1+Lam_pil_s, G_nj-Lam_pil_n)
         if (Grd_yinyang_L) communicate_S = "GRIDPEERS"
      end if
      call RPN_COMM_allreduce ( l_avg_8,g_avg_ps_8,&
                 2,"MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      g_avg_ps_8(:) = g_avg_ps_8(:) * scale_8

      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ""
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "DRY air mass",time_S,g_avg_ps_8(1),F_comment_S
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "WET air mass",time_S,g_avg_ps_8(2),F_comment_S

      !---------------------------------------------------------------

      return

 1000 format(1X,A12,1X,A7,1X,E19.12,1X,A16)

      end subroutine stat_psadj
