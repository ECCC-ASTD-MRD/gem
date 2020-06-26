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

      subroutine stat_psadj (F_time,F_comment_S)

      use geomh
      use glb_ld
      use gmm_pw
      use HORgrid_options
      use lun
      use psadjust
      use ptopo

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

      integer :: err,i,j
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0_dry_8
      real(kind=REAL64), pointer, dimension(:,:)   :: p0_wet_8
      real(kind=REAL64), pointer, dimension(:,:,:) :: pm_8
      real(kind=REAL64) :: l_avg_8(2),g_avg_8(2),g_avg_ps_wet_8,g_avg_ps_dry_8
      character(len= 9) :: communicate_S
      character(len= 1) :: in_S
      character(len= 7) :: time_S
!
!---------------------------------------------------------------------
!
      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"

      if (F_time==1) in_S = 'P'
      if (F_time==0) in_S = 'M'

      if (F_time==1) time_S = "TIME T1"
      if (F_time==0) time_S = "TIME T0"

      !Obtain Wet surface pressure/Pressure Momentum at appropriate time
      !-----------------------------------------------------------------
      if (F_time==1) p0_wet_8 => pw_p0_plus_8
      if (F_time==1) pm_8     => pw_pm_plus_8
      if (F_time==0) p0_wet_8 => pw_p0_moins_8
      if (F_time==0) pm_8     => pw_pm_moins_8

      !Compute Dry Surface pressure at appropriate time
      !------------------------------------------------
      call dry_sfc_pressure_8 (p0_dry_8,pm_8,p0_wet_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,in_S)

      !Estimate Wet/Dry air mass on CORE
      !---------------------------------
      l_avg_8(1:2) = 0.0d0

      do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e

            l_avg_8(1) = l_avg_8(1) + p0_wet_8(i,j) * geomh_area_mask_8(i,j)
            l_avg_8(2) = l_avg_8(2) + p0_dry_8(i,j) * geomh_area_mask_8(i,j)

         end do
      end do

      call RPN_COMM_allreduce (l_avg_8,g_avg_8,2,"MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      g_avg_ps_wet_8 = g_avg_8(1) * PSADJ_scale_8
      g_avg_ps_dry_8 = g_avg_8(2) * PSADJ_scale_8

      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ""
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "WET air mass",time_S,g_avg_ps_wet_8,F_comment_S
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "DRY air mass",time_S,g_avg_ps_dry_8,F_comment_S
!
!---------------------------------------------------------------------
!
      return

 1000 format(1X,A12,1X,A7,1X,E19.12,1X,A16)

      end subroutine stat_psadj
