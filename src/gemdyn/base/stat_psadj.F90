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

      use adz_mem
      use adz_options 
      use gem_options
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

      include 'mpif.h'
      include 'rpn_comm.inc'
      integer :: err,i,j,i0_sb,in_sb,j0_sb,jn_sb,comm
      logical :: do_subset_GY_L
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0_dry_8
      real(kind=REAL64), pointer, dimension(:,:)   :: p0_wet_8
      real(kind=REAL64), pointer, dimension(:,:,:) :: pm_8
      real(kind=REAL64) :: l_avg_8(2),g_avg_8(2),g_avg_ps_wet_8,g_avg_ps_dry_8
      real(kind=REAL64) :: gathV(2,Ptopo_numproc*Ptopo_ncolors)
      character(len= 9) :: communicate_S
      character(len= 1) :: in_S
      character(len= 7) :: time_S
!
!---------------------------------------------------------------------
!
      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"

      comm = RPN_COMM_comm ('MULTIGRID')

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

      if (Legacy_reduce_L) then

         call RPN_COMM_allreduce (l_avg_8,g_avg_8,2,"MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      else

         call MPI_Allgather(l_avg_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
         do i=1,2
            g_avg_8(i) = sum(gathV(i,:))
         end do

      end if

      g_avg_ps_wet_8 = g_avg_8(1) * PSADJ_scale_8
      g_avg_ps_dry_8 = g_avg_8(2) * PSADJ_scale_8

      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ""
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "C:WET air mass",time_S,g_avg_ps_wet_8,F_comment_S
      if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1000) "C:DRY air mass",time_S,g_avg_ps_dry_8,F_comment_S

      do_subset_GY_L = Grd_yinyang_L.and.adz_pil_sub_s_g /= -1

      !Estimate Wet/Dry air mass on SUBSET of YIN
      !------------------------------------------
      if (do_subset_GY_L) then

         i0_sb = 1+Adz_pil_sub_w ; j0_sb = 1+Adz_pil_sub_s ; in_sb = l_ni-Adz_pil_sub_e; jn_sb = l_nj-Adz_pil_sub_n

         l_avg_8(1:2) = 0.0d0

         do j=j0_sb,jn_sb
            do i=i0_sb,in_sb

               l_avg_8(1) = l_avg_8(1) + p0_wet_8(i,j) * geomh_area_mask_8(i,j)
               l_avg_8(2) = l_avg_8(2) + p0_dry_8(i,j) * geomh_area_mask_8(i,j)

            end do
         end do

         if (Ptopo_couleur/=0) l_avg_8(1:2) = 0.0d0 

         if (Legacy_reduce_L) then

            call RPN_COMM_allreduce (l_avg_8,g_avg_8,2,"MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

         else

            call MPI_Allgather(l_avg_8,2,MPI_DOUBLE_PRECISION,gathV,2,MPI_DOUBLE_PRECISION,comm,err)
            do i=1,2
               g_avg_8(i) = sum(gathV(i,:))
            end do

         end if

         g_avg_ps_wet_8 = g_avg_8(1) / Adz_gs_area_8 
         g_avg_ps_dry_8 = g_avg_8(2) / Adz_gs_area_8 

         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,*) ""
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1001) "SB:WET air mass",time_S,g_avg_ps_wet_8,F_comment_S
         if (Lun_out>0.and.Ptopo_couleur==0) write(Lun_out,1001) "SB:DRY air mass",time_S,g_avg_ps_dry_8,F_comment_S

      end if
!
!---------------------------------------------------------------------
!
      return

 1000 format(1X,A14,1X,A7,1X,E19.12,1X,A16)
 1001 format(1X,A15,1X,A7,1X,E19.12,1X,A16)

      end subroutine stat_psadj
