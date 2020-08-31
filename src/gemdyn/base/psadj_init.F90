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

!**s/r psadj_init - Estimate area for Yin-Yang/LAM and Store air mass at initial time for Yin-Yang

      subroutine psadj_init (F_kount)

      use cstv
      use dynkernel_options
      use dyn_fisl_options
      use init_options
      use geomh
      use glb_ld
      use gmm_geof
      use gmm_pw
      use HORgrid_options
      use lam_options
      use psadjust
      use rstr

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_kount

      !object
      !===================================================================================
      !     Estimate area for Yin-Yang/LAM and Store air mass at initial time for Yin-Yang
      !===================================================================================

      integer :: err,i,j,i0_c,in_c,j0_c,jn_c
      logical, save :: done_area_L       = .false.
      logical, save :: done_yy_initial_L = .false.
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0_dry_8,p0_1_8
      real(kind=REAL64), pointer, dimension(:,:)   :: p0_wet_8
      real(kind=REAL64), pointer, dimension(:,:,:) :: pm_8
      real(kind=REAL64) :: l_avg_8,x_avg_8
      character(len= 9) :: communicate_S
      logical :: almost_zero
!
!     ---------------------------------------------------------------
!
      if ( Schm_psadj < 0 .and. Schm_psadj > 2) &
         call gem_error(-1,'PSADJ_INIT','PSADJ OPTION NOT VALIDE')

      if ( Schm_psadj == 2 .and. trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') &
         call gem_error(-1,'PSADJ_INIT','PSADJ=2 NOT AVAILABLE in GEM-H')

      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"

      !Estimate area
      !-------------
      if (.not.done_area_L) then

         l_avg_8 = 0.0d0
         x_avg_8 = 0.0d0

         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e

               l_avg_8 = l_avg_8 + geomh_area_mask_8(i,j)

               if (fis0(i,j) > 1.) x_avg_8 = x_avg_8 + geomh_area_mask_8(i,j)

            end do
         end do

         call RPN_COMM_allreduce (l_avg_8, PSADJ_scale_8, 1, &
                  "MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

         call RPN_COMM_allreduce (x_avg_8, PSADJ_fact_8,  1, &
                  "MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

         PSADJ_fact_8  = PSADJ_fact_8/PSADJ_scale_8
         PSADJ_scale_8 = 1.0d0/PSADJ_scale_8

         if (.not.almost_zero(PSADJ_fact_8)) PSADJ_fact_8 = (1.0d0-(1.0d0-PSADJ_fact_8)*Cstv_psadj_8)/PSADJ_fact_8

         done_area_L = .true.

      end if

      if (Schm_psadj==0) goto 999

      if (.not.Grd_yinyang_L) goto 999

      if (Cstv_dt_8*F_kount <= Iau_period) goto 999

      if (done_yy_initial_L) goto 999

      if (Rstri_rstn_L) goto 999

      done_yy_initial_L = .true.

      !Obtain Wet surface pressure/Pressure Momentum at TIME P
      !-------------------------------------------------------
      p0_wet_8 => pw_p0_plus_8
      pm_8     => pw_pm_plus_8

      !Compute Dry surface pressure at TIME P (Schm_psadj==2)
      !------------------------------------------------------
      if (Schm_psadj==2) call dry_sfc_pressure_8 (p0_dry_8,pm_8,p0_wet_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,'P')

      i0_c = 1+pil_w ; j0_c = 1+pil_s ; in_c = l_ni-pil_e ; jn_c = l_nj-pil_n

      !Obtain Surface pressure minus Cstv_pref_8
      !-----------------------------------------
      if (Schm_psadj==1) p0_1_8(i0_c:in_c,j0_c:jn_c) = p0_wet_8(i0_c:in_c,j0_c:jn_c) - Cstv_pref_8
      if (Schm_psadj==2) p0_1_8(i0_c:in_c,j0_c:jn_c) = p0_dry_8(i0_c:in_c,j0_c:jn_c) - Cstv_pref_8

      !Store air mass at initial time
      !------------------------------
      l_avg_8 = 0.0d0

      do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e

            l_avg_8 = l_avg_8 + p0_1_8(i,j) * geomh_area_mask_8(i,j)

         end do
      end do

      call RPN_COMM_allreduce (l_avg_8,PSADJ_g_avg_ps_initial_8, 1, &
               "MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

      PSADJ_g_avg_ps_initial_8 = PSADJ_g_avg_ps_initial_8 * PSADJ_scale_8

  999 if (Schm_psadj_print_L) call stat_psadj (1,"BEFORE DYNSTEP")
!
!     ---------------------------------------------------------------
!
      return
      end
