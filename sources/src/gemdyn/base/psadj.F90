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

!**s/r psadj - Adjust surface pressure for conservation
!
      subroutine psadj (F_kount)
      use cstv
      use dyn_fisl_options
      use init_options
      use gem_options
      use geomh
      use glb_ld
      use glb_pil
      use gmm_geof
      use gmm_itf_mod
      use gmm_vt0
      use HORgrid_options
      use psadjust
      use rstr
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_kount

!author
!     Andre Plante from hzd_main
!
!revision
! v4_05 - Lepine M.         - VMM replacement with GMM
! v4_50 - Qaddouri-PLante   - YY version
! v4_70 - Tanguay M.        - dry air pressure conservation
! v4_80 - Tanguay M.        - REAL*8 with iterations and psadj LAM
!

      type(gmm_metadata) :: mymeta
      integer err,i,j,n,istat,MAX_iteration
      real*8,dimension(l_minx:l_maxx,l_miny:l_maxy,1:l_nk):: pr_m_8,pr_t_8
      real*8,dimension(l_minx:l_maxx,l_miny:l_maxy)       :: pr_p0_0_8,pr_p0_w_0_8
      real*8 l_avg_8(1),g_avg_ps_0_8
      real*8 ll_avg_8(l_ni,l_nj)
      character(len= 9) communicate_S
!
!     ---------------------------------------------------------------
!
      if ( Grd_yinyang_L ) then

         if ( Schm_psadj == 0 ) then
            if ( Schm_psadj_print_L ) then
               call stat_psadj (0,"AFTER DYNSTEP")
            end if

            return
         end if

         if ( Schm_psadj <= 2 .and. Cstv_dt_8*F_kount <= Iau_period ) then
            if (Schm_psadj_print_L) then
               call stat_psadj (0,"AFTER DYNSTEP")
            end if

            return
         end if

      else if ( .not.Grd_yinyang_L .and. Schm_psadj == 0 ) then

         if ( Schm_psadj_print_L ) then
            call stat_psadj (0,"AFTER DYNSTEP")
         end if

         return

      else

         if ( Rstri_rstn_L ) call gem_error(-1,'psadj', &
                        'PSADJ NOT AVAILABLE FOR LAMs in RESTART mode')

         if ( Cstv_dt_8*F_kount <= Iau_period ) then
            if (Schm_psadj_print_L) then
               call stat_psadj (0,"AFTER DYNSTEP")
            end if

            return
         end if

         if (.not.Init_mode_L) call adv_psadj_LAM_0

         if ( Schm_psadj_print_L ) then
            call stat_psadj (0,"AFTER DYNSTEP")
         end if

         return

      end if

      communicate_S = "GRID"
      if (Grd_yinyang_L) communicate_S = "MULTIGRID"

      istat = gmm_get(gmmk_fis0_s,fis0)
      istat = gmm_get(gmmk_st0_s,st0,mymeta)

      MAX_iteration = 3
      if (Schm_psadj==1) MAX_iteration = 1

      do n= 1,MAX_iteration

         !Obtain pressure levels
         !----------------------
         call calc_pressure_8 (pr_m_8,pr_t_8,pr_p0_0_8,st0,l_minx,l_maxx,l_miny,l_maxy,l_nk)

         !Compute dry surface pressure (- Cstv_pref_8)
         !--------------------------------------------
         if (Schm_psadj>=2) then

            call dry_sfc_pressure_8 (pr_p0_w_0_8,pr_m_8,pr_p0_0_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,'M')

         !Compute wet surface pressure (- Cstv_pref_8)
         !--------------------------------------------
         else if (Schm_psadj==1) then

            pr_p0_w_0_8(1:l_ni,1:l_nj) = pr_p0_0_8(1:l_ni,1:l_nj) - Cstv_pref_8

         end if

         l_avg_8 = 0.0d0
         ll_avg_8 = 0.0d0
         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            l_avg_8 = l_avg_8 + pr_p0_w_0_8(i,j) * geomh_area_8(i,j) * geomh_mask_8(i,j)
            ll_avg_8(i,j) =     pr_p0_w_0_8(i,j) * geomh_area_8(i,j) * geomh_mask_8(i,j)
         end do
         end do

         if ( Lctl_rxstat_S == 'GLB_8') then
         call glbsum8 (l_avg_8,ll_avg_8,1,l_ni,1,l_nj,1,          &
                   1+Lam_pil_w, G_ni-Lam_pil_e, 1+Lam_pil_s, G_nj-Lam_pil_n)
         if (Grd_yinyang_L) communicate_S = "GRIDPEERS"
         end if

         call RPN_COMM_allreduce (l_avg_8,g_avg_ps_0_8,1, &
                   "MPI_DOUBLE_PRECISION","MPI_SUM",communicate_S,err)

         g_avg_ps_0_8 = g_avg_ps_0_8 * PSADJ_scale_8

         !Correct surface pressure in order to preserve air mass
         !------------------------------------------------------
         do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            if(fis0(i,j) > 1.) then
               pr_p0_0_8(i,j) = pr_p0_0_8(i,j) + &
                 (PSADJ_g_avg_ps_initial_8 - g_avg_ps_0_8)*PSADJ_fact_8
            else
               pr_p0_0_8(i,j) = pr_p0_0_8(i,j) + &
                 (PSADJ_g_avg_ps_initial_8 - g_avg_ps_0_8)*Cstv_psadj_8
            end if
            st0(i,j) = log(pr_p0_0_8(i,j)/Cstv_pref_8)
         end do
         end do

      end do

      if (Schm_psadj_print_L) then
         call stat_psadj (0,"AFTER DYNSTEP")
      end if
!
!     ---------------------------------------------------------------
!
      return
      end
