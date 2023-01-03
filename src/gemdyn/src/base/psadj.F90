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

!**s/r psadj - Adjust surface pressure for conservation in Yin-Yang/LAM

      subroutine psadj (F_kount)

      use adz_mem
      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use omp_timing
      use geomh
      use gmm_geof
      use gmm_vt0
      use gmm_pw
      use HORgrid_options
      use init_options
      use mem_tracers
      use metric
      use psadjust
      use ptopo
      use rstr
      use tdpack, only : rgasd_8
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_kount

      !object
      !=============================================================
      !     Adjust surface pressure for conservation in Yin-Yang/LAM
      !=============================================================
      include 'mpif.h'
      include 'rpn_comm.inc'
      integer :: err,i,j,k,n,MAX_iteration,empty_i,n_reduce,i0_c,in_c,j0_c,jn_c,comm
      real :: empty
      real(kind=REAL64),dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0_dry_1_8,p0_dry_0_8
      real(kind=REAL64),dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0_1_8,p0_0_8,fl_0_8
      real(kind=REAL64), pointer, dimension(:,:)   :: p0_wet_1_8,p0_wet_0_8
      real(kind=REAL64), pointer, dimension(:,:,:) :: pm_8
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: qts,delps,pw_pm
!      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: sumq
      real(kind=REAL64) :: l_avg_8(2),g_avg_8(2),g_avg_ps_1_8,g_avg_ps_0_8,g_avg_fl_0_8,substract_8
      real(kind=REAL64) :: gathV(2,Ptopo_numproc*Ptopo_ncolors),gathS(Ptopo_numproc*Ptopo_ncolors)
      logical :: LAM_L
      real, pointer, dimension(:,:,:) :: tr
!     
!---------------------------------------------------------------------
!
      !Update Pressure PW_MOINS in accordance with TIME M
      !--------------------------------------------------
      call pressure ( pw_pm_moins,pw_pt_moins,pw_p0_moins,pw_log_pm,pw_log_pt, &
                      pw_pm_moins_8,pw_p0_moins_8, &
                      l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )

      LAM_L = .not.Grd_yinyang_L

      if ( Schm_psadj_print_L ) call stat_psadj (0,"AFTER DYNSTEP")

      if ( Schm_psadj == 0 ) return

      if ( Cstv_dt_8*F_kount <= Iau_period ) return

      if ( LAM_L .and. Init_mode_L ) return

      comm = RPN_COMM_comm ('MULTIGRID')

      substract_8 = 1.0d0
      if (LAM_L) substract_8 = 0.0d0

      !LAM: Estimate FLUX_out/FLUX_in and Evaluate water tracers at TIME M
      !-------------------------------------------------------------------
      if (LAM_L) then

         Adz_flux => Adz_flux_3CWP_PS

         !Estimate FLUX_out/FLUX_in using Tracer=1 based on Aranami et al. (2015)
         !-----------------------------------------------------------------------
         call adz_BC_LAM_Aranami (empty,Adz_pb,Adz_num_b,1,Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,empty_i,MAXTR3D+1)

         !Evaluate water tracers at TIME M (Schm_psadj==2)
         !------------------------------------------------
         if (Schm_psadj==2) then

!$omp parallel
            call sumhydro_hlt (sumq_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,Tr3d_ntr,trt0, Schm_wload_L)
!$omp end parallel
            tr=>tracers_M(Tr3d_hu)%pntr

            do k=Adz_k0t,l_nk
               sumq_8(1:l_ni,1:l_nj,k) = sumq_8(1:l_ni,1:l_nj,k) + tracers_M(Tr3d_hu)%pntr(1:l_ni,1:l_nj,k)
            end do

         else

            sumq_8 = 0.

         end if
!stop
      end if

      !Estimate or Recuperate air mass on CORE at TIME P
      !-------------------------------------------------
      if (LAM_L.or.Schm_psadj==2) then

         !Obtain Wet surface pressure/Pressure Momentum at TIME P
         !-------------------------------------------------------
         p0_wet_1_8 => pw_p0_plus_8
         pm_8       => pw_pm_plus_8

         !Compute Dry surface pressure at TIME P (Schm_psadj==2)
         !------------------------------------------------------
         if (Schm_psadj==2) call dry_sfc_pressure_8 (p0_dry_1_8,pm_8,p0_wet_1_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,Adz_k0t,'P')

         i0_c = 1+pil_w ; j0_c = 1+pil_s ; in_c = l_ni-pil_e ; jn_c = l_nj-pil_n

         !Obtain Surface pressure minus Cstv_pref_8
         !-----------------------------------------
         if (Schm_psadj==1.and.Adz_k0t==1) p0_1_8(i0_c:in_c,j0_c:jn_c) = p0_wet_1_8(i0_c:in_c,j0_c:jn_c)  - substract_8*Cstv_pref_8 
         if (Schm_psadj==1.and.Adz_k0t/=1) p0_1_8(i0_c:in_c,j0_c:jn_c) = pm_8(i0_c:in_c,j0_c:jn_c,l_nk+1) - pm_8(i0_c:in_c,j0_c:jn_c,Adz_k0t) - substract_8*Cstv_pref_8 
         if (Schm_psadj==2)                p0_1_8(i0_c:in_c,j0_c:jn_c) = p0_dry_1_8(i0_c:in_c,j0_c:jn_c)  - substract_8*Cstv_pref_8

         !Estimate air mass on CORE at TIME P
         !-----------------------------------
         l_avg_8(1) = 0.0d0

         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               l_avg_8(1) = l_avg_8(1) + p0_1_8(i,j) * geomh_area_mask_8(i,j)
            end do
         end do

         call MPI_Allgather(l_avg_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)
         g_avg_8(1) = sum(gathS)

         g_avg_ps_1_8 = g_avg_8(1) * PSADJ_scale_8

      else

         g_avg_ps_1_8 = PSADJ_g_avg_ps_initial_8

      end if

      MAX_iteration = 10 
      if (Schm_psadj==1) MAX_iteration = 1

      do n = 1,MAX_iteration

         if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H'.and.n==1) qts(:,:) = qt0(:,:,l_nk+1)

         !Obtain Wet surface pressure/Pressure Momentum at TIME M
         !-------------------------------------------------------
         p0_wet_0_8 => pw_p0_moins_8
         pm_8       => pw_pm_moins_8

         !Compute Dry surface pressure at TIME M (Schm_psadj==2)
         !------------------------------------------------------
         if (Schm_psadj==2) call dry_sfc_pressure_8 (p0_dry_0_8,pm_8,p0_wet_0_8,l_minx,l_maxx,l_miny,l_maxy,l_nk,Adz_k0t,'M')

         !!! The scope of validity of p0_dry_0_8 is 1+pil_s to
         !l_nj-pil_n in dry_sfc_pressure_8(), so we can't compute over
         !the full 1:l_nj here.  And anyways, only the interior
         !(non-pilot) points are used later.
         i0_c = 1+pil_w ; j0_c = 1+pil_s ; in_c = l_ni-pil_e ; jn_c = l_nj-pil_n

         !Obtain Surface pressure minus Cstv_pref_8
         !-----------------------------------------
         if (Schm_psadj==1.and.Adz_k0t==1) p0_0_8(i0_c:in_c,j0_c:jn_c) = p0_wet_0_8(i0_c:in_c,j0_c:jn_c)  - substract_8*Cstv_pref_8 
         if (Schm_psadj==1.and.Adz_k0t/=1) p0_0_8(i0_c:in_c,j0_c:jn_c) = pm_8(i0_c:in_c,j0_c:jn_c,l_nk+1) - pm_8(i0_c:in_c,j0_c:jn_c,Adz_k0t) - substract_8*Cstv_pref_8 
         if (Schm_psadj==2)                p0_0_8(i0_c:in_c,j0_c:jn_c) = p0_dry_0_8(i0_c:in_c,j0_c:jn_c)  - substract_8*Cstv_pref_8

         !Obtain FLUX on NEST+CORE at TIME M
         !----------------------------------
         if (LAM_L) then

            do j=Adz_j0b,Adz_jnb
               fl_0_8(:,j) = 0.0d0
               do k=1,l_nk
                  do i=Adz_i0b,Adz_inb
                     fl_0_8(i,j) = fl_0_8(i,j) + (1.-sumq_8(i,j,k)) * (pm_8(i,j,k+1) - pm_8(i,j,k)) &
                                   * (Adz_flux(1)%fi(i,j,k) - Adz_flux(1)%fo(i,j,k))
                 end do
               end do
            end do

         end if

         l_avg_8(1:2) = 0.0d0
         g_avg_8(1:2) = 0.0d0

         n_reduce = 1
         if (LAM_L) n_reduce = 2

         !Estimate air mass on CORE at TIME M
         !-----------------------------------
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e

               l_avg_8(1) = l_avg_8(1) + p0_0_8(i,j) * geomh_area_mask_8(i,j)

            end do
         end do

         !Estimate FLUX mass on NEST+CORE at TIME M
         !-----------------------------------------
         if (LAM_L) then

            do j=Adz_j0b,Adz_jnb
               do i=Adz_i0b,Adz_inb
                  l_avg_8(2) = l_avg_8(2) + fl_0_8(i,j) * geomh_area_mask_8(i,j)
               end do
            end do
            
            call MPI_Allgather(l_avg_8,n_reduce,MPI_DOUBLE_PRECISION,gathV,n_reduce,MPI_DOUBLE_PRECISION,comm,err)
            do i=1,n_reduce
               g_avg_8(i) = sum(gathV(i,:))
            end do

         else

            call MPI_Allgather(l_avg_8,n_reduce,MPI_DOUBLE_PRECISION,gathS,n_reduce,MPI_DOUBLE_PRECISION,comm,err)

            g_avg_8(1) = sum(gathS)

         end if

         g_avg_ps_0_8 = g_avg_8(1) * PSADJ_scale_8
         g_avg_fl_0_8 = g_avg_8(2) * PSADJ_scale_8

         !Correct surface pressure in order to preserve air mass on CORE (taking FLUX mass into account in LAM)
         !-----------------------------------------------------------------------------------------------------
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e

               if (fis0(i,j) > 1.) then
                  p0_wet_0_8(i,j) = p0_wet_0_8(i,j) + &
                  ((g_avg_ps_1_8 - g_avg_ps_0_8) + g_avg_fl_0_8)*PSADJ_fact_8
               else
                  p0_wet_0_8(i,j) = p0_wet_0_8(i,j) + &
                  ((g_avg_ps_1_8 - g_avg_ps_0_8) + g_avg_fl_0_8)*Cstv_psadj_8
               end if

               if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

                  st0(i,j) = log(p0_wet_0_8(i,j)/Cstv_pref_8)

               else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

                  qt0(i,j,l_nk+1) = rgasd_8*Cstv_Tstr_8 * &
                                    (log(p0_wet_0_8(i,j))-GVM%lg_pstar_8(i,j,l_nk+1))

               end if

            end do
         end do

         !GEM-H: Adjust q at the other vertical levels at TIME M
         !------------------------------------------------------
         if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e

                  delps(i,j) = exp(GVM%lg_pstar_8(i,j,l_nk+1))*&
                              (exp(qt0(i,j,l_nk+1)/(rgasd_8*Cstv_Tstr_8))-&
                               exp(qts(i,j)/(rgasd_8*Cstv_Tstr_8)))
               end do
            end do

            do k=l_nk,1,-1

               do j=1+pil_s,l_nj-pil_n
                  do i=1+pil_w,l_ni-pil_e

                     pw_pm(i,j) = exp(qt0(i,j,k)/(rgasd_8*Cstv_Tstr_8)+GVM%lg_pstar_8(i,j,k))

                     qt0(i,j,k) = rgasd_8*Cstv_Tstr_8*&
                                  (log(pw_pm(i,j)+delps(i,j)*exp(GVM%lg_pstar_8(i,j,k))/&
                                  exp(GVM%lg_pstar_8(i,j,l_nk+1)))-GVM%lg_pstar_8(i,j,k))

                  end do
               end do

            end do

         end if

         !Update Pressure PW_MOINS in accordance with TIME M
         !--------------------------------------------------
         call pressure ( pw_pm_moins,pw_pt_moins,pw_p0_moins,pw_log_pm,pw_log_pt, &
                         pw_pm_moins_8,pw_p0_moins_8, &
                         l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )

      end do !END MAX_iteration

      if ( Schm_psadj_print_L.and.Lun_out>0.and.Ptopo_couleur==0 ) then

         write(Lun_out,*)    ''
         if (Schm_psadj==1) write(Lun_out,*) 'PSADJ: Conservation WET air'
         if (Schm_psadj==2) write(Lun_out,*) 'PSADJ: Conservation DRY air'

         write(Lun_out,*)    ''
         write(Lun_out,*)    '------------------------------------------------------------------------------'
         write(Lun_out,1004) 'PSADJ STATS: N=',g_avg_ps_0_8,' O=',g_avg_ps_1_8,' F=',g_avg_fl_0_8
         write(Lun_out,*)    '------------------------------------------------------------------------------'

      end if

      if ( Schm_psadj_print_L ) call stat_psadj (0,"AFTER PSADJ")
!
!---------------------------------------------------------------------
!
      return

 1004 format(1X,A15,E19.12,A3,E19.12,A3,E19.12)

      end
