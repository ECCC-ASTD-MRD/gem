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

!**s/r Bermejo_Conde - Ensures conservation of interpolated field (Bermejo and Conde,2002)

      subroutine Bermejo_Conde (F_name_S,F_out,F_adv,F_low,F_min,F_max,F_old,F_for_flux_o,F_for_flux_i, &
                                F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0, &
                                F_BC_LAM_Aranami_L,F_BC_min_max_L,F_verbose_L)

      use adz_BC_deficit
      use adz_options
      use ctrl
      use dyn_fisl_options
      use glb_ld
      use gmm_itf_mod
      use gmm_tracers
      use HORgrid_options
      use lun
      use rstr
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      character (len=*),intent(in) :: F_name_S                        !Name of field to be ajusted
      integer,          intent(in) :: F_minx,F_maxx,F_miny,F_maxy     !Dimension H
      integer,          intent(in) :: F_nk                            !Number of vertical levels
      integer,          intent(in) :: F_i0,F_in,F_j0,F_jn,F_k0        !Scope of operator
      logical,          intent(in) :: F_BC_LAM_Aranami_L              !LAM: T=FLUX Aranami F=FLUX Zerroukat-Shipway (ZLF)
      logical,          intent(in) :: F_BC_min_max_L                  !Clipping Min Max + Proportional MF after BC IF T
      logical,          intent(in) :: F_verbose_L                     !VERBOSE IF T
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(out) :: F_out !Corrected (conservative) solution
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_adv !High-order SL solution
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_low !Low-order SL solution
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_min !MIN over cell
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_max !MAX over cell
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_old !Field at previous time step
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_for_flux_o !Advected mixing ratio with 0 in NEST
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk),intent(in)  :: F_for_flux_i !Advected mixing ratio with 0 in CORE

      !object
      !==============================================================================================================
      !     Ensures conservation of interpolated field.
      !     Based on Bermejo and Conde, 2002, MWR, 130, 423-430: A conservative quasi-monotone semi-Lagrangian scheme
      !     Based on Bermejo/Conde Template from GEM3 (Gravel and DeGrandpre 2012)
      !==============================================================================================================

      !------------------------------------------------------------------------------------------

      integer :: i,j,k,err,count(F_k0:F_nk,3),l_count(3),g_count(3),iprod,istat,tracer_id

      real(kind=REAL64) :: mass_old_8,mass_tot_old_8,mass_adv_8,mass_tot_adv_8, &
                mass_out_8,mass_tot_out_8,mass_wei_8,mass_bc_w_8,mass_bc_l_8,proportion_8, &
                mass_deficit_8,lambda_8,correction_8,p_exp_8,H_minus_L_8,ratio_8, &
                mass_flux_o_8,mass_flux_i_8,mass_bflux_8

      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: weight,bc_w,bc_l

      real :: bc_wrk

      logical :: LAM_L,almost_zero,reset_epsilon_L

      character(len= 9) :: communicate_S

      real(kind=REAL64), parameter :: ONE_8=1.d0

      real, pointer, dimension(:,:,:) :: air_mass_p,air_mass_m

      logical, save :: done_BC_deficit_init_L = .false.

      !-----------------------------------------------------------------------------------------

      if (Adz_BC_pexp_n>PEXP_LIST_MAX) call handle_error(-1,'BERMEJO-CONDE','BC PEXP not valid')

      if (Adz_BC_weight>3) call handle_error(-1,'BERMEJO-CONDE','BC WEIGHT not valid')

      if (Schm_psadj==0.and.Grd_yinyang_L.and..NOT.Ctrl_testcases_L) &
         call handle_error(-1,'BERMEJO-CONDE','Schm_psadj should be > 0 when YIN-YANG')

      if (MAXTR3D_/=MAXTR3D) call handle_error(-1,'BERMEJO-CONDE','MAXTR3D_/=MAXTR3D')

      !Assign a unique id to this tracer, so we can look up the corresponding mass deficit
      !-----------------------------------------------------------------------------------
      if (.not.Rstri_rstn_L) then

         if (.not.done_BC_deficit_init_L) then

            KEEP_mass_deficit_8(1:MAXTR3D_) = 0.0d0

            tracer_name(1:MAXTR3D_) = ''

            done_BC_deficit_init_L = .true.

         end if

         do tracer_id = 1,MAXTR3D
            if (tracer_name(tracer_id) == F_name_S(4:7)) exit
         end do
         if (tracer_id >= MAXTR3D) then
            do tracer_id = 1,MAXTR3D
               if (tracer_name(tracer_id) == '') exit
            end do
            tracer_name(tracer_id) = F_name_s(4:7)
         end if

      else

         do tracer_id = 1,MAXTR3D
            if (tracer_name(tracer_id) == F_name_S(4:7)) exit
         end do

      end if

      LAM_L = .not.Grd_yinyang_L

      !--------
      !PRINTING
      !--------
      if (F_verbose_L.and.Lun_out>0) then

                                                write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
                                                write(Lun_out,*) 'TRACERS: Bermejo-Conde: Restore Mass of High-order solution: ',F_name_S(4:6)
         if (LAM_L.and.     F_BC_LAM_Aranami_L) write(Lun_out,*) 'TRACERS: Bermejo-Conde: LAM: FLUX based on Aranami et al. (2015)'
         if (LAM_L.and..not.F_BC_LAM_Aranami_L) write(Lun_out,*) 'TRACERS: Bermejo-Conde: LAM: FLUX based on Zerroukat-Shipway (2017)'
                                                write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'

      end if

      !Reset AIR MASS at TIME_P and TIME_M
      !-----------------------------------
      istat = gmm_get(gmmk_airm1_s,airm1)
      istat = gmm_get(gmmk_airm0_s,airm0)

      air_mass_p => airm1
      air_mass_m => airm0

      if (     Adz_BC_LEGACY_L) p_exp_8 = 1.0
      if (.not.Adz_BC_LEGACY_L) p_exp_8 = Adz_BC_pexp_list(Adz_BC_pexp_n)

      !Default values if no Mass correction
      !------------------------------------
      F_out(F_i0:F_in,F_j0:F_jn,F_k0:F_nk) = F_adv(F_i0:F_in,F_j0:F_jn,F_k0:F_nk)

      call mass_tr (mass_old_8,F_old,air_mass_p,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)
      call mass_tr (mass_adv_8,F_adv,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)

      mass_bflux_8 = 0.0d0

      !Estimate mass of FLUX_out and mass of FLUX_in
      !---------------------------------------------
      if (LAM_L.and.F_BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then

         call mass_tr (mass_flux_o_8,F_for_flux_o,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,1,l_ni,1,l_nj,F_k0)
         call mass_tr (mass_flux_i_8,F_for_flux_i,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,1,l_ni,1,l_nj,F_k0)

         mass_bflux_8 = mass_flux_i_8 - mass_flux_o_8

      end if

      reset_epsilon_L = .TRUE.

      !Recall previous mass_deficit (EPSILON)
      !--------------------------------------
      if (reset_epsilon_L)  mass_old_8 = mass_old_8 - KEEP_mass_deficit_8(tracer_id)

      mass_tot_old_8 = mass_old_8 + mass_bflux_8
      mass_tot_adv_8 = mass_adv_8

      mass_deficit_8 = mass_tot_adv_8 - mass_tot_old_8

      !--------
      !PRINTING
      !--------
      if (F_verbose_L) then

         ratio_8 = 0.0d0
         if (.not.almost_zero(mass_tot_old_8)) ratio_8 = mass_deficit_8/mass_tot_old_8*100.

         if (Lun_out>0) then
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: P_exponent              =',p_exp_8
             if (Adz_BC_weight==1) &
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: BC_WEIGHT: ADDITIVE APPROACH'
             if (Adz_BC_weight==2) &
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: BC_WEIGHT: MULTIPLICATIVE APPROACH'
             if (Adz_BC_weight==3) &
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: BC_WEIGHT: ADDITIVE APPROACH + Factor pr_k/pr_s'
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: BC_min_max_L            =',F_BC_min_max_L
             write(Lun_out,1000) 'TRACERS: Bermejo-Conde: Mass BEFORE B.-C.       =',mass_tot_adv_8/Adz_gc_area_8
             write(Lun_out,1000) 'TRACERS: Bermejo-Conde: Mass to RESTORE         =',mass_tot_old_8/Adz_gc_area_8
             write(Lun_out,1001) 'TRACERS: Bermejo-Conde: Ori. Diff. of ',ratio_8
         end if

      end if

      weight = 0.0

      !---------------------------------------------------
      !Prepare Bermejo-Conde correction: Additive approach
      !---------------------------------------------------
      if (Adz_BC_weight==1) then

!$omp parallel private(k,i,j,H_minus_L_8) &
!$omp shared(weight,mass_deficit_8)
!$omp do
         do k=F_k0,F_nk

            do j=F_j0,F_jn
            do i=F_i0,F_in

               H_minus_L_8 = F_adv(i,j,k) - F_low(i,j,k)

               if (int(sign(ONE_8,mass_deficit_8)) == int(sign(ONE_8,H_minus_L_8))) then
                  weight(i,j,k) = abs(H_minus_L_8)**p_exp_8*sign(ONE_8,mass_deficit_8)
               end if

            end do
            end do

         end do
!$omp    end do
!$omp end parallel

      !---------------------------------------------------------
      !Prepare Bermejo-Conde correction: Multiplicative approach
      !---------------------------------------------------------
      else if (Adz_BC_weight==2) then

!$omp parallel private(k,i,j,H_minus_L_8) &
!$omp shared(weight,mass_deficit_8)
!$omp do
         do k=F_k0,F_nk

            do j=F_j0,F_jn
            do i=F_i0,F_in

               H_minus_L_8 = F_adv(i,j,k) - F_low(i,j,k)

               if (int(sign(ONE_8,mass_deficit_8)) == int(sign(ONE_8,H_minus_L_8))) then
                  weight(i,j,k) = abs(H_minus_L_8)**p_exp_8*sign(ONE_8,mass_deficit_8)*max(0.0,F_adv(i,j,k))
               end if

            end do
            end do

         end do
!$omp    end do
!$omp end parallel

      !--------------------------------------------------------------------------
      !Prepare Bermejo-Conde correction:  Additive approach with Factor pr_k/pr_s
      !--------------------------------------------------------------------------
      elseif (Adz_BC_weight==3) then

        istat = gmm_get(gmmk_pkps_s,pkps)

!$omp parallel private(k,i,j,H_minus_L_8) &
!$omp shared(weight,mass_deficit_8,pkps)
!$omp do
         do k=F_k0,F_nk

            do j=F_j0,F_jn
            do i=F_i0,F_in

               H_minus_L_8 = F_adv(i,j,k) - F_low(i,j,k)

               if (int(sign(ONE_8,mass_deficit_8)) == int(sign(ONE_8,H_minus_L_8))) then
                  weight(i,j,k) = abs(H_minus_L_8)**p_exp_8*sign(ONE_8,mass_deficit_8)*pkps(i,j,k)
               end if

            end do
            end do

         end do
!$omp    end do
!$omp end parallel

      end if

      call mass_tr (mass_wei_8,weight,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)

      if ( almost_zero(mass_wei_8) ) then

         if (F_verbose_L.and.Lun_out>0) &
         write(Lun_out,1002) 'TRACERS: Bermejo-Conde: Diff. too small             =', &
                             mass_tot_adv_8/Adz_gc_area_8,mass_tot_old_8/Adz_gc_area_8,(mass_tot_adv_8-mass_tot_old_8)/Adz_gc_area_8

         return

      end if

      lambda_8 = mass_deficit_8/mass_wei_8

      if (F_verbose_L.and.Lun_out>0) write(Lun_out,1003) 'TRACERS: Bermejo-Conde: LAMBDA                  = ',lambda_8

      !-----------------------------------------------------------------------------
      !Apply Bermejo-Conde correction without accumulating statistics (BC_MIN_MAX=F)
      !-----------------------------------------------------------------------------
      if (.NOT.F_verbose_L.and..NOT.F_BC_min_max_L) then

!$omp parallel do private(k,i,j,correction_8) &
!$omp shared(weight,lambda_8)
         do k=F_k0,F_nk

            do j=F_j0,F_jn
            do i=F_i0,F_in

               correction_8 = lambda_8 * weight(i,j,k)

               F_out(i,j,k) = F_adv(i,j,k) - correction_8

            end do
            end do

         end do
!$omp end parallel do

      !-----------------------------------------------------------------------------
      !Apply Bermejo-Conde correction without accumulating statistics (BC_MIN_MAX=T)
      !-----------------------------------------------------------------------------
      else if (.NOT.F_verbose_L.and.F_BC_min_max_L) then

!$omp parallel do private(k,i,j,correction_8,bc_wrk) &
!$omp shared(bc_l,weight,lambda_8)
         do k=F_k0,F_nk

            do j=F_j0,F_jn
            do i=F_i0,F_in

               correction_8 = lambda_8 * weight(i,j,k)

               bc_wrk      = F_adv(i,j,k) - correction_8

               bc_l(i,j,k) = bc_wrk

               !Reset Min-Max Monotonicity
               !--------------------------
               if (correction_8 > 0.d0 .and. bc_wrk < F_min(i,j,k)) then
                   bc_l(i,j,k) = F_min(i,j,k)
               end if
               if (correction_8 < 0.d0 .and. bc_wrk > F_max(i,j,k)) then
                   bc_l(i,j,k) = F_max(i,j,k)
               end if

            end do
            end do

         end do
!$omp end parallel do

         call mass_tr (mass_bc_l_8,bc_l,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)

      !---------------------------------------------------------------------------
      !Apply Bermejo-Conde correction while accumulating statistics (BC_MIN_MAX=F)
      !---------------------------------------------------------------------------
      else if (F_verbose_L.and..NOT.F_BC_min_max_L) then

!$omp parallel do private(k,i,j,correction_8) &
!$omp shared(weight,count,lambda_8)
         do k=F_k0,F_nk

            count(k,1) = 0.
            count(k,2) = 0.
            count(k,3) = 0.

            do j=F_j0,F_jn
            do i=F_i0,F_in

               correction_8 = lambda_8 * weight(i,j,k)

               if (.not.almost_zero(correction_8)) count(k,3) = count(k,3) + 1

               F_out(i,j,k) = F_adv(i,j,k) - correction_8

            end do
            end do

         end do
!$omp end parallel do

      !---------------------------------------------------------------------------
      !Apply Bermejo-Conde correction while accumulating statistics (BC_MIN_MAX=T)
      !---------------------------------------------------------------------------
      else if (F_verbose_L.and.F_BC_min_max_L) then

!$omp parallel do private(k,i,j,correction_8) &
!$omp shared(bc_w,bc_l,weight,count,lambda_8)
         do k=F_k0,F_nk

            count(k,1) = 0.
            count(k,2) = 0.
            count(k,3) = 0.

            do j=F_j0,F_jn
            do i=F_i0,F_in

               correction_8 = lambda_8 * weight(i,j,k)

               if (.not.almost_zero(correction_8)) count(k,3) = count(k,3) + 1

               bc_w(i,j,k) = F_adv(i,j,k) - correction_8

               bc_l(i,j,k) = bc_w(i,j,k)

               !Reset Min-Max Monotonicity
               !--------------------------
               if (correction_8 > 0.d0 .and. bc_w(i,j,k) < F_min(i,j,k)) then
                   count(k,1)  = count(k,1) + 1
                   bc_l(i,j,k) = F_min(i,j,k)
               end if
               if (correction_8 < 0.d0 .and. bc_w(i,j,k) > F_max(i,j,k)) then
                   count(k,2)  = count(k,2) + 1
                   bc_l(i,j,k) = F_max(i,j,k)
               end if

            end do
            end do

         end do
!$omp end parallel do

         call mass_tr (mass_bc_w_8,bc_w,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)

         call mass_tr (mass_bc_l_8,bc_l,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)

      end if

      !------------------------------------------------------------------
      !Apply Proportional mass-fixer due to change in mass (BC_MIN_MAX=T)
      !------------------------------------------------------------------
      if (F_BC_min_max_L) then

         proportion_8 = mass_tot_old_8/mass_bc_l_8

!$omp parallel do private(k,i,j) &
!$omp shared(bc_l,proportion_8)
         do k=F_k0,F_nk

            do j=F_j0,F_jn
            do i=F_i0,F_in

               F_out(i,j,k) = proportion_8 * bc_l(i,j,k)

            end do
            end do

         end do
!$omp end parallel do

      end if

      call mass_tr (mass_out_8,F_out,air_mass_m,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_i0,F_in,F_j0,F_jn,F_k0)

      mass_tot_out_8 = mass_out_8

      mass_deficit_8 = mass_tot_out_8 - mass_tot_old_8

      !Store current mass_deficit (EPSILON) for next time step
      !-------------------------------------------------------
      if (reset_epsilon_L) KEEP_mass_deficit_8(tracer_id) = mass_deficit_8

      !--------
      !PRINTING
      !--------
      if (F_verbose_L) then

         !Print Min-Max Monotonicity
         !--------------------------
         communicate_S = "GRID"
         if (Grd_yinyang_L) communicate_S = "MULTIGRID"

         iprod = 1
         if (Grd_yinyang_L) iprod = 2

         l_count = 0
         g_count = 0

         do k=F_k0,F_nk
            l_count(1) = count(k,1) + l_count(1)
            l_count(2) = count(k,2) + l_count(2)
            l_count(3) = count(k,3) + l_count(3)
         end do

         call rpn_comm_Allreduce (l_count,g_count,3,"MPI_INTEGER","MPI_SUM",communicate_S,err)

         !Print Masses
         !------------
         ratio_8 = 0.0d0
         if (.not.almost_zero(mass_tot_old_8)) ratio_8 = mass_deficit_8/mass_tot_old_8*100.

         if (Lun_out>0) then

             if (F_BC_min_max_L) then
             write(Lun_out,1000) 'TRACERS: Bermejo-Conde: Mass BEFORE MINMAX      =' ,mass_bc_w_8/Adz_gc_area_8
             write(Lun_out,1000) 'TRACERS: Bermejo-Conde: Mass  AFTER MINMAX      =' ,mass_bc_l_8/Adz_gc_area_8
             write(Lun_out,1000) 'TRACERS: Bermejo-Conde: Mass  AFTER PROPORT. MF =' ,mass_tot_out_8/Adz_gc_area_8
             end if
             write(Lun_out,1000) 'TRACERS: Bermejo-Conde: Mass    END B.-C.       =' ,mass_tot_out_8/Adz_gc_area_8
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: # pts treated by B.-C.  =' ,g_count(3),'over',G_ni*G_nj*F_nk*iprod
             if (F_BC_min_max_L) then
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: # pts RESET_MIN_BC      =' ,g_count(1)
             write(Lun_out,*)    'TRACERS: Bermejo-Conde: # pts RESET_MAX_BC      =' ,g_count(2)
             end if
             write(Lun_out,1001) 'TRACERS: Bermejo-Conde: Rev. Diff. of ',ratio_8

             if (LAM_L.and.F_BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then
             write(Lun_out,1004) 'TRACERS: Bermejo-Conde: STATS: N=',mass_out_8/Adz_gc_area_8,' O=',mass_old_8/Adz_gc_area_8, &
                                 ' FI=',mass_flux_i_8/Adz_gc_area_8,' FO=',mass_flux_o_8/Adz_gc_area_8
             else
             write(Lun_out,1005) 'TRACERS: Bermejo-Conde: STATS: N=',mass_out_8/Adz_gc_area_8,' O=',mass_old_8/Adz_gc_area_8
             end if

         end if

      end if

 1000 format(1X,A49,E20.12)
 1001 format(1X,A38,E11.4,'%')
 1002 format(1X,A49,3(E20.12,1X))
 1003 format(1X,A49,F20.2)
 1004 format(1X,A33,E19.12,A3,E19.12,A3,E19.12,A3,E19.12)
 1005 format(1X,A33,E19.12,A3,E19.12)

      return
      end
