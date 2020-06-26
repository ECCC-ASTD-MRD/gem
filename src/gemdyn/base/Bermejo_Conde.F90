!(n)---------------------------------- LICENCE BEGIN -------------------------------
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

      subroutine Bermejo_Conde ( F_bc, F_ntr_bc, F_i0, F_in, F_j0, F_jn )

      use adz_mem
      use adz_options
      use ctrl
      use dyn_fisl_options
      use gem_timing
      use gmm_tracers
      use HORgrid_options

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,  intent(in)                         :: F_i0,F_in,F_j0,F_jn,F_ntr_bc
      type(bc), intent(inout), dimension(F_ntr_bc) :: F_bc !Pointers for tracers using B-C

      !object
      !==============================================================================================================
      !     Ensures conservation of interpolated field.
      !     Based on Bermejo and Conde, 2002, MWR, 130, 423-430: A conservative quasi-monotone semi-Lagrangian scheme
      !     Based on Bermejo/Conde Template from GEM3 (Gravel and DeGrandpre 2012)
      !==============================================================================================================

      integer :: i,j,k,n

      real(kind=REAL64), dimension(F_ntr_bc) :: mass_p_8,mass_m_8,mass_fo_8,mass_fi_8, &
                                                mass_tot_p_8,mass_tot_m_8,mass_wei_8,  &
                                                mass_bc_8,mass_bc_r_8,mass_deficit_8

      integer,           dimension(F_ntr_bc) :: cycle_tr

      real(kind=REAL64) :: mass_tot_bc_8,mass_bflux_8,proportion_8,lambda_8,corr_8,p_exp_8,H_minus_L_8

      logical :: LAM_L,almost_zero,reset_epsilon_L,BC_LAM_Aranami_L

      real(kind=REAL64), parameter :: ONE_8=1.d0

      real, pointer, dimension(:,:,:) :: air_mass_p,air_mass_m

      real :: bc_wrk
!
!---------------------------------------------------------------------
!
      if (Schm_psadj==0.and.Grd_yinyang_L.and..NOT.Ctrl_testcases_L) &
         call gem_error(-1,'BERMEJO-CONDE','Schm_psadj should be > 0 when YIN-YANG')

      call gemtime_start (74, 'BermejoC', 38)

      LAM_L = .not.Grd_yinyang_L

      BC_LAM_Aranami_L = LAM_L.and.Adz_BC_LAM_flux==1

      !Recall AIR MASS at TIME_P/M
      !---------------------------
      air_mass_p => airm1
      air_mass_m => airm0

      !Evaluate Mass of Tracer TIME P/M and Mass of Flux_out/Flux_in for all tracers using Bermejo-Conde
      !-------------------------------------------------------------------------------------------------
      call mass_tr_PM ( mass_p_8, mass_m_8, mass_fo_8, mass_fi_8, F_bc, airm1, airm0, &
                        l_minx, l_maxx, l_miny, l_maxy, l_nk, F_i0, F_in, F_j0, F_jn, Adz_k0, F_ntr_bc )

      reset_epsilon_L = .TRUE.

      cycle_tr = 0

      do n=1,F_ntr_bc

         if (F_bc(n)%pexp_n>PEXP_LIST_MAX) call gem_error(-1,'BERMEJO-CONDE','BC PEXP   not valid')
         if (F_bc(n)%weight>3)             call gem_error(-1,'BERMEJO-CONDE','BC WEIGHT not valid')

         if (     F_bc(n)%LEGACY_L) p_exp_8 = 1.0
         if (.not.F_bc(n)%LEGACY_L) p_exp_8 = Adz_BC_pexp_list(F_bc(n)%pexp_n)

         mass_bflux_8 = mass_fi_8(n) - mass_fo_8(n)

         !Recall previous mass_deficit (EPSILON)
         !--------------------------------------
         if (reset_epsilon_L)  mass_p_8(n) = mass_p_8(n) - F_bc(n)%mass_deficit

         mass_tot_p_8(n) = mass_p_8(n) + mass_bflux_8
         mass_tot_m_8(n) = mass_m_8(n)

         mass_deficit_8(n) = mass_tot_m_8(n) - mass_tot_p_8(n)

         call gemtime_start (16, 'WEIGHT', 74)

         !---------------------------------------------------
         !Prepare Bermejo-Conde correction: Additive approach
         !---------------------------------------------------
         if (F_bc(n)%weight==1) then

            do k=Adz_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     F_bc(n)%w(i,j,k) = 0.

                     H_minus_L_8 = F_bc(n)%m(i,j,k) - F_bc(n)%lin(i,j,k)

                     if (int(sign(ONE_8,mass_deficit_8(n))) == int(sign(ONE_8,H_minus_L_8))) then
                        F_bc(n)%w(i,j,k) = abs(H_minus_L_8)**p_exp_8*sign(ONE_8,mass_deficit_8(n))
                     end if

                  end do
               end do
            end do

         !---------------------------------------------------------
         !Prepare Bermejo-Conde correction: Multiplicative approach
         !---------------------------------------------------------
         else if (F_bc(n)%weight==2) then

            do k=Adz_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     F_bc(n)%w(i,j,k) = 0.

                     H_minus_L_8 = F_bc(n)%m(i,j,k) - F_bc(n)%lin(i,j,k)

                     if (int(sign(ONE_8,mass_deficit_8(n))) == int(sign(ONE_8,H_minus_L_8))) then
                        F_bc(n)%w(i,j,k) = abs(H_minus_L_8)**p_exp_8*sign(ONE_8,mass_deficit_8(n))*max(0.0,F_bc(n)%m(i,j,k))
                     end if

                  end do
               end do
            end do

         !-------------------------------------------------------------------------
         !Prepare Bermejo-Conde correction: Additive approach with Factor pr_k/pr_s
         !-------------------------------------------------------------------------
         elseif (F_bc(n)%weight==3) then

            do k=Adz_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     F_bc(n)%w(i,j,k) = 0.

                     H_minus_L_8 = F_bc(n)%m(i,j,k) - F_bc(n)%lin(i,j,k)

                     if (int(sign(ONE_8,mass_deficit_8(n))) == int(sign(ONE_8,H_minus_L_8))) then
                        F_bc(n)%w(i,j,k) = abs(H_minus_L_8)**p_exp_8*sign(ONE_8,mass_deficit_8(n))*pkps(i,j,k)
                     end if

                  end do
               end do
            end do

         end if

         call gemtime_stop  (16)

      end do

      !Evaluate Mass of Tracer for all tracers using Bermejo-Conde (Component %w)
      !--------------------------------------------------------------------------
      call mass_tr_X (mass_wei_8,F_bc,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0,F_ntr_bc,1)

      do n=1,F_ntr_bc

         if (Adz_verbose>0) call Bermejo_Conde_write (1)
         if (Adz_verbose>0) call Bermejo_Conde_write (2)

         !-----------------------------------
         !Difference too small: Nothing to do
         !-----------------------------------
         if ( almost_zero(mass_wei_8(n)) ) then

            if (Adz_verbose>0) call Bermejo_Conde_write (3)

            cycle_tr(n) = 1

            cycle

         end if

         lambda_8 = mass_deficit_8(n)/mass_wei_8(n)

         if (Adz_verbose>0) call Bermejo_Conde_write (4)

         call gemtime_start (17, 'MINMAX', 74)

         !---------------------------------------------
         !Apply Bermejo-Conde correction (BC_MIN_MAX=F)
         !---------------------------------------------
         if (.NOT.Adz_BC_min_max_L) then

            do k=Adz_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     corr_8 = lambda_8 * F_bc(n)%w(i,j,k)

                     F_bc(n)%m(i,j,k) = F_bc(n)%m(i,j,k) - corr_8

                  end do
               end do
            end do

         !---------------------------------------------------
         !Apply Bermejo-Conde correction (BC_MIN_MAX=T) PART1
         !---------------------------------------------------
         else if (Adz_BC_min_max_L) then

            do k=Adz_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     corr_8 = lambda_8 * F_bc(n)%w(i,j,k)

                     bc_wrk = F_bc(n)%m(i,j,k) - corr_8

                     F_bc(n)%m(i,j,k) = bc_wrk

                     !Reset Min-Max Monotonicity
                     !--------------------------
                     if (corr_8 > 0.d0 .and. bc_wrk < F_bc(n)%min(i,j,k)) then
                        F_bc(n)%m(i,j,k) = F_bc(n)%min(i,j,k)
                     end if
                     if (corr_8 < 0.d0 .and. bc_wrk > F_bc(n)%max(i,j,k)) then
                        F_bc(n)%m(i,j,k) = F_bc(n)%max(i,j,k)
                     end if

                  end do
               end do
            end do

         end if

         call gemtime_stop  (17)

      end do

      call gemtime_start (17, 'MINMAX', 74)

      !---------------------------------------------------
      !Apply Bermejo-Conde correction (BC_MIN_MAX=T) PART2
      !---------------------------------------------------
      if (Adz_BC_min_max_L) then

         !Evaluate Mass of Tracer for all tracers using Bermejo-Conde (Component %m)
         !--------------------------------------------------------------------------
         call mass_tr_X (mass_bc_r_8,F_bc,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0,F_ntr_bc,2)

         !Apply Proportional mass-fixer due to change in mass
         !---------------------------------------------------
         do n=1,F_ntr_bc

            if (cycle_tr(n) == 1) cycle

            proportion_8 = mass_tot_p_8(n)/mass_bc_r_8(n)

            do k=Adz_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     F_bc(n)%m(i,j,k) = proportion_8 * F_bc(n)%m(i,j,k)

                  end do
               end do
            end do

         end do

      end if

      call gemtime_stop  (17)

      if (reset_epsilon_L.or.Adz_verbose>0) then

         !Evaluate Mass of Tracer for all tracers using Bermejo-Conde (Component %m)
         !--------------------------------------------------------------------------
         call mass_tr_X (mass_bc_8,F_bc,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0,F_ntr_bc,2)

         do n=1,F_ntr_bc

            if (cycle_tr(n) == 1) cycle

            mass_tot_bc_8  = mass_bc_8(n)

            mass_deficit_8(n) = mass_tot_bc_8 - mass_tot_p_8(n)

            !Store current mass_deficit (EPSILON) for next time step
            !-------------------------------------------------------
            if (reset_epsilon_L) F_bc(n)%mass_deficit = mass_deficit_8(n)

            if (Adz_verbose>0) call Bermejo_Conde_write (5)

         end do

      end if

      call gemtime_stop (74)
!
!---------------------------------------------------------------------
!
      return

contains

!
!---------------------------------------------------------------------
!

!**s/r Bermejo_Conde_write - Write Bermejo_Conde's diagnostics based on F_numero if verbose is activated

      subroutine Bermejo_Conde_write (F_numero)

      use lun
      use tr3d

      implicit none

      !arguments
      !---------
      integer, intent(in) :: F_numero

      integer :: err,count(Adz_k0:l_nk,3),iprod

      integer :: l_count(3,F_ntr_bc),g_count(3,F_ntr_bc)

      real(kind=REAL64), dimension(F_ntr_bc) :: mass_bc_w_8

      real(kind=REAL64) :: ratio_8

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: bc_w

      character(len=9) :: communicate_S
!
!---------------------------------------------------------------------
!
      if (F_numero==1.and.Lun_out>0) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Restore Mass of High-order solution  '

         if (LAM_L.and.     BC_LAM_Aranami_L) &
         write(Lun_out,*) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: LAM: FLUX based on Aranami et al. (2015)'

         if (LAM_L.and..not.BC_LAM_Aranami_L) &
         write(Lun_out,*) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: LAM: FLUX based on Zerroukat-Shipway (2017)'

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'

      elseif (F_numero==2) then

         ratio_8 = 0.0d0
         if (.not.almost_zero(mass_tot_p_8(n))) ratio_8 = mass_deficit_8(n)/mass_tot_p_8(n)*100.

         if (Lun_out>0) then
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: P_exponent              =',p_exp_8
             if (F_bc(n)%weight==1) &
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: BC_WEIGHT: ADDITIVE APPROACH'
             if (F_bc(n)%weight==2) &
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: BC_WEIGHT: MULTIPLICATIVE APPROACH'
             if (F_bc(n)%weight==3) &
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: BC_WEIGHT: ADDITIVE APPROACH + Factor pr_k/pr_s'
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: BC_min_max_L            =',Adz_BC_min_max_L
             write(Lun_out,1000) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Mass BEFORE BC          =',&
                                 mass_tot_m_8(n)/Adz_gc_area_8
             write(Lun_out,1000) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Mass to RESTORE         =',&
                                 mass_tot_p_8(n)/Adz_gc_area_8
             write(Lun_out,1001) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Ori. Diff. of ',ratio_8
         end if

      elseif (F_numero==3) then

         if (Lun_out>0) &
         write(Lun_out,1002) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Diff. too small             =', &
                             mass_tot_m_8(n)/Adz_gc_area_8,mass_tot_p_8(n)/Adz_gc_area_8,&
                             (mass_tot_m_8(n)-mass_tot_p_8(n))/Adz_gc_area_8

      elseif (F_numero==4) then

         if (Lun_out>0) write(Lun_out,1003) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: LAMBDA                  =',lambda_8

         l_count(1:3,n) = 0

         if (.NOT.Adz_BC_min_max_L) then

            do k=Adz_k0,l_nk

               count(k,3) = 0.

               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     corr_8 = lambda_8 * F_bc(n)%w(i,j,k)

                     if (.not.almost_zero(corr_8)) count(k,3) = count(k,3) + 1

                  end do
               end do

               l_count(3,n) = count(k,3) + l_count(3,n)

            end do

         else if (Adz_BC_min_max_L) then

            do k=Adz_k0,l_nk

               count(k,1) = 0.
               count(k,2) = 0.
               count(k,3) = 0.

               do j=F_j0,F_jn
                  do i=F_i0,F_in

                     corr_8 = lambda_8 * F_bc(n)%w(i,j,k)

                     if (.not.almost_zero(corr_8)) count(k,3) = count(k,3) + 1

                     bc_w(i,j,k) = F_bc(n)%m(i,j,k) - corr_8

                     !Reset Min-Max Monotonicity
                     !--------------------------
                     if (corr_8 > 0.d0 .and. bc_w(i,j,k) < F_bc(n)%min(i,j,k)) then
                        count(k,1)  = count(k,1) + 1
                     end if
                     if (corr_8 < 0.d0 .and. bc_w(i,j,k) > F_bc(n)%max(i,j,k)) then
                        count(k,2)  = count(k,2) + 1
                     end if

                  end do
               end do

               l_count(1,n) = count(k,1) + l_count(1,n)
               l_count(2,n) = count(k,2) + l_count(2,n)
               l_count(3,n) = count(k,3) + l_count(3,n)

            end do

            call mass_tr (mass_bc_w_8(n),bc_w,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0)

         end if

      elseif (F_numero==5) then

         call mass_tr (mass_bc_8(n),F_bc(n)%m,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0)

         mass_tot_bc_8  = mass_bc_8(n)

         mass_deficit_8(n) = mass_tot_bc_8 - mass_tot_p_8(n)

         !Print Min-Max Monotonicity
         !--------------------------
         communicate_S = "GRID"
         if (Grd_yinyang_L) communicate_S = "MULTIGRID"

         iprod = 1
         if (Grd_yinyang_L) iprod = 2

         g_count(1:3,n) = 0

         call rpn_comm_Allreduce (l_count(1,n),g_count(1,n),3,"MPI_INTEGER","MPI_SUM",communicate_S,err)

         !Print Masses
         !------------
         ratio_8 = 0.0d0
         if (.not.almost_zero(mass_tot_p_8(n))) ratio_8 = mass_deficit_8(n)/mass_tot_p_8(n)*100.

         if (Lun_out>0) then

             write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'

             if (Adz_BC_min_max_L) then
             write(Lun_out,1000) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Mass BEFORE MINMAX      =' ,&
                                 mass_bc_w_8(n)/Adz_gc_area_8
             write(Lun_out,1000) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Mass  AFTER MINMAX      =' ,&
                                 mass_bc_r_8(n)/Adz_gc_area_8
             write(Lun_out,1000) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Mass  AFTER PROPORT. MF =' ,&
                                 mass_tot_bc_8 /Adz_gc_area_8
             end if
             write(Lun_out,1000) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Mass    END BC          =' ,&
                                 mass_tot_bc_8/Adz_gc_area_8
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: # pts treated by BC ='     ,g_count(3,n),&
                                 'over',G_ni*G_nj*l_nk*iprod
             if (Adz_BC_min_max_L) then
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: # pts RESET_MIN_BC      =' ,g_count(1,n)
             write(Lun_out,*)    'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: # pts RESET_MAX_BC      =' ,g_count(2,n)
             end if
             write(Lun_out,1001) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: Rev. Diff. of ',ratio_8

             if (LAM_L.and.BC_LAM_Aranami_L.and..not.Ctrl_theoc_L) then
             write(Lun_out,1004) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: STATS: N=',mass_bc_8(n)/Adz_gc_area_8,' O=',&
                                 mass_p_8(n)/Adz_gc_area_8, &
                                 ' FI=',mass_fi_8(n)/Adz_gc_area_8,' FO=',mass_fo_8(n)/Adz_gc_area_8
             else
             write(Lun_out,1005) 'TRACERS: ',F_bc(n)%name_S(4:6),' Bermejo-Conde: STATS: N=',mass_bc_8(n)/Adz_gc_area_8,' O=',&
                                 mass_p_8(n)/Adz_gc_area_8
             end if

         end if

      end if

      return

 1000 format(1X,A9,A3,A41,E20.12)
 1001 format(1X,A9,A3,A30,E20.12,'%')
 1002 format(1X,A9,A3,A45,3(E20.12,1X))
 1003 format(1X,A9,A3,A41,F20.2)
 1004 format(1X,A9,A3,A25,E19.12,A3,E19.12,A3,E19.12,A3,E19.12)
 1005 format(1X,A9,A3,A25,E19.12,A3,E19.12)

      end subroutine Bermejo_Conde_write

      end subroutine Bermejo_Conde

