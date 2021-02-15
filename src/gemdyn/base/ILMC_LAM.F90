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

!**s/r ILMC_LAM - Ensures monotonicity of interpolated field while preserving mass (Sorenson et al.,2013)

      subroutine ILMC_LAM ( F_ns, F_name_S, F_i0, F_in, F_j0, F_jn )

      use adz_mem
      use adz_options
      use dynkernel_options
      use gem_options
      use gem_timing
      use glb_pil
      use gmm_tracers
      use HORgrid_options
      use ilmc_lam_array
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      character(len=*), intent(in)  :: F_name_S
      integer, intent(in) :: F_ns, F_i0, F_in, F_j0, F_jn

      !object
      !===============================================================================================================
      !     Ensures monotonicity of interpolated field while preserving mass.
      !     Based on Sorenson et al., 2013, Geosci. Model Dev., 6, 1029-1042:
      !     A mass conserving and multi-tracer efficient transport scheme in the online integrated Enviro-HIRLAM model
      !===============================================================================================================

      logical, save :: done_allocate_L=.false.

      integer :: i,j,k,sweep,i_rd,j_rd,k_rd,ii,ix,jx,kx,kx_m,kx_p,j1,j2, &
                 reset(Adz_k0:l_nk,5),w1,w2,size,n,il,ir,jl,jr,offi,offj,ext, &
                 g_i0,g_in,g_j0,g_jn,il_,ir_,jl_,jr_,il_cor,ir_cor,jl_cor,jr_cor,err

      real :: sweep_max(400),sweep_min(400),m_use_max,m_use_min,m_sweep_max,m_sweep_min, &
              o_shoot,u_shoot,ratio_max,ratio_min

      real, pointer, dimension(:,:,:) :: F_ilmc,air_mass_m

      logical :: limit_i_L

      logical, save :: done_print_L=.false.
!
!---------------------------------------------------------------------
!
      call gemtime_start (73, 'ILMC_', 38)

      !Localize Tracer TIME M requiring ILMC monotonicity
      !--------------------------------------------------
      F_ilmc => Adz_stack(F_ns)%dst

      if (Adz_verbose>0) call ilmc_lam_write (1)

      !Recall AIR MASS at TIME_M
      !-------------------------
      air_mass_m => airm0

      if (Adz_verbose>0) call ilmc_lam_write (2)

      !Horizontal grid: MIN/MAX and [F_i0,F_in] x [F_j0,F_jn]
      !------------------------------------------------------
      il = 1-G_halox ; ir = l_ni+G_halox ; jl = 1-G_haloy ; jr = l_nj+G_haloy 

      if (l_west)  il = F_i0
      if (l_east)  ir = F_in
      if (l_south) jl = F_j0
      if (l_north) jr = F_jn

      !Global equivalent to Adz_i0 Adz_in Adz_j0 Adz_jn at LAM boundary
      !----------------------------------------------------------------
      offi = Ptopo_gindx(1,Ptopo_myproc+1)-1
      offj = Ptopo_gindx(3,Ptopo_myproc+1)-1

      ext = 1
      if (Grd_yinyang_L) ext = 2

      g_i0 = 1    + Glb_pil_w - (ext+1) 
      g_in = G_ni - Glb_pil_e +  ext
      g_j0 = 1    + Glb_pil_s - (ext+1)
      g_jn = G_nj - Glb_pil_n +  ext

      !Check if Halo needs to be reduced for PE close to W/E/S/N LAM boundary
      !----------------------------------------------------------------------
      il_ = il ; ir_ = ir ; jl_ = jl ; jr_ = jr

      if (.not.l_west ) il = max(il_,g_i0-offi)
      if (.not.l_east ) ir = min(ir_,g_in-offi)
      if (.not.l_south) jl = max(jl_,g_j0-offj)
      if (.not.l_north) jr = min(jr_,g_jn-offj)

      !Print if Halo needs to be reduced for PE close to W/E/S/N LAM boundary 
      !----------------------------------------------------------------------
      if (.not.done_print_L) then

         call RPN_COMM_allreduce (il -il_,il_cor,1,"MPI_INTEGER","MPI_MAX","GRID",err)
         call RPN_COMM_allreduce (ir_-ir ,ir_cor,1,"MPI_INTEGER","MPI_MAX","GRID",err)
         call RPN_COMM_allreduce (jl -jl_,jl_cor,1,"MPI_INTEGER","MPI_MAX","GRID",err)
         call RPN_COMM_allreduce (jr_-jr ,jr_cor,1,"MPI_INTEGER","MPI_MAX","GRID",err)

         if (Lun_out>0 .and. (il_cor/=0.or.ir_cor/=0.or.jl_cor/=0.or.jr_cor/=0)) then

            if (il_cor/=0) write (Lun_out,*) 'CAUTION: ILMC_RC5 /= ILMC_RC4: Halo reduced for PE close to W by',il_cor
            if (ir_cor/=0) write (Lun_out,*) 'CAUTION: ILMC_RC5 /= ILMC_RC4: Halo reduced for PE close to E by',ir_cor
            if (jl_cor/=0) write (Lun_out,*) 'CAUTION: ILMC_RC5 /= ILMC_RC4: Halo reduced for PE close to S by',jl_cor
            if (jr_cor/=0) write (Lun_out,*) 'CAUTION: ILMC_RC5 /= ILMC_RC4: Halo reduced for PE close to N by',jr_cor

         else if (Lun_out>0) then 

            write (Lun_out,*) 'CAUTION: ILMC_RC5 == ILMC_RC4'

         end if

      end if

      done_print_L = .true.

      !Allocation and Define surrounding cells for each sweep
      !------------------------------------------------------
      if (.NOT.done_allocate_L) then

         call gemtime_start (93, 'ALLOC_', 73)

         !Allocation
         !----------
         allocate (sweep_rd(Adz_ILMC_sweep_max,il:ir,jl:jr,l_nk))

         do k=1,l_nk
            do j=jl,jr
            do i=il,ir
            do n=1,Adz_ILMC_sweep_max

               w1   = 2*n + 1
               w2   = w1-2
               size = 2*(w1**2 + w1*w2 + w2*w2)

               allocate (sweep_rd(n,i,j,k)%i_rd(size))
               allocate (sweep_rd(n,i,j,k)%j_rd(size))
               allocate (sweep_rd(n,i,j,k)%k_rd(size))

            end do
            end do
            end do
         end do

         !Define surrounding cells for each sweep
         !---------------------------------------
         do k=Adz_k0,l_nk

            kx_m = k ; kx_p = k

            do j=jl,jr
            do i=il,ir

               do sweep = 1,Adz_ILMC_sweep_max

                  sweep_rd(sweep,i,j,k)%cell = 0

                  if (.NOT.Schm_autobar_L) then

                     kx_m = max(k-sweep,Adz_k0) ; kx_p = min(k+sweep,l_nk)

                  end if

                  j1 = max(j-sweep,jl) ; j2 = min(j+sweep,jr)

                  do jx = j1,j2

                     do kx = kx_m,kx_p

                        limit_i_L = (jx/=j-sweep.and.jx/=j+sweep)

                        if (.NOT.Schm_autobar_L) limit_i_L = limit_i_L.and.(kx>k-sweep.and.kx<k+sweep)

                        do ix = max(i-sweep,il),min(i+sweep,ir)

                           if (limit_i_L.and.ix/=i-sweep.and.ix/=i+sweep) cycle

                           sweep_rd(sweep,i,j,k)%i_rd(sweep_rd(sweep,i,j,k)%cell + 1) = ix
                           sweep_rd(sweep,i,j,k)%j_rd(sweep_rd(sweep,i,j,k)%cell + 1) = jx
                           sweep_rd(sweep,i,j,k)%k_rd(sweep_rd(sweep,i,j,k)%cell + 1) = kx

                           sweep_rd(sweep,i,j,k)%cell = sweep_rd(sweep,i,j,k)%cell + 1

                        end do

                     end do

                  end do

               end do ! Do sweep

            end do
            end do

         end do

         done_allocate_L = .true.

         call gemtime_stop  (93)

      end if

      call gemtime_start (94, 'X_HALO', 73)

      !Fill Halos of F_ILMC/MIN/MAX (NOTE: Fill Halo of air_mass_m was DONE in set_post_tr)
      !------------------------------------------------------------------------------------
      if (Adz_set_post_tr==1) &
      call rpn_comm_xch_halo (air_mass_m,        l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,G_halox,G_haloy,.false.,.false.,l_ni,0)
      call rpn_comm_xch_halo (F_ilmc,            l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,G_halox,G_haloy,.false.,.false.,l_ni,0)
      call rpn_comm_xch_halo (Adz_post(F_ns)%min,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,G_halox,G_haloy,.false.,.false.,l_ni,0)
      call rpn_comm_xch_halo (Adz_post(F_ns)%max,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk,G_halox,G_haloy,.false.,.false.,l_ni,0)

      Adz_set_post_tr = 0

      call gemtime_stop  (94)

      reset = 0 ! Use to accumulate DIAGNOSTICS in LOOP

      !-----------------------------------------------------------------------------------
      !LOOP: Compute ILMC solution while preserving mass: USE ELEMENTS inside/outside CORE
      !-----------------------------------------------------------------------------------
      call gemtime_start (96, 'LOOP__', 73)
      call ilmc_lam_loop ()
      call gemtime_stop  (96)

      call gemtime_start (98, 'MINMAX', 73)

      !Reset Min-Max Monotonicity if requested
      !---------------------------------------
      if (Adz_ILMC_min_max_L) then

         do k=Adz_k0,l_nk
            do j=F_j0,F_jn
            do i=F_i0,F_in

               if (F_ilmc(i,j,k) < Adz_post(F_ns)%min(i,j,k)) then

                   reset (k,3)   = reset(k,3) + 1
                   F_ilmc(i,j,k) = Adz_post(F_ns)%min(i,j,k)

               end if

               if (F_ilmc(i,j,k) > Adz_post(F_ns)%max(i,j,k)) then

                   reset(k,4)    = reset(k,4) + 1
                   F_ilmc(i,j,k) = Adz_post(F_ns)%max(i,j,k)

               end if

            end do
            end do
         end do

      end if

      call gemtime_stop  (98)

      if (Adz_verbose>0) call ilmc_lam_write (3)

      call gemtime_stop  (73)
!
!---------------------------------------------------------------------
!
      return

contains

      include 'ilmc_lam_loop.inc'
!
!---------------------------------------------------------------------
!
!**s/r ILMC_LAM_write - Write ILMC_LAM diagnostics based on F_numero if verbose is activated

      subroutine ILMC_LAM_write (F_numero)

      use HORgrid_options
      use lun
      use tr3d

      implicit none

      !arguments
      !---------
      integer, intent(in) :: F_numero

      integer :: err,iprod,l_reset(5),g_reset(5)

      real(kind=REAL64), save :: mass_adv_8

      real(kind=REAL64) :: mass_ilmc_8,ratio_8,mass_deficit_8

      logical :: almost_zero

      character(len=9) :: communicate_S
!
!---------------------------------------------------------------------
!
      if (F_numero==1.and.Lun_out>0) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) 'TRACERS: ILMC: Reset Monotonicity without changing Mass of SL advection: ',F_name_S(4:6)

      elseif (F_numero==2) then

         call mass_tr (mass_adv_8,F_ilmc,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0)

         if (Lun_out>0) then

            write(Lun_out,*)    'TRACERS: ILMC: ILMC_min_max_L          =',Adz_ILMC_min_max_L
            write(Lun_out,*)    'TRACERS: ILMC: ILMC_sweep_max          =',Adz_ILMC_sweep_max
            write(Lun_out,1000) 'TRACERS: ILMC: Mass BEFORE ILMC        =',mass_adv_8/Adz_gc_area_8

         end if

      elseif (F_numero==3) then

         !Print Min-Max Monotonicity
         !--------------------------
         communicate_S = "GRID"
         if (Grd_yinyang_L) communicate_S = "MULTIGRID"

         iprod = 1
         if (Grd_yinyang_L) iprod = 2

         l_reset = 0
         g_reset = 0

         do k=Adz_k0,l_nk
            l_reset(1) = reset(k,1) + l_reset(1)
            l_reset(2) = reset(k,2) + l_reset(2)
            l_reset(3) = reset(k,3) + l_reset(3)
            l_reset(4) = reset(k,4) + l_reset(4)
            l_reset(5) = reset(k,5) + l_reset(5)
         end do

         call RPN_COMM_allreduce (l_reset,g_reset,5,"MPI_INTEGER","MPI_SUM",communicate_S,err)

         !Print Masses
         !------------
         call mass_tr (mass_ilmc_8,F_ilmc,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0)

         mass_deficit_8 = mass_ilmc_8 - mass_adv_8

         ratio_8 = 0.
         if (.not.almost_zero(mass_adv_8)) ratio_8 = mass_deficit_8/mass_adv_8*100.

         if (Lun_out>0) then
            write(Lun_out,1000) 'TRACERS: ILMC: Mass    END ILMC        =',mass_ilmc_8/Adz_gc_area_8
            write(Lun_out,*)    'TRACERS: ILMC: # pts OVER/UNDER SHOOT  =',g_reset(5),'over',G_ni*G_nj*l_nk*iprod
            write(Lun_out,*)    'TRACERS: ILMC: # pts RESET MIN many PE =',g_reset(1)
            write(Lun_out,*)    'TRACERS: ILMC: # pts RESET_MAX many PE =',g_reset(2)
            if (Adz_ILMC_min_max_L) then
               write(Lun_out,*)    'TRACERS: ILMC: # pts RESET_MIN_ILMC    =',g_reset(3)
               write(Lun_out,*)    'TRACERS: ILMC: # pts RESET_MAX_ILMC    =',g_reset(4)
            end if
            write(Lun_out,1001) 'TRACERS: ILMC: Rev. Diff. of ',ratio_8
         end if

      end if
!
!----------------------------------------------------------------
!
      return

 1000 format(1X,A40,E20.12)
 1001 format(1X,A29,E11.4,'%')

      end subroutine ILMC_LAM_write

      end subroutine ILMC_LAM
