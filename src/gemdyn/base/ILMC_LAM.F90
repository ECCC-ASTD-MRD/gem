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

      subroutine ILMC_LAM ( F_ns, F_itr, F_i0, F_in, F_j0, F_jn )

      use adz_mem
      use adz_options
      use array_ilmc
      use dynkernel_options
      use gem_options
      use gmm_itf_mod
      use gmm_tracers
      use HORgrid_options
      use lun
      use tr3d

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(IN) :: F_ns, F_itr, F_i0, F_in, F_j0, F_jn

      !object
      !===============================================================================================================
      !     Ensures monotonicity of interpolated field while preserving mass.
      !     Based on Sorenson et al., 2013, Geosci. Model Dev., 6, 1029-1042:
      !     A mass conserving and multi-tracer efficient transport scheme in the online integrated Enviro-HIRLAM model
      !===============================================================================================================

      logical, save :: done_allocate_L=.false.

      integer :: SIDE_max_0,i,j,k,err,sweep,i_rd,j_rd,k_rd,ii,ix,jx,kx,kx_m,kx_p,iprod, &
                 shift,j1,j2,reset(Adz_k0:l_nk,3),l_reset(3),g_reset(3),w1,w2,size,n, &
                 il,ir,jl,jr,il_c,ir_c,jl_c,jr_c,istat

      real :: sweep_max(400),sweep_min(400),m_use_max,m_use_min,m_sweep_max,m_sweep_min, &
              o_shoot,u_shoot,ratio_max,ratio_min

      real(kind=REAL64) :: mass_adv_8,mass_ilmc_8,ratio_8,mass_deficit_8

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: high,copy,F_min,F_max

      logical :: limit_i_L,almost_zero

      character(len=9) :: communicate_S

      real, pointer, dimension(:,:,:) :: air_mass_m

      character(len=8) :: name_S

      real, pointer, dimension(:,:,:) :: F_ilmc
!
!---------------------------------------------------------------------
!
      name_S = 'TR/'//trim(Tr3d_name_S(F_itr))//':M'

      F_ilmc => Adz_stack(F_ns)%dst

      F_min(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) = Adz_post(l_minx:l_maxx,l_miny:l_maxy,1:l_nk,(F_ns-1)*3+2) !USE POINTER
      F_max(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) = Adz_post(l_minx:l_maxx,l_miny:l_maxy,1:l_nk,(F_ns-1)*3+3) !USE POINTER

      if (Adz_verbose>0.and.Lun_out>0) then

         write(Lun_out,*) 'TRACERS: ----------------------------------------------------------------------'
         write(Lun_out,*) 'TRACERS: ILMC: Reset Monotonicity without changing Mass of SL advection: ',name_S(4:6)

      end if

      !Reset AIR MASS at TIME_M
      !------------------------
      istat = gmm_get(gmmk_airm0_s,airm0)

      air_mass_m => airm0

      high(:,:,:) = F_ilmc(:,:,:)

      call mass_tr (mass_adv_8,high,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0)

      if (Adz_verbose>0.and.Lun_out>0) then

         write(Lun_out,*)    'TRACERS: ILMC: ILMC_min_max_L          =',Adz_ILMC_min_max_L
         write(Lun_out,*)    'TRACERS: ILMC: ILMC_sweep_max          =',Adz_ILMC_sweep_max
         write(Lun_out,1000) 'TRACERS: ILMC: Mass BEFORE ILMC        =',mass_adv_8/Adz_gc_area_8

      end if

      !CORE
      !----
      il_c = F_i0
      ir_c = F_in
      jl_c = F_j0
      jr_c = F_jn

      !Horizontal grid: MIN/MAX and CORE limits
      !----------------------------------------
      il = l_minx
      ir = l_maxx
      jl = l_miny
      jr = l_maxy
      if (l_west)  il = il_c
      if (l_east)  ir = ir_c
      if (l_south) jl = jl_c
      if (l_north) jr = jr_c

      !Allocation and Define surrounding cells for each sweep
      !------------------------------------------------------
      if (.NOT.done_allocate_L) then

         allocate (sweep_rd(Adz_ILMC_sweep_max,l_ni,l_nj,l_nk))

         do k=1,l_nk
            do j=jl_c,jr_c
            do i=il_c,ir_c
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

         do k=Adz_k0,l_nk

            kx_m = k
            kx_p = k

            do j=jl_c,jr_c
            do i=il_c,ir_c

               do sweep = 1,Adz_ILMC_sweep_max

                  sweep_rd(sweep,i,j,k)%cell = 0

                  if (.NOT.Schm_autobar_L) then
                     kx_m = max(k-sweep,Adz_k0)
                     kx_p = min(k+sweep,l_nk)
                  end if

                  j1 = max(j-sweep,jl)
                  j2 = min(j+sweep,jr)

                  do jx = j1,j2

                     do kx = kx_m,kx_p


                        limit_i_L = (jx/=j-sweep.and.jx/=j+sweep) .and. &
                                    ((.NOT.Schm_autobar_L.and.(kx>k-sweep.and.k<k+sweep)).or.Schm_autobar_L)

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

      end if

      !Fill Halos of MIN/MAX (NOTE: Fill Halo of air_mass_m DONE in adv_tracers)
      !-------------------------------------------------------------------------
      call rpn_comm_xch_halo (F_min,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      call rpn_comm_xch_halo (F_max,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !Evaluate admissible distance between threads
      !--------------------------------------------
      SIDE_max_0 = (2*Adz_ILMC_sweep_max + 1)*2+1

      reset = 0

      !-----------------------------------------------------------------------------------
      !LOOP1: Compute ILMC solution F_ilmc while preserving mass: USE ELEMENTS INSIDE CORE
      !-----------------------------------------------------------------------------------
      call ilmc_lam_loop_core_y ()

      !Fill Halo of HIGH
      !-----------------
      call rpn_comm_xch_halo (high,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !Set NESTING boundary of HIGH to ZERO since Adjoint is done after LOOP2
      !----------------------------------------------------------------------
      if (l_west)  high(l_minx:il-1,l_miny:l_maxy,Adz_k0:l_nk) = 0.
      if (l_east)  high(ir+1:l_maxx,l_miny:l_maxy,Adz_k0:l_nk) = 0.
      if (l_south) high(l_minx:l_maxx,l_miny:jl-1,Adz_k0:l_nk) = 0.
      if (l_north) high(l_minx:l_maxx,jr+1:l_maxy,Adz_k0:l_nk) = 0.

      copy = high

      !-------------------------------------------------------------------------------------------------
      !LOOP2: Compute ILMC solution F_ilmc while preserving mass: USE ELEMENTS INSIDE HALO (internal PE)
      !-------------------------------------------------------------------------------------------------
      call ilmc_lam_loop_core_n ()

      !Recover perturbation HIGH in HALO
      !---------------------------------
      do k=Adz_k0,l_nk
         do j=jl,jr
            do i=il,ir
               if ((i>=il_c.and.i<=ir_c).and.(j>=jl_c.and.j<=jr_c)) cycle
               high(i,j,k) = high(i,j,k) - copy(i,j,k)
            end do
         end do
      end do

      copy = high

      !Adjoint of Fill Halo of HIGH
      !----------------------------
      call rpn_comm_adj_halo (copy,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,l_nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !Initialize F_ilmc while maintaining NESTING values
      !--------------------------------------------------
      F_ilmc(il_c:ir_c,jl_c:jr_c,Adz_k0:l_nk) = copy(il_c:ir_c,jl_c:jr_c,Adz_k0:l_nk)

      !Reset Min-Max Monotonicity if requested
      !---------------------------------------
      if (Adz_ILMC_min_max_L) then

         do k=Adz_k0,l_nk
            do j=jl_c,jr_c
            do i=il_c,ir_c

               if (F_ilmc(i,j,k) < F_min(i,j,k)) then

                   reset (k,1)   = reset(k,1) + 1
                   F_ilmc(i,j,k) = F_min(i,j,k)

               end if

               if (F_ilmc(i,j,k) > F_max(i,j,k)) then

                   reset(k,2)    = reset(k,2) + 1
                   F_ilmc(i,j,k) = F_max(i,j,k)

               end if

            end do
            end do
         end do

      end if

      if (Adz_verbose>0) then

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
         end do

         call RPN_COMM_allreduce (l_reset,g_reset,3,"MPI_INTEGER","MPI_SUM",communicate_S,err)

         !Print Masses
         !------------
         call mass_tr (mass_ilmc_8,F_ilmc,air_mass_m,l_minx,l_maxx,l_miny,l_maxy,l_nk,F_i0,F_in,F_j0,F_jn,Adz_k0)

         mass_deficit_8 = mass_ilmc_8 - mass_adv_8

         ratio_8 = 0.
         if (.not.almost_zero(mass_adv_8)) ratio_8 = mass_deficit_8/mass_adv_8*100.

         if (Lun_out>0) then
         write(Lun_out,1000) 'TRACERS: ILMC: Mass    END ILMC        =',mass_ilmc_8/Adz_gc_area_8
         write(Lun_out,*)    'TRACERS: ILMC: # pts OVER/UNDER SHOOT  =',g_reset(3),'over',G_ni*G_nj*l_nk*iprod
         if (Adz_ILMC_min_max_L) then
         write(Lun_out,*)    'TRACERS: ILMC: # pts RESET_MIN_ILMC    =',g_reset(1)
         write(Lun_out,*)    'TRACERS: ILMC: # pts RESET_MAX_ILMC    =',g_reset(2)
         end if
         write(Lun_out,1001) 'TRACERS: ILMC: Rev. Diff. of ',ratio_8
         end if

      end if

 1000 format(1X,A40,E20.12)
 1001 format(1X,A29,E11.4,'%')
!
!---------------------------------------------------------------------
!
      return

contains

      include 'ilmc_lam_loop_core_y.inc'
      include 'ilmc_lam_loop_core_n.inc'

      end
