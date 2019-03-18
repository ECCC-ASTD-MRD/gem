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
!
      subroutine adv_cubic (F_name, fld_out,  fld_in, F_capx, F_capy, F_capz, &
                            F_ni, F_nj, F_nk, F_minx, F_maxx, F_miny, F_maxy, &
                            F_nind, F_ii, F_i0, F_in, F_j0, F_jn, F_k0,       &
                            F_lev_S, F_mono_kind, F_mass_kind)

      use adv_grid
      use adv_options
      use glb_ld
      use gmm_itf_mod
      use gmm_tracers
      use HORgrid_options
      implicit none
#include <arch_specific.hf>
      character(len=1) :: F_lev_S
      character(len=*), intent(in) :: F_name
      integer, intent(in) :: F_ni,F_nj,F_nk ! dims of position fields
      integer, intent(in) :: F_minx,F_maxx,F_miny, F_maxy ! wind fields array bounds
      integer :: F_k0             !I, vertical scope k0 to F_nk
      integer :: F_i0, F_j0, F_in, F_jn
      integer ,  intent(in) :: F_nind
      integer ,  dimension(F_nind*4) :: F_ii    ! precomputed indices to be used in tricubic lagrangian interp
      real, intent(in)::  F_capx(*), F_capy(*), F_capz(*) !I, upstream positions at t1
      integer, intent(in) :: F_mono_kind !I, Kind of Shape-Preserving
      integer, intent(in) :: F_mass_kind !I, Kind of Mass conservation
!     @objective  prepare for  cubic interpolation of RHS and Tracers
!     @author RPN-A Model Infrastructure Group (based on adx_interp_gmm , adx_int_rhs , adx_interp ) June 2015
!     @arguments

      logical :: mono_L, BC_LAM_Aranami_L, BC_LAM_zlf_L
      integer :: err, k, nbpts, i0_2, in_2, j0_2, jn_2, i0_0, in_0, j0_0, jn_0
      real, dimension(adv_lminx:adv_lmaxx,adv_lminy:adv_lmaxy,F_nk) :: fld_adv,adv_rho,adv_o,adv_i
      real, dimension(F_ni,F_nj,F_nk) :: wrkc,w_mono,w_lin,w_min,w_max,w_cub_o,w_cub_i
      real, dimension(F_minx:F_maxx, F_miny:F_maxy ,F_nk), intent(inout)  :: Fld_in
      real, dimension(F_minx:F_maxx, F_miny:F_maxy ,F_nk), intent(out) :: Fld_out
      real, dimension(F_minx:F_maxx, F_miny:F_maxy ,F_nk) :: cub_o,cub_i,store_pilot
      real, dimension(1,1,1), target :: empty
!---------------------------------------------------------------------
!
      mono_L = .false.

      if (F_name(1:3) == 'TR/') mono_L = (F_mono_kind==1)

      if ( adv_rhst_mono_L .and. (F_name=='RHST_S')) mono_L=.true.

      !Check if Mass Conservation/Shape-Preserving is activated for current tracer
      !---------------------------------------------------------------------------
      Adv_Mass_Cons_tr_L = F_mono_kind>1.or.F_mass_kind>0

      !Initializations Bermejo-Conde LAM
      !---------------------------------
      BC_LAM_Aranami_L = F_mass_kind==1.and..not.Grd_yinyang_L.and.Adv_BC_LAM_flux==1
      BC_LAM_zlf_L     = F_mass_kind==1.and..not.Grd_yinyang_L.and.Adv_BC_LAM_flux==2

      if (BC_LAM_Aranami_L) then

         Adv_BC_LAM_flux_n = 1

         call adv_get_ij0n_ext (i0_2,in_2,j0_2,jn_2,2) !EXTENSION (CFL)

         fld_adv = 0. ; cub_o = 0. ; cub_i = 0.

      end if

      if (BC_LAM_zlf_L) call adv_get_ij0n_ext (i0_0,in_0,j0_0,jn_0,0) !CORE

      !Initializations LCSL advection
      !------------------------------
      Adv_LCSL_n      = max(F_mass_kind - 1,0)
      Adv_LCSL_mono_n = F_mono_kind

      nbpts = F_ni*F_nj*F_nk

      !Bermejo-Conde LAM Flux ZLF: Set ZERO piloting conditions outside EXTENSION (CFL)
      !--------------------------------------------------------------------------------
      if (BC_LAM_zlf_L) call adv_BC_LAM_zlf (fld_in,empty,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,1)

      call rpn_comm_xch_halox( fld_in, F_minx, F_maxx,F_miny, F_maxy , &
       F_ni, F_nj, F_nk, adv_halox, adv_haloy, G_periodx, G_periody  , &
       fld_adv, adv_lminx,adv_lmaxx,adv_lminy,adv_lmaxy, F_ni, 0)

      !Preparation Mass Conservation/Shape-Preserving
      !----------------------------------------------
      if (Adv_Mass_Cons_tr_L) then

         nullify(fld_cub,fld_mono,fld_lin,fld_min,fld_max)
         err = gmm_get(gmmk_cub_s ,fld_cub )
         err = gmm_get(gmmk_mono_s,fld_mono)
         err = gmm_get(gmmk_lin_s ,fld_lin )
         err = gmm_get(gmmk_min_s ,fld_min )
         err = gmm_get(gmmk_max_s ,fld_max )

         fld_cub (:,:,:) = fld_out(:,:,:) !Copy piloting conditions
         fld_mono(:,:,:) = fld_out(:,:,:) !Copy piloting conditions

         !Bermejo-Conde LAM Flux Aranami: Apply mask_o/mask_i to Tracer for Flux calculations
         !-----------------------------------------------------------------------------------
         if (BC_LAM_Aranami_L) call adv_BC_LAM_mask (adv_o,adv_i,fld_adv,adv_lminx,adv_lmaxx,adv_lminy,adv_lmaxy,F_nk)

         !LCSL advection: Convert Mixing ratio to Density (on Advection grid)
         !-------------------------------------------------------------------
         if (Adv_LCSL_n>0.and.Adv_LCSL_n/=3) then
            call adv_LCSL_mix_2_den (adv_rho,fld_in,adv_lminx,adv_lmaxx,adv_lminy,adv_lmaxy, &
                                     F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0,1)
         end if

      else
         fld_cub  => empty
         fld_mono => empty
         fld_lin  => empty
         fld_min  => empty
         fld_max  => empty
      end if

      call adv_tricub_lag3d (wrkc, w_mono, w_lin, w_min, w_max,    &
                             fld_adv, w_cub_o,adv_o,w_cub_i,adv_i, &
                             F_capx, F_capy, F_capz,               &
                             nbpts, F_nind, F_ii, F_k0, F_nk,      &
                             mono_L, F_lev_S)

      !----------------------------
      !Standard cubic interpolation
      !----------------------------
      if (.not.Adv_Mass_Cons_tr_L) then

!$omp parallel do private(k) shared(wrkc)
         do k = F_k0, F_nk
            Fld_out(F_i0:F_in,F_j0:F_jn,k) = wrkc(F_i0:F_in,F_j0:F_jn,k)
         end do
!$omp end parallel do

      !----------------------------------
      !Mass Conservation/Shape-Preserving
      !----------------------------------
      else

!$omp parallel do private(k) shared(wrkc,w_mono,w_lin,w_min,w_max)
         do k = F_k0, F_nk
            Fld_cub (F_i0:F_in,F_j0:F_jn,k) = wrkc  (F_i0:F_in,F_j0:F_jn,k)
            Fld_mono(F_i0:F_in,F_j0:F_jn,k) = w_mono(F_i0:F_in,F_j0:F_jn,k)
            Fld_lin (F_i0:F_in,F_j0:F_jn,k) = w_lin (F_i0:F_in,F_j0:F_jn,k)
            Fld_min (F_i0:F_in,F_j0:F_jn,k) = w_min (F_i0:F_in,F_j0:F_jn,k)
            Fld_max (F_i0:F_in,F_j0:F_jn,k) = w_max (F_i0:F_in,F_j0:F_jn,k)
         end do
!$omp end parallel do

         if (BC_LAM_Aranami_L) then

!$omp parallel do private(k) &
!$omp shared(cub_o,cub_i,w_cub_o,w_cub_i,i0_2,in_2,j0_2,jn_2)
         do k = F_k0, F_nk
            cub_o(i0_2:in_2,j0_2:jn_2,k)= w_cub_o(i0_2:in_2,j0_2:jn_2,k)
            cub_i(i0_2:in_2,j0_2:jn_2,k)= w_cub_i(i0_2:in_2,j0_2:jn_2,k)
         end do
!$omp end parallel do

         end if

         !Bermejo-Conde LAM Flux ZLF
         !--------------------------
         if (BC_LAM_zlf_L) then

             !Keep piloting conditions outside CORE
             !-------------------------------------
             call adv_BC_LAM_zlf (fld_out,store_pilot,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,2)

             !ZERO piloting conditions outside EXTENSION (BCS_BASE)
             !-----------------------------------------------------
             call adv_BC_LAM_zlf (fld_cub,empty,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,3)

         end if

         !LCSL advection: Convert Density to Mixing ratio
         !-----------------------------------------------
         if (Adv_LCSL_n>0.and.Adv_LCSL_n/=3) call adv_LCSL_mix_2_den (empty,fld_cub,adv_lminx,adv_lmaxx,adv_lminy,adv_lmaxy, &
                                                                      F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0,2)

         !Apply a posteriori Mass-fixer/Shape-Preserving schemes or Finalization
         !----------------------------------------------------------------------
         call adv_a_posteriori (F_name, fld_out, fld_cub, fld_mono, fld_lin, fld_min, fld_max, &
                                fld_in, cub_o, cub_i, F_minx, F_maxx, F_miny, F_maxy, F_nk,    &
                                F_i0, F_in, F_j0, F_jn, F_k0, F_mono_kind, F_mass_kind)

         !Bermejo-Conde LAM Flux ZLF: Reset piloting conditions outside CORE
         !------------------------------------------------------------------
         if (BC_LAM_zlf_L) call adv_BC_LAM_zlf (fld_out,store_pilot,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk,4)

      end if

      Adv_Mass_Cons_tr_L = .false. !Reset to DEFAULT
      Adv_LCSL_n         = 0       !Reset to DEFAULT
      Adv_LCSL_mono_n    = 0       !Reset to DEFAULT
      Adv_BC_LAM_flux_n  = 0       !Reset to DEFAULT
!
!---------------------------------------------------------------------
!
      return
      end subroutine adv_cubic
