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

      subroutine adv_tricub_lag3d (F_cub, F_mono, F_lin, F_min, F_max, F_in,        &
                                   F_cub_o, F_in_o, F_cub_i, F_in_i, F_x, F_y, F_z, &
                                   F_num, F_nind, ii, F_k0, F_nk, F_mono_L,  F_lev)
      use glb_ld
      use adv
      use adv_grid
      use adv_interp
      use adv_options
      use gem_timing
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_lev ! m/t : Momemtum/thermo level
      integer, intent(in) :: F_num ! number points
      integer, intent(in) :: F_nk ! number of vertical levels
      integer, intent(in) :: F_k0 ! scope of operator
      logical, intent(in) :: F_mono_L ! .true. monotonic interpolation
      real,dimension(F_num), intent(in)  :: F_x, F_y, F_z ! interpolation target x,y,z coordinates
      real,dimension(*),     intent(in)  :: F_in          ! field to interpolate
      real,dimension(*),     intent(in)  :: F_in_o,F_in_i ! field to interpolate (FLUX_out/FLUX_in)
      real,dimension(F_num), intent(out) :: F_cub ! High-order SL solution
      real,dimension(F_num), intent(out) :: F_mono! High-order monotone SL solution
      real,dimension(F_num), intent(out) :: F_lin ! Low-order SL solution
      real,dimension(F_num), intent(out) :: F_min ! MIN over cell
      real,dimension(F_num), intent(out) :: F_max ! MAX over cell
      real,dimension(F_num), intent(out) :: F_cub_o ! High-order SL solution FLUX_out
      real,dimension(F_num), intent(out) :: F_cub_i ! High-order SL solution FLUX_in
      integer , intent(in) :: F_nind
      integer , dimension(F_nind*4), intent(in)  :: ii            ! pre-computed indices to be used in: adv_tricub_lag3d_loop

   !@revisions
   !  2012-05,  Stephane Gaudreault: code optimization
   !  2016-01,  Monique Tanguay    : GEM4 Mass-Conservation
   !@objective Tri-cubic interp: Lagrange 3d (Based on adx_tricub v3.1.1) (MASS-CONSERVATION)

      real*8, dimension(:), pointer, contiguous :: p_bsz_8, p_zbc_8, p_zabcd_8, &
                                                   p_zbacd_8, p_zcabd_8, p_zdabc_8, &
                                                   p_zxabcde_8, p_zaxbcde_8, p_zbxacde_8, &
                                                   p_zcxabde_8, p_zdxabce_8, p_zexabcd_8
      character(len=12) intp_S
!
!---------------------------------------------------------------------
!
      if ( trim(Adv_component_S) == 'TRAJ' ) then
         call gemtime_start (37, 'ADV_LAG3D', 34)
      else if ( trim(Adv_component_S) == 'INTP_RHS' ) then
         call gemtime_start (38, 'ADV_LAG3D', 31)
      else
         call gemtime_start (39, 'ADV_LAG3D', 27)
      end if

      intp_S = 'CUBIC'
      if ( trim(Adv_component_S) == 'INTP_TR') intp_S = adv_intp_S

      if (F_lev/='t'.and.intp_S=='QUINTIC') call handle_error(-1,'adv_tricub_lag3d','p_zxabcde NOT completed')

      if (F_lev == 'm') then
         p_bsz_8   => adv_bsz_8%m
         p_zabcd_8 => adv_zabcd_8%m
         p_zbacd_8 => adv_zbacd_8%m
         p_zcabd_8 => adv_zcabd_8%m
         p_zdabc_8 => adv_zdabc_8%m
         p_zbc_8   => adv_zbc_8%m
      else if (F_lev  == 't') then
         p_bsz_8   => adv_bsz_8%t
         p_zabcd_8 => adv_zabcd_8%t
         p_zbacd_8 => adv_zbacd_8%t
         p_zcabd_8 => adv_zcabd_8%t
         p_zdabc_8 => adv_zdabc_8%t
         p_zbc_8   => adv_zbc_8%t
         !
         p_zxabcde_8 => adv_zxabcde_8%t
         p_zaxbcde_8 => adv_zaxbcde_8%t
         p_zbxacde_8 => adv_zbxacde_8%t
         p_zcxabde_8 => adv_zcxabde_8%t
         p_zdxabce_8 => adv_zdxabce_8%t
         p_zexabcd_8 => adv_zexabcd_8%t
      else if (F_lev == 'x') then
         p_bsz_8   => adv_bsz_8%x
         p_zabcd_8 => adv_zabcd_8%x
         p_zbacd_8 => adv_zbacd_8%x
         p_zcabd_8 => adv_zcabd_8%x
         p_zdabc_8 => adv_zdabc_8%x
         p_zbc_8   => adv_zbc_8%x
      end if

      !---------------------------------------------------------------------------
      !Bermejo-Conde LAM: Estimate FLUX_out/FLUX_in based on Aranami et al. (2015)
      !---------------------------------------------------------------------------
      if (Adv_BC_LAM_flux_n>0) call adv_BC_LAM_Aranami (F_cub_o, F_in_o, F_cub_i, F_in_i, &
                                                        F_x, F_y, F_z, F_num, F_k0, F_nk, F_lev)

      if (Adv_BC_LAM_flux_n == 2) return

      !----------------------------
      !Standard Cubic interpolation
      !----------------------------
      if (.not.Adv_Mass_Cons_tr_L.and.intp_S=='CUBIC') then

         if (F_mono_L) then
            call adv_tricub_mono_loop()
         else
            call adv_tricub_loop()
         end if

      !--------------------------------
      !Standard Quintic_V interpolation
      !--------------------------------
      elseif (.not.Adv_Mass_Cons_tr_L.and.intp_S=='QUINTIC') then

         if (F_mono_L) then
            call adv_triqut_mono_loop()
         else
            call adv_triqut_loop()
         end if

      !--------------------------------------------------------------------
      !Standard Cubic interpolation with MIN/MAX/LIN for Bermejo-Conde/ILMC
      !--------------------------------------------------------------------
      elseif (Adv_Mass_Cons_tr_L.and.Adv_LCSL_n==0.and.intp_S=='CUBIC') then

         call adv_tricub_conserv_loop()

      !------------------------------------------------------------------------
      !Standard Quintic_V interpolation with MIN/MAX/LIN for Bermejo-Conde/ILMC
      !------------------------------------------------------------------------
      elseif (Adv_Mass_Cons_tr_L.and.Adv_LCSL_n==0.and.intp_S=='QUINTIC') then

         call adv_triqut_conserv_loop()

      !---------------------------------------------------
      !Local Conservative Semi-Lagrangian (LCSL) Advection
      !---------------------------------------------------
      else

         call adv_LCSL (F_cub, l_ni, l_nj, F_nk)

      end if

      if ( trim(Adv_component_S) == 'TRAJ' ) then
         call gemtime_stop (37)
      else if ( trim(Adv_component_S) == 'INTP_RHS' ) then
         call gemtime_stop (38)
      else
         call gemtime_stop (39)
      end if
!
!---------------------------------------------------------------------
!
      return


contains

      include 'adv_tricub_loop.inc'
      include 'adv_tricub_mono_loop.inc'
      include 'adv_tricub_conserv_loop.inc'
      include 'adv_triqut_loop.inc'
      include 'adv_triqut_mono_loop.inc'
      include 'adv_triqut_conserv_loop.inc'

      end subroutine adv_tricub_lag3d
