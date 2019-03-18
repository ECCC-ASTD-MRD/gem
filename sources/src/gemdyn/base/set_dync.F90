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

!**   s/r set_dync - initialize the dynamics model configuration

      subroutine set_dync (F_check_and_stop_L, F_errcode)
      use dcst
      use dynkernel_options
      use matvec, only: matvec_init
      use gmm_vt1
      use HORgrid_options
      use dyn_fisl_options
      use lam_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use ver
      use gmm_itf_mod
      use sol
      use opr
      implicit none
#include <arch_specific.hf>

      logical F_check_and_stop_L
      integer F_errcode


      integer k,err,k0
      real tmean(G_nk)
      real*8  w1, w2
      real*8, parameter :: zero=0.d0, one=1.d0, half=.5d0
!
!     ---------------------------------------------------------------

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H') then
         return
      end if

      if( Cstv_Tstr_8 < 0. ) then
         ! TSTAR variable in the vertical
         err = gmm_get(gmmk_tt1_s,tt1)
         call estimate_tmean (tmean,tt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         Ver_Tstar_8%t(1:G_nk) = tmean(1:G_nk)
      else
         do k=1,G_nk
            Ver_Tstar_8%t(k) = Cstv_Tstr_8
         end do
      end if

      Ver_Tstar_8%m(1) = Ver_Tstar_8%t(1)
      do k=2,G_nk
         Ver_Tstar_8%m(k) = Ver_wp_8%m(k)*Ver_Tstar_8%t(k) + &
                            Ver_wm_8%m(k)*Ver_Tstar_8%t(k-1)
      end do
      Ver_Tstar_8%m(G_nk+1) = Ver_Tstar_8%t(G_nk)

      Ver_fistr_8(G_nk+1)= 0.d0
      do k = G_nk, 1, -1
         Ver_fistr_8(k) = Ver_fistr_8(k+1) - &
           Rgasd_8*Ver_Tstar_8%t(k)*(Ver_z_8%m(k)-Ver_z_8%m(k+1))
      end do

      do k=1,G_nk
         Ver_epsi_8(k)=Rgasd_8*Ver_Tstar_8%t(k)*Ver_igt2_8
         Ver_gama_8(k)=Cstv_invT_8*Cstv_invT_m_8/ &
             (Rgasd_8*Ver_Tstar_8%t(k)*(cappa_8+Ver_epsi_8(k)))
      end do

      Cstv_hco0_8 = Dcst_rayt_8**2
      Cstv_hco1_8 = zero
      Cstv_hco2_8 = one

      Ver_alfat_8 = one
      Ver_cst_8   = zero
      Ver_cstp_8  = zero

      k0=1+Lam_gbpil_T
      if(Schm_opentop_L) then
         w1 = Ver_idz_8%t(k0-1)* &
               (Ver_idz_8%m(k0)-Ver_wm_8%m(k0)*(one+Ver_epsi_8(k0-1))) &
            + half*Ver_epsi_8(k0-1)* &
               (Ver_idz_8%m(k0)-Ver_wm_8%m(k0)*(one-cappa_8))
         w2 = one/(Ver_idz_8%t(k0-1)+Ver_epsi_8(k0-1)*half)
         Ver_alfat_8 = (Ver_idz_8%t(k0-1) - Ver_epsi_8(k0-1)*half) * w2
         Ver_cst_8   =                      one / Ver_gama_8(k0-1) * w2
         Ver_cstp_8  =                                          w1 * w2
      end if

      Ver_css_8   = one/Ver_gama_8(G_nk) &
                   /(Ver_idz_8%t(G_nk)+cappa_8*Ver_wpstar_8(G_nk))
      w1= Ver_wmstar_8(G_nk)*half*(Ver_gama_8(G_nk  )*Ver_epsi_8(G_nk) &
                                  -Ver_gama_8(G_nk-1)*Ver_epsi_8(G_nk-1))
      w2 = Ver_wmstar_8(G_nk)*Ver_gama_8(G_nk-1)*Ver_idz_8%t(G_nk-1)
      Ver_alfas_8 = Ver_css_8*Ver_gama_8(G_nk)*Ver_idz_8%t(G_nk) &
                  + Ver_css_8 * ( w1 + w2 )
      Ver_betas_8 = Ver_css_8 * ( w1 - w2 )
      w1=Ver_gama_8(G_nk)*Ver_idz_8%t(G_nk)*(Ver_idz_8%m(G_nk) + &
         Ver_wp_8%m(G_nk))/Ver_wpstar_8(G_nk)
      w2=((one-cappa_8)*Ver_wp_8%m(G_nk) + Ver_idz_8%m(G_nk) - &
         Ver_wpA_8(G_nk)*Ver_idz_8%t(G_nk)) &
        *Ver_gama_8(G_nk)*Ver_epsi_8(G_nk)
      Ver_cssp_8  = Ver_css_8 * ( w1 - w2 )

      Cstv_bar0_8 = zero
      Cstv_bar1_8 = one
      if(Schm_autobar_L) then
         Cstv_bar0_8 = Cstv_invT_8**2/Ver_FIstr_8(1)
         Cstv_bar1_8 = zero
         Ver_alfas_8 = one
         Ver_css_8   = zero
         Ver_cssp_8  = zero
         Cstv_hco1_8 = Cstv_bar0_8
         Cstv_hco2_8 = zero
      end if

      call set_oprz (F_errcode)

      if (F_check_and_stop_L) then
         call gem_error (F_errcode,'set_dync',&
              'VERTICAL LAYERING and TIMESTEP INCOMPATIBILITY')
         call set_sol
         if (Sol_type_S == 'ITERATIVE_3D') call matvec_init
      end if
!
!     ---------------------------------------------------------------
!
      return
      end
