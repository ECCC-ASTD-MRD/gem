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

!**s/r set_oprz - Computes vertical operators

      subroutine set_oprz (F_errcode)
      use dyn_fisl_options
      use glb_ld
      use lam_options
      use lun
      use tdpack
      use ver
      use opr
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer :: F_errcode

      real(kind=REAL64), parameter :: one = 1.d0
      integer :: k, k0, AA, BB, CC
      real(kind=REAL64) :: Falfas_8, Fbetas_8, r_8(G_nk)
      real(kind=REAL64), dimension (G_nk*G_nk) :: br_8, bl_8
!     __________________________________________________________________
!
!     Compute the vertical operators: tri-diagnonal matrices
!                  O(AA+k): lower diagonal
!                  O(BB+k):       diagonal
!                  O(CC+k): upper diagonal
!           O(AA+k)*P(k-1)+O(BB+k)*P(k)+O(CC+k)*P(k+1)
!
      if (Lun_out > 0) write(Lun_out,1000)

      AA=0
      BB=G_nk
      CC=G_nk*2
      k0=1+Lam_gbpil_T
!
      Falfas_8=.5d0*Ver_wmstar_8(G_nk)+Ver_wpstar_8(G_nk)*Ver_alfas_8
      Fbetas_8=.5d0*Ver_wmstar_8(G_nk)+Ver_wpstar_8(G_nk)*Ver_betas_8

!     ~~~~~~~~~~~~~~~~~
!     Diagonal Operator
!     ~~~~~~~~~~~~~~~~~
      do k = 1, G_nk
         Opr_opszp0_8(AA+k) = 0.d0
         Opr_opszp0_8(BB+k) = one
         Opr_opszp0_8(CC+k) = 0.d0
         !        Zero in the others to start
         Opr_opszp2_8(AA+k) = 0.d0
         Opr_opszp2_8(BB+k) = 0.d0
         Opr_opszp2_8(CC+k) = 0.d0
         Opr_opszpm_8(AA+k) = 0.d0
         Opr_opszpm_8(BB+k) = 0.d0
         Opr_opszpm_8(CC+k) = 0.d0
         Opr_opszpl_8(AA+k) = 0.d0
         Opr_opszpl_8(BB+k) = 0.d0
         Opr_opszpl_8(CC+k) = 0.d0
      end do

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Second Derivative Operator
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~

      Opr_opszp2_8(AA+k0)   = 0.d0
      Opr_opszp2_8(BB+k0)   =-Ver_idz_8%t(k0 )*gama_8*Ver_idz_8%m(k0)
      Opr_opszp2_8(CC+k0)   =+Ver_idz_8%t(k0 )*gama_8*Ver_idz_8%m(k0)
      do k = k0+1, G_nk-1
         Opr_opszp2_8(AA+k) =+Ver_idz_8%t(k-1)*gama_8*Ver_idz_8%m(k)
         Opr_opszp2_8(BB+k) =-Ver_idz_8%t(k-1)*gama_8*Ver_idz_8%m(k) &
                             -Ver_idz_8%t(k)*gama_8*Ver_idz_8%m(k)
         Opr_opszp2_8(CC+k) =+Ver_idz_8%t(k)*gama_8*Ver_idz_8%m(k)
      end do
      Opr_opszp2_8(AA+G_nk) =+Ver_idz_8%m(G_nk)/Ver_wpstar_8(G_nk)* &
                              (Ver_idz_8%t(G_nk-1)*gama_8 &
                              +Ver_idz_8%t(G_nk  )*gama_8*Ver_betas_8 )
      Opr_opszp2_8(BB+G_nk) =-Ver_idz_8%m(G_nk)/Ver_wpstar_8(G_nk)* &
                             (Ver_idz_8%t(G_nk-1)*gama_8 &
                             +Ver_idz_8%t(G_nk  )*gama_8*(one-Ver_alfas_8) )
      Opr_opszp2_8(CC+G_nk) = 0.d0

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     First Derivative Operator (non symetric) Average of gamma*(1+epsi) X Derivative
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      Opr_opszpl_8(AA+k0)   = 0.d0
      Opr_opszpl_8(BB+k0)   =-Ver_wp_8%m(k0)*Ver_idz_8%t(k0)*gama_8*(one+epsi_8)
      Opr_opszpl_8(CC+k0)   =+Ver_wp_8%m(k0)*Ver_idz_8%t(k0)*gama_8*(one+epsi_8)
      do k = k0+1, G_nk-1
         Opr_opszpl_8(AA+k) =-Ver_wm_8%m(k)*Ver_idz_8%t(k-1)*gama_8*(one+epsi_8)
         Opr_opszpl_8(BB+k) =+Ver_wm_8%m(k)*Ver_idz_8%t(k-1)*gama_8*(one+epsi_8) &
                             -Ver_wp_8%m(k)*Ver_idz_8%t(k)*gama_8*(one+epsi_8)
         Opr_opszpl_8(CC+k) =+Ver_wp_8%m(k)*Ver_idz_8%t(k)*gama_8*(one+epsi_8)
      end do
      Opr_opszpl_8(AA+G_nk) =-Ver_wmA_8(G_nk)*Ver_idz_8%t(G_nk-1)*gama_8*(one+epsi_8) &
                             +Ver_betas_8 &
                             *Ver_wpA_8(G_nk)*Ver_idz_8%t(G_nk)*gama_8*(one+epsi_8)
      Opr_opszpl_8(BB+G_nk) =+Ver_wmA_8(G_nk)*Ver_idz_8%t(G_nk-1)*gama_8*(one+epsi_8) &
                             -(one-Ver_alfas_8) &
                             *Ver_wpA_8(G_nk)*Ver_idz_8%t(G_nk)*gama_8*(one+epsi_8)
      Opr_opszpl_8(CC+G_nk) = 0.d0
!
!     substracting Derivative of gamma*epsi X Average
      Opr_opszpl_8(BB+k0) = Opr_opszpl_8(BB+k0) - .5d0*Ver_idz_8%m(k0)*gama_8*epsi_8
      Opr_opszpl_8(CC+k0) = Opr_opszpl_8(CC+k0) - .5d0*Ver_idz_8%m(k0)*gama_8*epsi_8
      do k = k0+1, G_nk-1
         Opr_opszpl_8(AA+k) = Opr_opszpl_8(AA+k) + .5d0*Ver_idz_8%m(k)*gama_8*epsi_8
         Opr_opszpl_8(BB+k) = Opr_opszpl_8(BB+k) + .5d0*Ver_idz_8%m(k)*gama_8*epsi_8 &
                                                 - .5d0*Ver_idz_8%m(k)*gama_8*epsi_8
         Opr_opszpl_8(CC+k) = Opr_opszpl_8(CC+k) - .5d0*Ver_idz_8%m(k)*gama_8*epsi_8
      end do
      Opr_opszpl_8(AA+G_nk) = Opr_opszpl_8(AA+G_nk) &
             - Ver_idz_8%m(G_nk)/Ver_wpstar_8(G_nk) &
             * (gama_8*epsi_8*Fbetas_8 &
               -.5d0*gama_8*epsi_8)
      Opr_opszpl_8(BB+G_nk) = Opr_opszpl_8(BB+G_nk) &
             - Ver_idz_8%m(G_nk)/Ver_wpstar_8(G_nk) &
             * (gama_8*epsi_8*Falfas_8 &
               -.5d0*gama_8*epsi_8)

!     ~~~~~~~~~~~~~~~~~~~~~~~
!     Double average operator
!     ~~~~~~~~~~~~~~~~~~~~~~~
!
      Opr_opszpm_8(AA+k0)   = 0.d0
      Opr_opszpm_8(BB+k0)   = Ver_wp_8%m(k0)*.5d0*gama_8*epsi_8
      Opr_opszpm_8(CC+k0)   = Ver_wp_8%m(k0)*.5d0*gama_8*epsi_8
      do k = k0+1, G_nk-1
         Opr_opszpm_8(AA+k) = Ver_wm_8%m(k)*.5d0*gama_8*epsi_8
         Opr_opszpm_8(BB+k) = Ver_wm_8%m(k)*.5d0*gama_8*epsi_8 &
                            + Ver_wp_8%m(k)*.5d0*gama_8*epsi_8
         Opr_opszpm_8(CC+k) = Ver_wp_8%m(k)*.5d0*gama_8*epsi_8
      end do
      Opr_opszpm_8(AA+G_nk) = Ver_wmA_8(G_nk)*.5d0*gama_8*epsi_8 &
                            + Ver_wpA_8(G_nk)*Fbetas_8*gama_8*epsi_8
      Opr_opszpm_8(BB+G_nk) = Ver_wmA_8(G_nk)*.5d0*gama_8*epsi_8 &
                            + Ver_wpA_8(G_nk)*Falfas_8*gama_8*epsi_8
      Opr_opszpm_8(CC+G_nk) = 0.d0

      if(Schm_opentop_L) then
         Opr_opszp2_8(BB+k0) = Opr_opszp2_8(BB+k0) - Ver_idz_8%t(k0-1)*gama_8*Ver_idz_8%m(k0)*(one-Ver_alfat_8)
         Opr_opszpl_8(BB+k0) = Opr_opszpl_8(BB+k0) &
                             + Ver_wm_8%m(k0)*Ver_idz_8%t(k0-1)*gama_8*(one+epsi_8)*(one-Ver_alfat_8)
         Opr_opszpl_8(BB+k0) = Opr_opszpl_8(BB+k0) &
                             + .5d0*Ver_idz_8%m(k0)*gama_8*epsi_8*(one+Ver_alfat_8)
         Opr_opszpm_8(BB+k0) = Opr_opszpm_8(BB+k0) + Ver_wm_8%m(k0) &
                             * (Ver_wp_8%t(k0-1)+Ver_wm_8%t(k0-1)*Ver_alfat_8)*gama_8*epsi_8
      end if
!
!     multiplying by 1-cappa
      do k = k0, G_nk
         Opr_opszpm_8(AA+k) = Opr_opszpm_8(AA+k)*(one-cappa_8)
         Opr_opszpm_8(BB+k) = Opr_opszpm_8(BB+k)*(one-cappa_8)
         Opr_opszpm_8(CC+k) = Opr_opszpm_8(CC+k)*(one-cappa_8)
      end do

!     ---------------------------------------------------
!     Compute eigenvalues and eigenvector in the vertical
!     ---------------------------------------------------

      call preverln ( r_8, bl_8, br_8, G_nk, G_nk, F_errcode)
!
!     transfer results back in Opr_* output arrays
!
      Opr_zevec_8 = br_8 ; Opr_lzevec_8 = bl_8 ; Opr_zeval_8 = r_8

 1000 format( &
      /,'COMPUTE EIGENVALUES AND EIGENVECTOR IN THE VERTICAL (S/R SET_OPRZ)', &
      /,'===========================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
