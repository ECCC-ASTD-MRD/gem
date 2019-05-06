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

!**s/r set_oprz_H - Computes vertical operators

      subroutine fislh_set_oprz (F_errcode)
      use dcst
      use glb_ld
      use lun
      use opr
      use ver
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_errcode

!Author: Claude Girard, July 2017

      integer :: k, AA, BB, CC
      real*8, dimension(G_nk) :: r_8
      real*8, dimension(G_nk*G_nk) :: br_8, bl_8
      real*8, parameter :: one = 1.d0, half = 0.5d0
!     __________________________________________________________________
!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Compute the vertical operators: tri-diagnonal matrices
!                  O(AA+k): lower diagonal
!                  O(BB+k):       diagonal
!                  O(CC+k): upper diagonal
!           O(AA+k)*P(k-1)+O(BB+k)*P(k)+O(CC+k)*P(k+1)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (Lun_out > 0) write(Lun_out,1000)

      AA = 0
      BB = G_nk
      CC = G_nk*2
!
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

!     ~~~~~~~~~~~~~~~~~~~~~~
!     Second Derivative TERM: D gama D
!     ~~~~~~~~~~~~~~~~~~~~~~

      Opr_opszp2_8(AA+1)    =0.d0
      Opr_opszp2_8(BB+1)    =-gama_8*Ver_idz_8%t(1)*Ver_idz_8%m(1)
      Opr_opszp2_8(CC+1)    =+gama_8*Ver_idz_8%t(1)*Ver_idz_8%m(1)
      do k = 2, G_nk-1
         Opr_opszp2_8(AA+k) =+gama_8*Ver_idz_8%t(k-1)*Ver_idz_8%m(k)
         Opr_opszp2_8(BB+k) =-gama_8*Ver_idz_8%t(k-1)*Ver_idz_8%m(k) &
                             -gama_8*Ver_idz_8%t(k  )*Ver_idz_8%m(k)
         Opr_opszp2_8(CC+k) =+gama_8*Ver_idz_8%t(k  )*Ver_idz_8%m(k)
      end do
      Opr_opszp2_8(AA+G_nk) =+gama_8*Ver_idz_8%t(G_nk-1)*Ver_idz_8%m(G_nk)
      Opr_opszp2_8(BB+G_nk) =-gama_8*Ver_idz_8%t(G_nk-1)*Ver_idz_8%m(G_nk)
      Opr_opszp2_8(CC+G_nk) = 0.d0

!     ~~~~~~~~~~~~~~~~~~~~~~
!     First Derivative TERMS: - M gama*epsi D - D gama*mu M
!     ~~~~~~~~~~~~~~~~~~~~~~

      Opr_opszpl_8(AA+1)    = 0.d0
      Opr_opszpl_8(BB+1)    =+gama_8*(epsi_8*Ver_wp_8%m(1)*Ver_idz_8%t(1) - half*mu_8*Ver_idz_8%m(1))
      Opr_opszpl_8(CC+1)    =-gama_8*(epsi_8*Ver_wp_8%m(1)*Ver_idz_8%t(1) + half*mu_8*Ver_idz_8%m(1))
      do k = 2, G_nk-1
         Opr_opszpl_8(AA+k) =+gama_8*(epsi_8*Ver_wm_8%m(k)*Ver_idz_8%t(k-1) + half*mu_8*Ver_idz_8%m(k))
         Opr_opszpl_8(BB+k) =-gama_8*(epsi_8*Ver_wm_8%m(k)*Ver_idz_8%t(k-1) - half*mu_8*Ver_idz_8%m(k)) &
                             +gama_8*(epsi_8*Ver_wp_8%m(k)*Ver_idz_8%t(k)   - half*mu_8*Ver_idz_8%m(k))
         Opr_opszpl_8(CC+k) =-gama_8*(epsi_8*Ver_wp_8%m(k)*Ver_idz_8%t(k)   + half*mu_8*Ver_idz_8%m(k))
      end do
      Opr_opszpl_8(AA+G_nk) =+gama_8*(epsi_8*Ver_wmA_8(G_nk)*Ver_idz_8%t(G_nk-1)+half*mu_8*Ver_idz_8%m(G_nk))
      Opr_opszpl_8(BB+G_nk) =-gama_8*(epsi_8*Ver_wmA_8(G_nk)*Ver_idz_8%t(G_nk-1)-half*mu_8*Ver_idz_8%m(G_nk))
      Opr_opszpl_8(CC+G_nk) = 0.d0
!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Double Average - Constant TERMS: epsi M gama*mu M - gg
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      Opr_opszpm_8(AA+1)    = 0.d0
      Opr_opszpm_8(BB+1)    =+half*gama_8*epsi_8*mu_8*Ver_wp_8%m(1) - gg_8
      Opr_opszpm_8(CC+1)    =+half*gama_8*epsi_8*mu_8*Ver_wp_8%m(1)
      do k = 2, G_nk-1
         Opr_opszpm_8(AA+k) =+half*gama_8*epsi_8*mu_8*Ver_wm_8%m(k)
         Opr_opszpm_8(BB+k) =+half*gama_8*epsi_8*mu_8*Ver_wm_8%m(k) &
                             +half*gama_8*epsi_8*mu_8*Ver_wp_8%m(k) - gg_8
         Opr_opszpm_8(CC+k) =+half*gama_8*epsi_8*mu_8*Ver_wp_8%m(k)
      end do
      Opr_opszpm_8(AA+G_nk) =+half*gama_8*epsi_8*mu_8*Ver_wmA_8(G_nk)
      Opr_opszpm_8(BB+G_nk) =+half*gama_8*epsi_8*mu_8*Ver_wmA_8(G_nk) - gg_8
      Opr_opszpm_8(CC+G_nk) = 0.d0

!     ---------------------------------------------------
!     Compute eigenvalues and eigenvector in the vertical
!     ---------------------------------------------------

      call preverln ( r_8, bl_8, br_8, G_nk, G_nk, F_errcode)
!
!     transfer results back in Opr_* output arrays
!
      Opr_zevec_8 = br_8 ; Opr_lzevec_8 = bl_8 ; Opr_zeval_8 = r_8

 1000 format( &
      /,'COMPUTE EIGENVALUES AND EIGENVECTOR IN THE VERTICAL (S/R SET_OPRZ_H)', &
      /,'===========================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
