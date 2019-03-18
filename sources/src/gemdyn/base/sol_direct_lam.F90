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
!*s/r sol_direct_lam - Solution of the elliptic problem for LAM grids
!
      subroutine sol_direct_lam ( F_sol_8, F_rhs_8, F_ni, F_nj, F_nk )
      use fft
      use gem_options
      use glb_ld
      use glb_pil
      use ldnh
      use opr
      use ptopo
      use sol
      use trp
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_ni, F_nj, F_nk
      real*8, dimension(F_ni,F_nj,F_nk), intent(in) :: F_rhs_8
      real*8, dimension(F_ni,F_nj,F_nk), intent(out) :: F_sol_8

!author
!     Desgagne   Spring 2013
!
      integer :: Gni, i, j, dimens
      real*8, dimension(:), allocatable :: wk_evec_8
      real*8, dimension((ldnh_maxy -ldnh_miny +1)*(trp_12smax-trp_12smin+1)*(G_ni+Ptopo_npex  )) :: fdg1
      real*8, dimension((trp_12smax-trp_12smin+1)*(trp_22max -trp_22min +1)*(G_nj+Ptopo_npey  )) :: fdg2
      real*8, dimension((ldnh_maxy -ldnh_miny +1)*(trp_12smax-trp_12smin+1)*(G_ni+2+Ptopo_npex)) :: fdwfft
!
!     ---------------------------------------------------------------
!
      if (Fft_fast_L) then

         call sol_fft_lam ( F_sol_8, F_rhs_8,                               &
                           ldnh_maxx, ldnh_maxy, ldnh_nj,                   &
                           trp_12smax, trp_12sn, trp_22max, trp_22n,        &
                           G_ni, G_nj, G_nk, sol_nk, Ptopo_npex, Ptopo_npey,&
                           Sol_ai_8, Sol_bi_8, Sol_ci_8, fdg2, fdwfft )

      else

         Gni    = G_ni - Lam_pil_w - Lam_pil_e
         dimens = Gni * Gni

         allocate ( wk_evec_8(dimens) )

         do j=1,Gni
            do i=1,Gni
               wk_evec_8((j-1)*Gni+i) = Opr_xevec_8((j+Lam_pil_w-1)*G_ni+i+Lam_pil_w)
            end do
         end do

         call sol_mxma ( F_sol_8, F_rhs_8, wk_evec_8,                      &
                         ldnh_maxx, ldnh_maxy, ldnh_nj, dimens,            &
                         trp_12smax, trp_12sn, trp_22max, trp_22n,         &
                         G_ni, G_nj, G_nk, sol_nk, Ptopo_npex, Ptopo_npey, &
                         Sol_ai_8, Sol_bi_8, Sol_ci_8, fdg1, fdg2, fdwfft)

         deallocate (wk_evec_8)

      end if
!
!     ---------------------------------------------------------------
!
      return
      end
