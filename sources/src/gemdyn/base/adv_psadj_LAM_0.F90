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

!**s/r adv_psadj_LAM_0 - Estimate FLUX_out/FLUX_in based on Aranami et al. (2015) and
!                        Call psadj_LAM to adjust surface pressure

      subroutine adv_psadj_LAM_0 ()

      use adv_pos
      use lam_options
      use glb_ld
      use adv_options

      implicit none

#include <arch_specific.hf>

      !object
      !=====================================================================================
      !     Estimate FLUX_out/FLUX_in based on Aranami et al.,QJRMS,141,1795-1803 (2015) and
      !     Call psadj_LAM to adjust surface pressure
      !=====================================================================================

      !-----------------------------------------------------------------

      integer :: k, k0, nbpts, i0_2, in_2, j0_2, jn_2, empty_i1
      real,   dimension(l_ni,l_nj,l_nk) :: w_cub_o,w_cub_i
      real,   dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: cub_o,cub_i
      real,   dimension(1,1,1) :: empty
      integer,dimension(1)     :: empty_i2

      !-----------------------------------------------------------------

      k0 = Lam_gbpil_t+1

      call adv_get_ij0n_ext (i0_2,in_2,j0_2,jn_2,2) !EXTENSION (CFL)

      nbpts = l_ni*l_nj*l_nk

      cub_o = 0. ; cub_i = 0.

      adv_pos_reset(1) = 1

      adv_BC_LAM_flux_n = 2

      !Estimate FLUX_out/FLUX_in based on Aranami et al. (2015)
      !--------------------------------------------------------
      call adv_tricub_lag3d (empty, empty, empty, empty, empty, empty, &
                             w_cub_o, adv_BC_LAM_mask_o, w_cub_i, adv_BC_LAM_mask_i, pxt, pyt, pzt, nbpts, &
                             empty_i1, empty_i2,  k0, l_nk, .false., 't')

      adv_BC_LAM_flux_n = 0 !Reset to DEFAULT

!$omp parallel do private(k) &
!$omp shared(k0,i0_2,in_2,j0_2,jn_2,cub_o,cub_i,w_cub_o,w_cub_i)
      do k = k0,l_nk
         cub_o(i0_2:in_2,j0_2:jn_2,k)= w_cub_o(i0_2:in_2,j0_2:jn_2,k)
         cub_i(i0_2:in_2,j0_2:jn_2,k)= w_cub_i(i0_2:in_2,j0_2:jn_2,k)
      end do
!$omp end parallel do

      !Adjust surface pressure
      !-----------------------
      call psadj_LAM (cub_o,cub_i,l_minx,l_maxx,l_miny,l_maxy,l_nk,k0)

      return

      end subroutine adv_psadj_LAM_0
