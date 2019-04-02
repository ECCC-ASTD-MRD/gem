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

!**s/r adz_BC_LAM_mask - Apply mask_o/mask_i to Tracer for Flux calculations based on Aranami et al. (2015)

      subroutine adz_BC_LAM_mask (F_adv_o,F_adv_i,F_adv,F_aminx,F_amaxx,F_aminy,F_amaxy,F_nk)

      use adz_options

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_aminx,F_amaxx,F_aminy,F_amaxy,F_nk
      real, dimension(F_aminx:F_amaxx,F_aminy:F_amaxy,F_nk), intent(in)  :: F_adv
      real, dimension(F_aminx:F_amaxx,F_aminy:F_amaxy,F_nk), intent(out) :: F_adv_o,F_adv_i

#include <arch_specific.hf>

      !object
      !=======================================================================================
      !     Apply mask_o/mask_i to Tracer for Flux calculations based on Aranami et al. (2015)
      !=======================================================================================

      integer :: i,j,k

!$omp parallel do private (i,j,k)
      do k=1,F_nk
      do j=F_aminy,F_amaxy
      do i=F_aminx,F_amaxx
         F_adv_o(i,j,k) = Adz_BC_LAM_mask_o(i,j,k) * F_adv(i,j,k)
         F_adv_i(i,j,k) = Adz_BC_LAM_mask_i(i,j,k) * F_adv(i,j,k)
      end do
      end do
      end do
!$omp end parallel do

      return

      end subroutine adz_BC_LAM_mask
