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

!**s/r adz_od_BC_LAM_setup - Set F_mask_o/F_mask_i for Flux calculations based on Aranami et al. (2015)

      subroutine adz_od_BC_LAM_setup (F_mask_o,F_mask_i,F_minx,F_maxx,F_miny,F_maxy,F_ni,F_nj,F_nk)

      use adz_mem
      use gem_options 

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy !Dimension H
      integer, intent(in) :: F_ni,F_nj,F_nk              !Dimension No Halo
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_mask_o,F_mask_i

      !object
      !===================================================================================================
      !     Set F_mask_o/F_mask_i for Flux calculations based on Aranami et al.,QJRMS,141,1795-1803 (2015)
      !===================================================================================================

      !--------------------------------------------------------------------------

      integer :: i,j,k
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: fld_ONE,work

      !--------------------------------------------------------------------------

      fld_ONE = 0.

      fld_ONE(1:F_ni,1:F_nj,1:F_nk) = 1.

      !For Flux_out: Set mixing ratio = 0 on NEST-HV
      !---------------------------------------------
      work = 0.

      do k=Adz_k0t,F_nk
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               work(i,j,k) = fld_ONE(i,j,k)
            end do
         end do
      end do

      F_mask_o = work

      !For Flux_in: Set mixing ratio = 0 on CORE-HV
      !--------------------------------------------
      work = fld_ONE

      do k=Adz_k0t,F_nk
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               work(i,j,k) = 0.
            end do
         end do
      end do

      F_mask_i = work

      !Fill up HALO zone
      !-----------------
      call rpn_comm_xch_halo( F_mask_o,F_minx,F_maxx,F_miny,F_maxy,&
                              F_ni,F_nj,F_nk,G_halox,G_haloy,&
                              G_periodx,G_periody,G_ni,0)

      call rpn_comm_xch_halo( F_mask_i,F_minx,F_maxx,F_miny,F_maxy,&
                              F_ni,F_nj,F_nk,G_halox,G_haloy,&
                              G_periodx,G_periody,G_ni,0)

      return

      end subroutine adz_od_BC_LAM_setup
