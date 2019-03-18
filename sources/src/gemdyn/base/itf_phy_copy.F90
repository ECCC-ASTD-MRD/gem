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

!**s/r itf_phy_copy

      subroutine itf_phy_copy
      use gmm_pw
      use glb_ld
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

!author
!     Michel Desgagne - Summer 2013
!
!revision
! v4_70 - Desgagne, M.     - initial version


      integer istat,k
      real, pointer, dimension (:,:,:) :: uu_copy,vv_copy  => null()
!     ________________________________________________________________
!
      istat = gmm_get(gmmk_pw_uu_plus_s,pw_uu_plus)
      istat = gmm_get(gmmk_pw_uu_copy_s,uu_copy   )
      istat = gmm_get(gmmk_pw_vv_plus_s,pw_vv_plus)
      istat = gmm_get(gmmk_pw_vv_copy_s,vv_copy   )

!$omp parallel  private(k)
!$omp do
      do k= 1, G_nk
         uu_copy(:,:,k) = pw_uu_plus(:,:,k)
         vv_copy(:,:,k) = pw_vv_plus(:,:,k)
      end do
!$omp end do
!$omp end parallel
!     ________________________________________________________________
!
      return
      end
