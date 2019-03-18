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

!**s/r pospers - initialise upstream positions at time th as grid point positions
!

!
      subroutine pospers
!
      use gmm_vth
      use geomh
      use glb_ld
      use ver
      use type_mod
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>
!
!author
!     Alain Patoine
!
!revision
! v2_00 - Desgagne M.       - initial MPI version
! V4_10 - Plante A.         - Thermo upstream positions
!
!object
!
!arguments
!     none
!

!
      type(gmm_metadata) :: mymeta
      integer i, j, k, ijk, nij,istat
!*
!
!     ---------------------------------------------------------------
!
      nij = l_ni * l_nj
!
      istat = gmm_get(gmmk_xth_s,xth,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'pospers ERROR at gmm_get(xth)'
      istat = gmm_get(gmmk_yth_s,yth,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'pospers ERROR at gmm_get(yth)'
      istat = gmm_get(gmmk_zth_s,zth,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'pospers ERROR at gmm_get(zth)'
!
!
!
      do k = 1, l_nk
      do j = 1, l_nj
      do i = 1, l_ni
         ijk=(k-1)*nij+(j-1)*l_ni+i
         xth(ijk)  = geomh_x_8(i)
         yth(ijk)  = geomh_y_8(j)
         zth(ijk)  = Ver_z_8%m(k)
      end do
      end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end
