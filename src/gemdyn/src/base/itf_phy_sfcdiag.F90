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

      subroutine itf_phy_sfcdiag (F_dest, minx, maxx, miny, maxy, &
                                  F_var_S, F_status, F_quiet_L)
      use phy_itf, only: phy_get
      use rmn_gmm
      use HORgrid_options
      use glb_ld
      use outp
      implicit none
#include <arch_specific.hf>

      character(len=*) :: F_var_S
      logical :: F_quiet_L
      integer :: minx, maxx, miny, maxy, F_status
      real, dimension(minx:maxx,miny:maxy) :: F_dest

      real,dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk+1), target :: w4
      real, dimension(:,:  ), pointer :: qdiag
      real, dimension(:,:,:), pointer :: ptr3d
      integer, dimension(3), save :: lijk = [ -1, -1, -1 ]
      integer, dimension(3), save :: uijk = [ -1, -1, -1 ]
!
!-----------------------------------------------------------------
!
      if (trim(F_var_S)=='TR/HU:P') then

         F_status = gmm_get(gmmk_diag_hu_s,qdiag)
         if (F_status == 0) F_dest(:,:) = qdiag(:,:)

      else

         ptr3d => w4(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:G_nk+1)
         F_status = phy_get (ptr3d, F_var_S, F_npath='VO', F_bpath='D'    ,&
                       F_start=lijk, F_end=uijk, &
                       F_quiet=F_quiet_L)
         if (F_status == 0) F_dest(:,:) = w4(:,:,G_nk+1)

      end if
!
!-----------------------------------------------------------------
!
      return
      end subroutine itf_phy_sfcdiag
