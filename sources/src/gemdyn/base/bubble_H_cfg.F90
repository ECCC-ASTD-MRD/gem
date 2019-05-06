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
!**function bubble_H_cfg - to configure bubble case for Height Coord.
!
      integer function bubble_H_cfg()
      use VERgrid_options
      use bubble_options
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer k
      real*8 height_8,htop_8
!
!     ---------------------------------------------------------------
!
      bubble_H_cfg = -1

      ! establish vertical grid configuration

      G_nk   = bubble_nk
      htop_8 = G_nk*bubble_dz

      if ( hyb_H(1) < 0 ) then
         do k=1,G_nk
            height_8=htop_8*(1.d0-(dble(k)-.5d0)/G_nk)
            hyb_H(k)=height_8
         end do
      else
         do k=1024,1,-1
            if(hyb_H(k) < 0) G_nk=k-1
         end do
      end if

      bubble_H_cfg = 1

      return
      end function bubble_H_cfg

