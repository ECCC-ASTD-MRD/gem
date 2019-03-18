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

!** s/r cal_div - Computes horizontal divergence

      subroutine cal_div ( F_DD, F_uu, F_vv  , &
                           F_filtdd, F_coefdd, &
                           Minx,Maxx,Miny,Maxy,Nk )

      use geomh
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer  F_filtdd, Minx,Maxx,Miny,Maxy,Nk
      real     F_DD (Minx:Maxx,Miny:Maxy,Nk), &
               F_uu (Minx:Maxx,Miny:Maxy,Nk), &
               F_vv (Minx:Maxx,Miny:Maxy,Nk), F_coefdd
!author
!    Michel Desgagne   - summer 2015
!
!revision
! v4_80 - Desgagne M.       - initial version


      integer i, j, k, i0, in, j0, jn
!
!----------------------------------------------------------------------
!
      i0 = 1
      in = l_niu
      j0 = 1
      jn = l_njv
      if (l_west ) i0 = 2
      if (l_south) j0 = 2

      do k = 1 , Nk
         do j = j0, jn
         do i = i0, in
            F_DD(i,j,k) = &
            ((F_uu(i,j,k) - F_uu(i-1,j,k)) * geomh_invDX_8(j)) &
          + ( (F_vv(i,j  ,k)*geomh_cyv_8(j  )   &
             - F_vv(i,j-1,k)*geomh_cyv_8(j-1))  &
             * geomh_invDY_8 * geomh_invcy_8(j) )
         end do
         end do
         F_DD(1:i0-1,:,k) = 0. ; F_DD(in+1:l_ni,:,k)= 0.
         F_DD(:,1:j0-1,k) = 0. ; F_DD(:,jn+1:l_nj,k)= 0.
      end do

      if (F_filtdd > 0) call filter2 ( F_DD, F_filtdd,F_coefdd, &
                                  l_minx,l_maxx,l_miny,l_maxy,Nk )
!
!----------------------------------------------------------------------
!
      return
      end
