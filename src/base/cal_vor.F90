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

!** s/r cal_vor - Computes horizontal vorticity and relative vorticity

      subroutine cal_vor ( F_QR,F_QQ, F_uu,F_vv          , &
                           F_filtqq, F_coefqq, F_absvor_L, &
                           Minx,Maxx,Miny,Maxy,Nk )
      use dcst
      use dynkernel_options
      use geomh
      use glb_ld
      use tdpack
      implicit none
#include <arch_specific.hf>

      logical  F_absvor_L
      integer  F_filtqq, Minx,Maxx,Miny,Maxy,Nk
      real     F_QR (Minx:Maxx,Miny:Maxy,Nk), &
               F_QQ (Minx:Maxx,Miny:Maxy,Nk), &
               F_uu (Minx:Maxx,Miny:Maxy,Nk), &
               F_vv (Minx:Maxx,Miny:Maxy,Nk), F_coefqq

      integer i, j, k, i0, in, j0, jn
      real deg2rad
!
!----------------------------------------------------------------------
!
      i0 = 1
      in = l_niu
      j0 = 1
      jn = l_njv

      do k = 1 , Nk
         do j = j0, jn
            do i = i0, in
               F_QR(i,j,k) = ((F_vv(i+1,j,k) - F_vv(i,j,k)) * geomh_invDXv_8(j)) &
                           - ( (F_uu(i,j+1,k)*geomh_cy_8(j+1) - F_uu(i,j  ,k)*geomh_cy_8(j  )) &
                           * geomh_invDY_8 * geomh_invcyv_8(j))
            end do
         end do
         F_QR(1:i0-1,:,k) = 0. ; F_QR(in+1:l_ni,:,k)= 0.
         F_QR(:,1:j0-1,k) = 0. ; F_QR(:,jn+1:l_nj,k)= 0.
      end do

      if (F_filtqq > 0) call filter ( F_QR, F_filtqq,F_coefqq, &
                                  l_minx,l_maxx,l_miny,l_maxy,Nk )

      if (F_absvor_L)then
         deg2rad= pi_8/180.d0
         do k =  1, Nk
            do j = j0, jn
               do i = i0, in
                  F_QQ(i,j,k)= F_QR(i,j,k) + 2.0*Dcst_omega_8 &
                             * sin(geomh_latrx(i,j)*deg2rad)
               end do
            end do
            F_QQ(1:i0-1,:,k) = 0. ; F_QQ(in+1:l_ni,:,k)= 0.
            F_QQ(:,1:j0-1,k) = 0. ; F_QQ(:,jn+1:l_nj,k)= 0.
         end do
      end if
!
!----------------------------------------------------------------------
!
      return
      end
