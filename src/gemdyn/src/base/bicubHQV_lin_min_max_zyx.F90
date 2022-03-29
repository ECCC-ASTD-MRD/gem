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

!**s/r bicubHQV_lin_min_max_zyx - SL Bi-Cubic H interpolation and Quintic V interpolation except 
!                                 Cubic or Linear V interpolation near vertical boundaries
!                                 (Based on tricublin_beta.c): Interpolation ZYX 

      subroutine bicubHQV_lin_min_max_zyx (xyz, li, mi, ma, src, cx, cy, qz, lx, ly, lz, & 
                                           ni, ninj, zqutic_L, zcubic_L)

      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: ni, ninj
      real, dimension(0:*), intent(in) :: src
      real(kind=REAL64), dimension(0:3), intent(in) :: cx,cy
      real(kind=REAL64), dimension(0:5), intent(in) :: qz 
      real(kind=REAL64), dimension(1:2), intent(in) :: lx,ly,lz 
      logical :: zqutic_L, zcubic_L 
      real, intent(out) :: xyz,li,mi,ma

      integer :: i, ni2, ni3, ninj1, ninj2, ninj3, ninj4, ninj5, ninjc, ninjl
      real(kind=REAL64), dimension(0:3) :: va4, vb4, vc4, vd4
      real, dimension(0:3) :: dst,dsl
      real :: mi_,ma_
!
!---------------------------------------------------------------------
!
      ninjc = 0
      if (zqutic_L) ninjc = ninj 

      ninjl = 0
      if (zqutic_L.or.zcubic_L) ninjl = ninj 

      ninj1 =     0 + ninjc
      ninj2 = ninj1 + ninjl
      ninj3 = ninj2 + ninj
      ninj4 = ninj3 + ninjl
      ninj5 = ninj4 + ninjc
      ni2 = ni  + ni
      ni3 = ni2 + ni

      !Bi-Cubic H interpolation and Quintic V interpolation except 
      !Cubic or Linear V interpolation near vertical boundaries
      !----------------------------------------------------------- 
      do i = 0,3

         va4(i) = src(i    )*qz(0) + src(i    +ninj1)*qz(1) +  src(i    +ninj2)*qz(2) + src(i    +ninj3)*qz(3) + src(i    +ninj4)*qz(4) + src(i    +ninj5)*qz(5)
         vb4(i) = src(i+ni )*qz(0) + src(i+ni +ninj1)*qz(1) +  src(i+ni +ninj2)*qz(2) + src(i+ni +ninj3)*qz(3) + src(i+ni +ninj4)*qz(4) + src(i+ni +ninj5)*qz(5)
         vc4(i) = src(i+ni2)*qz(0) + src(i+ni2+ninj1)*qz(1) +  src(i+ni2+ninj2)*qz(2) + src(i+ni2+ninj3)*qz(3) + src(i+ni2+ninj4)*qz(4) + src(i+ni2+ninj5)*qz(5)
         vd4(i) = src(i+ni3)*qz(0) + src(i+ni3+ninj1)*qz(1) +  src(i+ni3+ninj2)*qz(2) + src(i+ni3+ninj3)*qz(3) + src(i+ni3+ninj4)*qz(4) + src(i+ni3+ninj5)*qz(5)

         dst(i) = va4(i)*cy(0) + vb4(i)*cy(1) + vc4(i)*cy(2) + vd4(i)*cy(3)

      end do

      xyz = dst(0)*cx(0) + dst(1)*cx(1) + dst(2)*cx(2) + dst(3)*cx(3)

      !Linear interpolation over 2x2x2 inner box
      !-----------------------------------------
      do i = 1,2

         vb4(i) = src(i+ni +ninj2)*lz(1) + src(i+ni +ninj3)*lz(2)
         vc4(i) = src(i+ni2+ninj2)*lz(1) + src(i+ni2+ninj3)*lz(2)

         dsl(i) = vb4(i)*ly(1) + vc4(i)*ly(2)

      end do

      li = dsl(1)*lx(1) + dsl(2)*lx(2)

      !MIN/MAX over 2x2x2 inner box
      !---------------------------- 
      ma_ = src(1+ni +ninj2); mi_ = ma_

      do i = 1,2
         ma_ = max(ma_ , src(i+ni +ninj2)); ma_ = max(ma_ , src(i+ni +ninj3));
         mi_ = min(mi_ , src(i+ni +ninj2)); mi_ = min(mi_ , src(i+ni +ninj3));

         ma_ = max(ma_ , src(i+ni2+ninj2)); ma_ = max(ma_ , src(i+ni2+ninj3));
         mi_ = min(mi_ , src(i+ni2+ninj2)); mi_ = min(mi_ , src(i+ni2+ninj3));
      end do

      ma = ma_ ; mi = mi_

      return
      end
