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
!
      subroutine adv_setgrid ( )
      use HORgrid_options
      use glb_ld
      use cstv
      use ver
      use adv_grid
      use outgrid
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

!     @objective:  set parmaters of the advection grid ( from  adx_set )

#include "constants.h"

      integer  i, j
      real(kind=REAL64) :: prhxmn, prhymn
!
!     ---------------------------------------------------------------
!
      adv_maxcfl = max(1,Grd_maxcfl)
      adv_halox = max(1,adv_maxcfl + 1)
      adv_haloy = adv_halox

      adv_int_i_off = l_i0 - 1
      adv_int_j_off = l_j0 - 1

      adv_gminx = 1 - adv_halox
      adv_gmaxx = G_ni + adv_halox
      adv_gminy = 1 - adv_haloy
      adv_gmaxy = G_nj + adv_haloy

      adv_lminx = 1 - adv_halox
      adv_lmaxx = l_ni + adv_halox
      adv_lminy = 1 - adv_haloy
      adv_lmaxy = l_nj + adv_haloy

      adv_iimax = G_ni+2*adv_halox-2
      adv_jjmax = G_nj+adv_haloy
      adv_nit = adv_lmaxx - adv_lminx + 1
      adv_njt = adv_lmaxy - adv_lminy + 1
      adv_nijag = adv_nit * adv_njt

      allocate ( adv_xg_8(adv_gminx:adv_gmaxx), &
                 adv_yg_8(adv_gminy:adv_gmaxy), &
                 adv_xx_8(adv_lminx:adv_lmaxx), &
                 adv_yy_8(adv_lminy:adv_lmaxy), &
                 adv_cy_8(l_nj) )

      do i = 1,G_ni
         adv_xg_8(i) = G_xg_8(i)
      end do

      do j = 1,G_nj
         adv_yg_8(j) = G_yg_8(j)
      end do

      prhxmn =  adv_xg_8(2)-adv_xg_8(1)
      do i = 0,adv_gminx,-1
         adv_xg_8(i) = adv_xg_8(i+1)  - prhxmn
      end do
      do i = G_ni+1,adv_gmaxx
         adv_xg_8(i) = adv_xg_8(i-1) + prhxmn
      end do

      prhymn =  adv_yg_8(2)-adv_yg_8(1)
      do j = 0,adv_gminy,-1
         adv_yg_8(j) = adv_yg_8(j+1) - prhymn
      end do
      do j = G_nj+1,adv_gmaxy
         adv_yg_8(j) = adv_yg_8(j-1) + prhymn
      end do

!- advection grid
      do i = adv_lminx,adv_lmaxx
         adv_xx_8(i) = adv_xg_8(l_i0-1+i)
      end do
      do j = adv_lminy,adv_lmaxy
         adv_yy_8(j) = adv_yg_8(l_j0-1+j)
      end do

      do j = 1,l_nj
         adv_cy_8(j) = cos(adv_yy_8(j))
      end do
!
!---------------------------------------------------------------------
!
      return
      end subroutine adv_setgrid
