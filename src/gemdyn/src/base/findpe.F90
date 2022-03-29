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
!*** function find_col
!
      integer function find_col (F_posx)
      use glb_ld
      use adz_mem
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN) ::  F_posx

      integer pe,ipos
!
!     ---------------------------------------------------------------
!
      find_col=-1
      if (F_posx<Adz_glbminpos) return      
      ipos=floor(F_posx)
      do pe=0,Ptopo_npex
         find_col= pe
         if (Adz_gindx_alongX(pe)>=ipos) exit
      end do
!
!     ---------------------------------------------------------------
!
      return
      end function find_col
      
      integer function find_row (F_posy)
      use glb_ld
      use adz_mem
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN) ::  F_posy

      integer pe,ipos
!
!     ---------------------------------------------------------------
!
      find_row=-1
      if (F_posy<Adz_glbminpos) return      
      ipos=floor(F_posy)
      do pe=0,Ptopo_npey
         find_row= pe
         if (Adz_gindx_alongY(pe)>=ipos) exit
      end do
!
!     ---------------------------------------------------------------
!
      return
      end function find_row

!***  function findpeyy - to find PE given global I,J index in other GRID
      
      integer function findpeyy(F_x_8,F_y_8,F_xa_8,F_ya_8)
      use yyg_param
      use glb_ld
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN ) :: F_x_8,F_y_8
      real(kind=REAL64), intent(OUT) :: F_xa_8,F_ya_8

!author V.Lee - November 2020
!objective Get point on other grid (couleur)

      integer, external :: find_col,find_row
      real(kind=REAL64) :: x_d, y_d, s(2,2), x_a, y_a
      integer ecol, erow
!
!     ---------------------------------------------------------------
!
! Convert index positions to radian for smat
      x_d = YYG_xg_8(1) + (F_x_8-1.0d0)*geomh_hx_8 - pi_8
      y_d = YYG_yg_8(1) + (F_y_8-1.0d0)*geomh_hy_8
      call smat(s,x_a,y_a,x_d,y_d)
      x_a =x_a+pi_8
      
! Convert new positions from radian to global index
      F_xa_8 = (x_a-YYG_xg_8(1))/geomh_hx_8 + 1.0d0
      F_ya_8 = (y_a-YYG_yg_8(1))/geomh_hy_8 + 1.0d0

      ecol= find_col(F_xa_8)
      erow= find_row(F_ya_8)
      findpeyy= Ptopo_colrow(Ptopo_couleur,ecol,erow) + (1-Ptopo_couleur)*Ptopo_numproc
!
!     ---------------------------------------------------------------
!
      return
      end function findpeyy

