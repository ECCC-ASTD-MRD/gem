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

module adv_grid
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
! GRIDS USED IN THE ADVECTION CONTEXT                                  |
!______________________________________________________________________|
!                                                                      |
! 3 different grids are refered to throughout the advection process:   |
!                                                                      |
! GLOBAL GRID:    Global domain                                        |
!                                                                      |
! ADVECTION GRID: In the general case, the tiles adjacent to the poles |
!                 would have an advection source grid periodic in x.   |
!                 For the other tiles, the advection grid would be the |
!                 same as the local grid.                              |
!                                                                      |
! LOCAL GRID:     Local domain                                         |
!                                                                      |
! All of these grids could include halos or not                        |
!______________________________________________________________________|

   !- Global grid
   integer :: adv_maxcfl             !Max CFL defined by user
   integer :: adv_gminx, adv_gmaxx   !min/max x dims for model global fields
   integer :: adv_gminy, adv_gmaxy   !min/max y dims for model global fields

   !- Advection grid
   integer :: adv_halox, adv_haloy   !Advection x/y halo size
   integer :: adv_lminx, adv_lmaxx , adv_lminy, adv_lmaxy !min/max on x/y dims for adv fields
   integer :: adv_nit, adv_njt, adv_nijag, adv_int_i_off, adv_int_j_off
   integer :: adv_iimax,adv_jjmax

   real*8, dimension(:), allocatable :: adv_xg_8, adv_yg_8 !glb grid x,y coor
   real*8, dimension(:), allocatable :: adv_xx_8, adv_yy_8 !local grid x,y coor
   real*8, dimension(:), allocatable :: adv_cy_8           ! = cos(adv_yy_8)

end module adv_grid
