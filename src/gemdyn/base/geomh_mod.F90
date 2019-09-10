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
module geomh
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save
!______________________________________________________________________
!                                                                      |
!  HORIZONTAL GRID COORDINATES and other related parameters            |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! geomh_x_8          | longitude                                       |
! geomh_xu_8         | longitude of u points                           |
! geomh_y_8          | latitude                                        |
! geomh_yv_8         | latitude  of v points                           |
!--------------------|-------------------------------------------------|
! geomh_hx_8         | distance between grid points in x direction     |
! geomh_hy_8         | distance between grid points in y direction     |
! geomh_inv_hx_8     | inverse of geomh_hx_8                           |
! geomh_inv_hy_8     | inverse of geomh_hy_8                           |
!----------------------------------------------------------------------
! geomh_sx_8         | sine of longitude                               |
! geomh_sy_8         | sine of latitude                                |
! geomh_cx_8         | cosine of longitude                             |
! geomh_cy_8         | cosine of latitude                              |
! geomh_cy2_8        | cosine of latitude squared                      |
! geomh_cyv_8        | cosine of latitude for v points                 |
! geomh_cyv2_8       | cosine of latitude squared for v points         |
! geomh_cyM_8     ** | SPECIAL for MODEL core only                     |
!--------------------|-------------------------------------------------|
! geomh_invcy_8      | inverse of geomh_cy_8                           |
! geomh_invcyv_8     | inverse of geomh_cyv_8                          |
! geomh_invcy2_8     | inverse of geomh_cy2_8                          |
! geomh_invcyv2_8    | inverse of geomh_cyv2_8                         |
!--------------------|-------------------------------------------------|
! geomh_tyoa_8       | tangent of latitude divided by rayt             |
! geomh_tyoav_8      | tangent of latitude divided by rayt (v points)  |
!--------------------|-------------------------------------------------|
! geomh_invDY_8      | inverse DELTAY = cos(theta)2/dsin(theta)        |
! geomh_invDYM_8  ** | SPECIAL for MODEL core only                     |
! geomh_invDYMv_8 ** | SPECIAL for MODEL core only                     |
! geomh_invDX_8      | inverse DELTAX = 1/d(lambda)                    |
! geomh_invDXM_8  ** | SPECIAL for MODEL core only                     |
! geomh_invDXMu_8 ** | SPECIAL for MODEL core only                     |
! geomh_invDXv_8     | inverse DELTAX for v points (ex: for vorticity) |
!----------------------------------------------------------------------
! geomh_area_8       | Area = cos(y)*hx*hy                             |
! geomh_mask_8       | Mask to avoid double count of area in YY overlap|
!----------------------------------------------------------------------
   real(kind=REAL64) :: geomh_invDY_8, geomh_hx_8, geomh_hy_8, geomh_inv_hx_8, geomh_inv_hy_8

   real(kind=REAL64), dimension(:), allocatable :: &
          geomh_x_8    ,geomh_xu_8    ,geomh_y_8     ,geomh_yv_8     ,&
          geomh_sx_8   ,geomh_sy_8    ,geomh_cx_8    ,geomh_cy_8     ,&
          geomh_cy2_8  ,geomh_cyv_8   ,geomh_cyv2_8  ,geomh_cyM_8    ,&
          geomh_invcy_8,geomh_invcyv_8,geomh_invcy2_8,geomh_invcyv2_8,&
          geomh_tyoa_8 ,geomh_tyoav_8 ,&
          geomh_invDYM_8,geomh_invDYMv_8,&
          geomh_invDX_8 ,geomh_invDXM_8,geomh_invDXMu_8,geomh_invDXv_8

   real(kind=REAL64), dimension(:,:), allocatable :: geomh_area_8, geomh_mask_8

!______________________________________________________________________
!                                                                      |
!  LATITUDES AND LONGITUDES IN DEGREES                                 |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! geomh_latgs        | phi latitudes in degrees for model output       |
! geomh_longs        | phi longitudes in degrees for model output      |
! geomh_latgv        | V   latitudes in degrees for model output       |
! geomh_longu        | U   longitudes in degrees for model output      |
! geomh_latrx        | phi geographical latitudes in degrees           |
! geomh_lonrx        | phi geographical longitudes in degrees          |
!----------------------------------------------------------------------
      integer geomh_minx, geomh_maxx, geomh_miny, geomh_maxy
      real, dimension(:  ), pointer, contiguous :: geomh_lonQ, geomh_latQ, &
                                                   geomh_lonF, geomh_latF
      real, dimension(:  ), pointer, contiguous :: geomh_latgs, geomh_longs, &
                                                   geomh_latgv, geomh_longu
      real, dimension(:,:), pointer, contiguous :: geomh_latrx, geomh_lonrx, &
                                                   geomh_latij, geomh_lonij

end module geomh
