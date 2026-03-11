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

!> Horizontal grid coordinates and other related parameters
module geomh
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !> Inverse DELTAY = cos(theta)2/dsin(theta)
   real(kind=REAL64) :: geomh_invDY_8
   !> Distance between grid points in x direction
   real(kind=REAL64) :: geomh_hx_8
   !> Distance between grid points in y direction
   real(kind=REAL64) :: geomh_hy_8
   !> Inverse of geomh_hx_8 (Distance between grid points in x direction)
   real(kind=REAL64) :: geomh_inv_hx_8
   !> Inverse of geomh_hy_8 (Distance between grid points in y direction)
   real(kind=REAL64) :: geomh_inv_hy_8

   !> Longitude
   real(kind=REAL64), dimension(:), allocatable :: geomh_x_8
   !> Longitude of u points
   real(kind=REAL64), dimension(:), allocatable :: geomh_xu_8
   !> Latitude
   real(kind=REAL64), dimension(:), allocatable :: geomh_y_8
   !> Latitude  of v points
   real(kind=REAL64), dimension(:), allocatable :: geomh_yv_8
   !> Sine of longitude
   real(kind=REAL64), dimension(:), allocatable :: geomh_sx_8
   !> Sine of latitude
   real(kind=REAL64), dimension(:), allocatable :: geomh_sy_8
   !> Cosine of longitude
   real(kind=REAL64), dimension(:), allocatable :: geomh_cx_8
   !> Cosine of latitude
   real(kind=REAL64), dimension(:), allocatable :: geomh_cy_8
   !> Cosine of latitude squared
   real(kind=REAL64), dimension(:), allocatable :: geomh_cy2_8
   !> Cosine of latitude for v points
   real(kind=REAL64), dimension(:), allocatable :: geomh_cyv_8
   !> Cosine of latitude squared for v points
   real(kind=REAL64), dimension(:), allocatable :: geomh_cyv2_8
   !> SPECIAL for MODEL core only - Copy of geomh_cyv_8(1-G_haloy:l_nj+G_haloy)
   !> The whole array is never initialized.
   !> The only assignment is in gemdyn/base/set_geomh.F90
   real(kind=REAL64), dimension(:), allocatable :: geomh_cyM_8
   !> Inverse of geomh_cy_8 (Cosine of latitude)
   real(kind=REAL64), dimension(:), allocatable :: geomh_invcy_8
   !> Inverse of geomh_cyv_8 (Cosine of latitude for v points)
   real(kind=REAL64), dimension(:), allocatable :: geomh_invcyv_8
   !> Inverse of geomh_cy2_8 (Cosine of latitude squared)
   real(kind=REAL64), dimension(:), allocatable :: geomh_invcy2_8
   !> Inverse of geomh_cyv2_8 (Cosine of latitude squared for v points)
   real(kind=REAL64), dimension(:), allocatable :: geomh_invcyv2_8
   !> Tangent of latitude divided by rayt
   real(kind=REAL64), dimension(:), allocatable :: geomh_tyoa_8
   !> Tangent of latitude divided by rayt (v points)
   real(kind=REAL64), dimension(:), allocatable :: geomh_tyoav_8
   !> SPECIAL for MODEL core only - geomh_invDY_8*geomh_invcy_8(1-G_haloy:l_nj+G_haloy)
   !> The whole array is never initialized.
   !> The only assignment is in gemdyn/base/set_geomh.F90
   real(kind=REAL64), dimension(:), allocatable :: geomh_invDYM_8
   !> SPECIAL for MODEL core only - Copy of geomh_invDY_8(1-G_haloy:l_nj+G_haloy)
   !> The whole array is never initialized.
   !> The only assignment is in gemdyn/base/set_geomh.F90
   real(kind=REAL64), dimension(:), allocatable :: geomh_invDYMv_8
   !> Inverse DELTAX = 1/d(lambda)
   real(kind=REAL64), dimension(:), allocatable :: geomh_invDX_8
   !> SPECIAL for MODEL core only -
   real(kind=REAL64), dimension(:), allocatable :: geomh_invDXM_8
   !> SPECIAL for MODEL core only -
   real(kind=REAL64), dimension(:), allocatable :: geomh_invDXMu_8
   !> Inverse DELTAX for v points (ex: for vorticity)
   real(kind=REAL64), dimension(:), allocatable :: geomh_invDXv_8

   !> Area = cos(y)*hx*hy
   real(kind=REAL64), dimension(:,:), allocatable :: geomh_area_8
   !> Mask to avoid double count of area in YY overlap
   real(kind=REAL64), dimension(:,:), allocatable :: geomh_mask_8
   !> Area_Mask = Area*Mask
   real(kind=REAL64), dimension(:,:), allocatable :: geomh_area_mask_8

   integer :: geomh_minx
   integer :: geomh_maxx
   integer :: geomh_miny
   integer :: geomh_maxy

   real, dimension(:), pointer, contiguous :: geomh_lonQ
   real, dimension(:), pointer, contiguous :: geomh_latQ
   real, dimension(:), pointer, contiguous :: geomh_lonF
   real, dimension(:), pointer, contiguous :: geomh_latF

   real, dimension(:), pointer, contiguous :: geomh_4output
 
   !> Phi latitudes in degrees for model output
   real, dimension(:), pointer, contiguous :: geomh_latgs
   !> Phi longitudes in degrees for model output
   real, dimension(:), pointer, contiguous :: geomh_longs
   !> V latitudes in degrees for model output
   real, dimension(:), pointer, contiguous :: geomh_latgv
   !> U longitudes in degrees for model output
   real, dimension(:), pointer, contiguous :: geomh_longu

   !> Phi geographical latitudes in degrees
   real, dimension(:,:), pointer, contiguous :: geomh_latrx
   !> Phi geographical longitudes in degrees
   real, dimension(:,:), pointer, contiguous :: geomh_lonrx
   real, dimension(:,:), pointer, contiguous :: geomh_latij
   real, dimension(:,:), pointer, contiguous :: geomh_lonij

end module geomh
