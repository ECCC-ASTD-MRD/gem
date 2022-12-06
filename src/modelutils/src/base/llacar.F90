!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

subroutine llacar(F_xyz_8,F_lon,F_lat,ni,nj)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@object transformation from a set of points (F_lat,F_lon) in
   !        the spherical coordinate system to cartesian space
   !@arguments
   ! F_xyz_8      O    - coordinates in cartesian space
   ! F_lon        I    - longitudes of the grid in spherical coordinates
   ! F_lat        I    - latitudes of the grid in spherical coordinates
   integer ni,nj
   real(REAL64) :: F_xyz_8(3,ni*nj)
   real F_lon(ni),F_lat(nj)
   !@author Michel Roch - April 90
   !@revision
   ! v2_00 - Lee V.            - initial MPI version (from llacar v1_03)
   ! v2_30 - Dugas B.          - output real(REAL64) :: cartesian coordinates

   integer i,j,k
   real(REAL64), parameter :: ONE = 1.0d0, PI = 180.0d0
   real(REAL64) :: deg2rad_8
   !----------------------------------------------------------------
   deg2rad_8 = acos(-ONE)/PI

   k=0
   do j=1,nj
      do i=1,ni
         k=k+1
         F_xyz_8(1,K) = cos(deg2rad_8*F_lat(j))*cos(deg2rad_8*F_lon(i))
         F_xyz_8(2,K) = cos(deg2rad_8*F_lat(j))*sin(deg2rad_8*F_lon(i))
         F_xyz_8(3,K) = sin(deg2rad_8*F_lat(j))
      enddo
   enddo
   !----------------------------------------------------------------
   return
end subroutine llacar


subroutine llacar_8(F_xyz_8,F_lon,F_lat,ni,nj )
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@object transformation from a set of points (F_lat,F_lon) in
   !        the spherical coordinate system to cartesian space
   !@arguments
   ! F_xyz_8      O    - coordinates in cartesian space
   ! F_lon        I    - longitudes of the grid in spherical coordinates
   ! F_lat        I    - latitudes of the grid in spherical coordinates
   integer ni,nj
   real(REAL64) :: F_xyz_8(3,ni*nj),F_lon(ni),F_lat(nj)
   !@author Michel Roch - April 90
   !@revision
   ! v2_00 - Lee V.            - initial MPI version (from llacar v1_03)
   ! v2_30 - Dugas B.          - output real(REAL64) :: cartesian coordinates

   integer i,j,k
   real(REAL64), parameter :: ONE = 1.0d0, PI = 180.0d0
   real(REAL64) ::  deg2rad_8
   !----------------------------------------------------------------
   deg2rad_8 = acos(-ONE)/PI

   k=0
   do j=1,nj
      do i=1,ni
         k=k+1
         F_xyz_8(1,K) = cos(deg2rad_8*F_lat(j))*cos(deg2rad_8*F_lon(i))
         F_xyz_8(2,K) = cos(deg2rad_8*F_lat(j))*sin(deg2rad_8*F_lon(i))
         F_xyz_8(3,K) = sin(deg2rad_8*F_lat(j))
      enddo
   enddo
   !----------------------------------------------------------------
   return
end subroutine llacar_8
