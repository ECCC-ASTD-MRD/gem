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

!*s/r cartall_8 - Computes the F_lon,F_lat positions for a rotated system
!
      subroutine cartall_8 ( F_lon,F_lat,F_xyz_8,n )
      implicit none

      integer n
      real*8 F_lon(*),F_lat(*), F_xyz_8(3,*)
!
!author
!     unknown
!
!revision
! v2_00 - Lee V.            - initial MPI version (from cartall v1_03)
! v2_30 - Dugas B.          - input real*8 cartesian coordinates
!
!object
!     see above ID
!     
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_lon        O   F_longitude positions for a rotated coordinate system
! F_lat        O   F_latitude  positions for a rotated coordinate system
! F_xyz_8      I   rotation matrix
!
!*
      integer i
      real*8, parameter :: ONE = 1.0d0, PI = 180.0d0
      real*8 rad_8
!
!     ---------------------------------------------------------------
!
      rad_8 = PI/acos(-ONE)
      do i=1,n
         F_lat(i) = asin(max(-1.d0,min(1.d0,F_xyz_8(3,i)))) * rad_8
         F_lon(i) = atan2( F_xyz_8(2,i),F_xyz_8(1,i) ) * rad_8
         F_lon(i) = mod( F_lon(i) , 360.0d0 )
         if (F_lon(i).lt.0.0) F_lon(i) = F_lon(i)+360.0
      enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine cartall_8



!*s/r cartall - Computes the F_lon,F_lat positions for a rotated system
!
      subroutine cartall( F_lon,F_lat,F_xyz_8,n)
      implicit none
#include <arch_specific.hf>
      integer n
      real*8 F_xyz_8(3,*)
      real F_lon(*),F_lat(*)
!
!author
!     unknown
!
!revision
! v2_00 - Lee V.            - initial MPI version (from cartall v1_03)
! v2_30 - Dugas B.          - input real*8 cartesian coordinates
!
!object
!     see above ID
!     
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_lon        O   F_longitude positions for a rotated coordinate system
! F_lat        O   F_latitude  positions for a rotated coordinate system
! F_xyz_8      I   rotation matrix
!
!*
      integer i
      real*8, parameter :: ONE = 1.0d0, PI = 180.0d0
      real*8 rad_8
!
!     ---------------------------------------------------------------
!
      rad_8 = PI/acos(-ONE)
      do i=1,n
         F_lat(i) = asin(max(-1.d0,min(1.d0,F_xyz_8(3,i)))) * rad_8
         F_lon(i) = atan2( F_xyz_8(2,i),F_xyz_8(1,i) ) * rad_8
         F_lon(i) = amod( F_lon(i) , 360.0 )
         if (F_lon(i).lt.0.0) F_lon(i) = F_lon(i)+360.0
      enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine cartall
