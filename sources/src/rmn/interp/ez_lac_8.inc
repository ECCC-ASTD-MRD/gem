!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!**s/r ez_lac transformation from a set of points (lat,lon) in the spherical 
!             coordinate system to cartesian space
!
      subroutine ez_lac_8( xyz,lon,lat,n)
      implicit none
      integer n
      real*8 xyz(3,n)
      real lon(n),lat(n) 
!
!author michel roch - april 90
!
!revision
!       001 - michel roch - feb 94 - repackaging for integration in rmnlib
!
!arguments
!    out    xyz     - coordinates in cartesian space
!    in     lon     - longitudes in spherical coordinates
!           lat     - latitudes  in spherical coordinates
!           n       - dimension of the fields
!
!*
      
      integer i
      real dar, cosdar
      dar = acos(-1.)/180.
!     
      do 30 i=1,n
!!!!JP
         cosdar=cos(dar*lat(i))
         xyz(1,i)=cosdar*cos(dar*lon(i))
         xyz(2,i)=cosdar*sin(dar*lon(i))
         xyz(3,i)=sin(dar*lat(i))
 30   continue
!     
      return
      end 
      
      
