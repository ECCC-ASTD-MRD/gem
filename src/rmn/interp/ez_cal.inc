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
!**s/r qqqcal - compute spherical coordinates from cartesian coordinates
!
      subroutine ez_cal( lon,lat,xyz,n)
      implicit none
      integer n
      real xyz(3,n),lon(n),lat(n)
!
!author michel roch - april 1990
!
!revision
!       001 - michel roch - feb 94 - repackaging for integration in rmnlib
!
!*
!
      integer i
      real rad
!
      rad = 180./acos(-1.00)
!
         do 30 i=1,n
            lat(i) =asin(max(-1.00,min(1.0,xyz(3,i)))) * rad
            lon(i) =atan2( xyz(2,i),xyz(1,i) ) * rad
            lon(i) =amod( lon(i) , 360.0 )
            if(lon(i).lt.0.0) lon(i)=lon(i)+360.0
 30      continue
!
      return
      end 
