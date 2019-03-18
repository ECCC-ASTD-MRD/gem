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
!**s/r ez_crot - compute the rotation matrix that allows the transformation
!             from the spherical coordinate system to the rotated spherical
!             coordinate system.
!
      subroutine ez_crot( r, ri, lon1, lat1, lon2, lat2 )
      implicit none
      real r(3,3), ri(3,3)
      real lon1,lat1,lon2,lat2
!
!author michel roch - april 1990
!
!revision
!	001 - michel roch - feb 94 - repackaging for integration in rmnlib
!
!arguments
!    out    r       - rotation matrix from true to rotated coordinate system
!           ri      - rotation matrix from rotated to true coordinate system
!    in     lon1    - longitude on the unrotated coordinate system corresponding to
!                     the point (lat,lon)=(0,180) of the rotated coordinate system
!           lat1    - latitude on the unrotated coordinate system corresponding to
!                     the point (lat,lon)=(0,180) of the rotated coordinate system
!           lon2    - longitude on the unrotated coordinate system corresponding to
!                     a point (lat,lon) located on the equator of the rotated 
!                     coordinate system
!           lat2    - latitude on the unrotated coordinate system corresponding to
!                     a point (lat,lon) located on the equator of the rotated 
!                     coordinate system
!
!
!modules
      external ez_lac
!
!*
      integer i,j
      real a,b,c,d
      real xyz1(3),xyz2(3)
!
!     calcul des coordonnees cartesiennes des points r1 et r2
!
      call ez_lac(xyz1,lon1,lat1,1)
      call ez_lac(xyz2,lon2,lat2,1)
!
!     calcul de a=cos(alpha)
!
      a=(xyz1(1)*xyz2(1)) + (xyz1(2)*xyz2(2)) + (xyz1(3)*xyz2(3))
!
!     calcul de b=sin(alpha)
!
      b=sqrt(  ( (xyz1(2)*xyz2(3)) - (xyz2(2)*xyz1(3)) )**2 +          ( (xyz2(1)*xyz1(3)) - (xyz1(1)*xyz2(3)) )**2 +          ( (xyz1(1)*xyz2(2)) - (xyz2(1)*xyz1(2)) )**2  )
!
!     calcul de c=norm(-r1)
!
      c=sqrt( xyz1(1)**2 + xyz1(2)**2 + xyz1(3)**2 )
!
!     calcul de d=norm(r4)
!
      d=sqrt(  ( ( (a*xyz1(1)) - xyz2(1) ) / b )**2 +          ( ( (a*xyz1(2)) - xyz2(2) ) / b )**2 +          ( ( (a*xyz1(3)) - xyz2(3) ) / b )**2  )
!
      r(1,1)=     -xyz1(1)/c
      r(1,2)=     -xyz1(2)/c
      r(1,3)=     -xyz1(3)/c
      r(2,1)=  ( ((a*xyz1(1)) - xyz2(1)) / b)/d
      r(2,2)=  ( ((a*xyz1(2)) - xyz2(2)) / b)/d
      r(2,3)=  ( ((a*xyz1(3)) - xyz2(3)) / b)/d
      r(3,1)=  ( ( xyz1(2)*xyz2(3)) - (xyz2(2)*xyz1(3)) ) / b
      r(3,2)=  ( ( xyz2(1)*xyz1(3)) - (xyz1(1)*xyz2(3)) ) / b
      r(3,3)=  ( ( xyz1(1)*xyz2(2)) - (xyz2(1)*xyz1(2)) ) / b
!
!     matrice de rotation inverse
!
      do 10 i=1,3
      do 10 j=1,3
 10      ri(i,j)=r(j,i)
!
      return
      end 
