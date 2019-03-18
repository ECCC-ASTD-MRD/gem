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
!**s/r ez_vxyfll - computes the grid co-ordinates of a point
!
      subroutine ez_vxyfll(x,y,dlat,dlon,npts,d60,dgrw,pi,pj,nhem)
      implicit none

!author   - j.d. henderson  -  feb 75
!        
!revision 001   c. thibeault  -  nov 79  documentation
!revision 002   y. chartier   -  feb 91  vectorization
!                                                     
!language - fortran
!                  
!object(vxyfll)
!     - computes the grid co-ordinates measured from the pole of a
!       point, given the latitude and longitude in degrees.
!
!libraries
!         - source  rmnsourcelib,id=rmnp     deck=vxyfll
!         - object  rmnlib,id=rmnp
!
!usage    - call vxyfll(x,y,dlat,dlon,npts,d60,dgrw,pi,pj,nhem)
!
!arguments
!   out   - x    - x-co-ordinate of the point as measured with pole
!                  as origin
!         - y    - y-co-ordinate of the point as measured with pole
!                  as origin
!   in    - dlat - latitude in degrees (-90 to +90, positive n)
!         - dlon - longitude in degrees (-180 to +180, positive e)
!         - d60  - grid length (in metres) of the polar stereographic
!                  grid at 60 degrees
!         - dgrw - orientation of greenwich meridian with respect to
!                  the grid (in degrees)
!         - nhem - 1 for northern hemisphere
!                  2 for southern hemisphere
!                                                                      

!notes    - the companion routine llfxy, which computes the latitude
!         - and longitude given the grid-coordinates, 
!         - is also available.
!*--------------------------------------------------------------------
#include "pi.cdk"
#include "qqqpar.cdk"

      integer npts, nhem
      real x(npts), y(npts), dlat(npts), dlon(npts)
      real d60, dgrw, pi, pj
      real*8 re,rlon,rlat,sinlat,r
      integer i
      
      re=1.866025d0*6.371d+6/d60
!     
      if (nhem.eq.NORD) then
         do 10 i=1,npts
            rlon=dgtord*(dlon(i)+dgrw)
            rlat=dgtord*dlat(i)
            sinlat=sin(rlat)
            r=re*sqrt((1.d0-sinlat)/(1.d0+sinlat))
            x(i)=r*cos(rlon) + pi
            y(i)=r*sin(rlon) + pj
 10      continue
         return
      endif

      if (nhem.eq.SUD) then
         do 20 i=1,npts
            rlon = dlon(i)
            if (rlon.gt.180.0d0) rlon = rlon - 360.0d0
            rlon=dgtord*(-rlon+dgrw)
            rlat=dgtord*(-dlat(i))
            sinlat=sin(rlat)
            r=re*sqrt((1.d0-sinlat)/(1.d0+sinlat))
            x(i)=r*cos(rlon)+pi
            y(i)=r*sin(rlon)+pj
 20      continue
      endif
      
      return
      end
!-----------------------------------------------------------------
