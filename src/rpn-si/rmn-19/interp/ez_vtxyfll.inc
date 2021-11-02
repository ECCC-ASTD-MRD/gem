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
      subroutine testprojt

      implicit none

      integer i,j,ni,nj
      real clat, clon, d60, dgrw
      real lat(11,11), lon(11,11)

      real tlat(5,5), tlon(5,5)
      real tlat2(5,5), tlon2(5,5)
      real x(5,5),y(5,5)

      ni = 11
      nj = 11
      clat = 45.0
      clon = 260.0
      d60 = 100000.0
      dgrw = 0.0

      do j=1,5
         do i=1,5
            tlat(i,j) = 35.0 + (j-1)*5
            tlon(i,j) = 250.0 +(i-1)*5
         enddo
      enddo

      call ez_vtxyfll(x, y, tlat, tlon, clat, clon, d60, dgrw, ni, nj, 25)
      
      do j=1,5
         do i=1,5
            print *, i,j,tlat(i,j),tlon(i,j),x(i,j),y(i,j)
         enddo
      enddo
      
      print *, '***************************************************'

      call ez_vtllfxy(tlat2, tlon2, x, y, clat, clon, d60, dgrw, ni, nj, 25)
      
      do j=1,5
         do i=1,5
            print *, i,j,tlat2(i,j)-tlat(i,j),tlon2(i,j) - tlon(i,j)
         enddo
      enddo
      
      stop
      end


      subroutine ez_vtxyfll(x, y, lat, lon, clat, clon, d60, dgrw, ni, nj, n)

      implicit none
#include "pi.cdk"
      
      integer i, n, ni, nj
      real x(n), y(n), lat(n), lon(n), clat, clon, d60, dgrw
      real r,k,lat0,lon0
      real offsetx, offsety,sinclat,cosclat,sinclon,cosclon

      r = 6371000.0 
      sinclat = sin (clat * dgtord)
      cosclat = cos (clat * dgtord)
      sinclon = sin (clon * dgtord)
      cosclon = cos (clon * dgtord)

      offsetx = (-ni-1) * 0.5
      offsety = (-nj-1) * 0.5

      do i=1,n
         k = 2.0 / (1.0 + sinclat*sin(lat(i) * dgtord)+ cosclat* cos(lat(i)*dgtord)*cos(dgtord*(lon(i)-clon)))
         x(i) = r * k * cos(lat(i)*dgtord) * sin(dgtord*(lon(i)-clon))
         y(i) = r * k * (cosclat*sin(lat(i) * dgtord) -  (sinclat * cos(lat(i)*dgtord)*cos(dgtord*(lon(i)-clon))))
         x(i) = x(i) / d60 - offsetx
         y(i) = y(i) / d60 - offsety
      enddo

      return
      end

      subroutine ez_vtllfxy(lat, lon, x, y, clat, clon, d60, dgrw, ni, nj, n)

      implicit none
#include "pi.cdk"
      
      integer i, n, ni, nj
      real x(n), y(n), lat(n), lon(n), clat, clon, d60, dgrw
      real r,k,lat0,lon0
      real offsetx, offsety,sinclat,cosclat,sinclon,cosclon
      real rho, c, a, b, temp


      r = 6371000.0 
      sinclat = sin (clat * dgtord)
      cosclat = cos (clat * dgtord)
      sinclon = sin (clon * dgtord)
      cosclon = cos (clon * dgtord)

      offsetx = (-ni-1) * 0.5
      offsety = (-nj-1) * 0.5

      do i=1,n
         x(i) = (x(i) + offsetx) * d60
         y(i) = (y(i) + offsety) * d60

         rho = sqrt(x(i)*x(i) + y(i)*y(i))
         if (rho.eq.0.0) then
            lat(i) = clat
            lon(i) = clon
         else
            c = 2.0 * atan(rho/(2.0*r))
            
            temp = cos(c)*sinclat + (y(i) * sin(c) * cosclat)/rho
            lat(i) = rdtodg * asin(max(-1.0, min(temp, 1.0)))
            lon(i) = clon + rdtodg * atan(((x(i)*sin(c))/            (rho*cosclat*cos(c)-y(i)*sinclat*sin(c))))
            a = x(i) * sin(c)
            b = rho*cosclat*cos(c)-y(i)*sinclat*sin(c)
            lon(i) = clon + rdtodg * atan2(a,b)

            
         endif
         lon(i) = mod(mod(lon(i),360.0)+360.0,360.0)
      enddo
      
      return
      end
