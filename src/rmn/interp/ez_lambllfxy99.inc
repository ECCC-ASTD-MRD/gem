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
      subroutine ez_lambllfxy99(lat,lon,x,y,n,      latin1,latin2,yaxislat,yaxislon)
      implicit none
      integer n
      real x(n),y(n),lat(n),lon(n)
      real lat11,lon11,latninj,lonninj,latin1,latin2,yaxislat,yaxislon
      
      real pi,pisur4,d2r,rn,rphi1,rphi2,r,f,rtan,rho,theta,rhozero,dlon
      integer i
      
      pisur4 = atan(1.0)
      pi     = 4.0 * pisur4
      d2r    = pi / 180.0
      r = 6370997.0
      
      rphi1 = d2r * latin1
      rphi2 = d2r * latin2

      if (rphi1.eq.rphi2) then
         rn = sin(rphi1)
      else
         rn = log(cos(rphi1)/cos(rphi2))
         rn = rn / log((tan(pisur4+0.5*rphi2))/tan(pisur4+0.5*rphi1))
      endif

      rtan = tan(pisur4 + rphi1 * 0.5)
      f = (cos(rphi1)*(rtan ** rn))/rn
      rhozero = r * f / (tan(pisur4+yaxislat*d2r*.5)**rn)

      do i = 1,n
         rho = sign(1.0,rn) * sqrt(x(i)*x(i)+((rhozero-y(i))**2))
         theta = atan(x(i)/(rhozero - y(i)))

         lat(i) = (2.0 * atan((r*f/rho)**(1.0/rn)) - 0.5*pi)/d2r
         lon(i) = theta/(d2r*rn) + yaxislon
      enddo
      
      return
      end
