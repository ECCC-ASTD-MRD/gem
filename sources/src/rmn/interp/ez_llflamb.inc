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
      subroutine ez_llflamb(xlat,xlon,x,y,npts,grtyp,ig1,ig2,ig3,ig4)
      implicit none
      integer npts
      real xlat(npts),xlon(npts)
      real x(npts),y(npts) 
      character*1 grtyp,gtypout
      integer ig1,ig2,ig3,ig4,nxg

      real xg(15)

      real xlat11,xlon11,xlatninj,xlonninj,dx,dy,latin1,latin2
      real yaxislat,yaxislon
      real x11,y11
      integer i
      
      nxg =15
      if (grtyp.eq.'!') then
         call igaxg95(gtypout,xg,15,grtyp,ig1,ig2,ig3,ig4)
         if (gtypout.eq.'H') then
            xlat11 =   xg( 1)
            xlon11 =   xg( 2)
            xlatninj = xg(10)
            xlonninj = xg(11)
            yaxislon = xg(5)
            yaxislat = 0.5 * (xlat11 + xlatninj)
            latin1   = xg(6)
            latin2   = xg(7)
            dx       = xg(3)*1000.0
            dy       = xg(4)*1000.0

            call ez_lambxyfll99(x11,y11,xlat11,xlon11,1,            latin1,latin2,yaxislat,yaxislon)
            
            do i=1,npts
               x(i) = x11 + dx * (x(i) - 1.0)
               y(i) = y11 + dy * (y(i) - 1.0)
            enddo

            call ez_lambllfxy99(xlat,xlon,x,y,npts,            latin1,latin2,yaxislat,yaxislon)
            
         endif
      endif

      return
      end
      
