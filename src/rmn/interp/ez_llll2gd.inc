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
!**s/r ez_llll2gd - computes the grid co-ordinates of a point
!
      subroutine ez_llll2gd(x,y,dlat,dlon,npts,xlat0,xlon0, dellat, dellon, lonref)
      implicit none

!*--------------------------------------------------------------------

      integer npts
      real x(npts), y(npts), dlat(npts), dlon(npts)
      real xlat0, xlon0, dellat, dellon, lonref
      integer i

      if (lonref.eq.-180.0) then
         do i=1,npts
            if (dlon(i).gt.180.0) then
               dlon(i) = dlon(i) - 360.0
            endif
         enddo
      else
!    mai 2015 -correction by JP Gauthier
         do i=1,npts
            if (dlon(i).lt.0.0) then
               dlon(i) = dlon(i) + 360.0
            endif
         enddo
      endif

      do 10 i=1,npts
         x(i) = (dlon(i) - xlon0)/dellon + 1.0
         y(i) = (dlat(i) - xlat0)/dellat + 1.0
 10   continue

      return
      end
!-----------------------------------------------------------------

