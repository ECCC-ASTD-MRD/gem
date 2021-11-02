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
!**s/r ez_ggll2gd - computes the grid co-ordinates of a point on
!                a gaussian grid
!
      subroutine ez_ggll2gd(x,y,xlat,xlon,npts,ni,nj,hem,lroots)
      implicit none

!*--------------------------------------------------------------------

#include "qqqpar.cdk"

      integer npts, ni, nj
      real x(npts), y(npts), xlat(npts), xlon(npts)
      real lroots(nj)
      real abs,del,epsphi
      integer ez_cherche
      external ez_cherche

      integer i,j,guess,hem,indy
      real tmplat, dellat, dellon, xlat0, xlon0
      
      dellon = 360.0 / real(ni)
      xlon0 = 0.0
      
      do 10 i = 1, npts
         x(i) = (xlon(i) - xlon0)/dellon + 1.0
 10   continue
      
      do i=1,npts
         indy = ez_cherche(xlat(i),lroots,nj)
         if (indy .ge. nj) indy = nj - 1
         
         y(i)= real(indy)+(xlat(i)-lroots(indy))/         (lroots(indy+1)-lroots(indy))
      enddo
      return
      end
!-----------------------------------------------------------------
