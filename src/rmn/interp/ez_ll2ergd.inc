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
!**s/r ll2ergd - conversion de coordonnees lat-lon a pts de grille

      subroutine ez_ll2ergd(px,py,xlat,xlon,npts,ni,nj,grtyp,      ig1,ig2,ig3,ig4)
      implicit none
#include "qqqpar.cdk"
#include "pi.cdk"

      external cigaxg,ez_llll2gd,permut

      integer npts, ni, nj
      real px(npts),py(npts),xlat(npts),xlon(npts)
      character*1 grtyp
      integer ig1,ig2,ig3,ig4
      real dellat, dellon, xlat0, xlon0
      real xlat1,xlat2,xlon1,xlon2
      real xlatgf(npts),xlongf(npts)
      integer ier


      integer i,j
      real rinc

      call cigaxg(grtyp,xlat1,xlon1,xlat2,xlon2,ig1,ig2,ig3,ig4)
      call ez_gfxyfll(xlon,xlat,xlongf,xlatgf,npts,      xlat1,xlon1,xlat2,xlon2)

      if (grtyp.eq.'E') then
         dellon = 360.0 / real(ni-1)
         xlon0 = 0.0

         dellat = 180.0 / real(nj)
         xlat0 = -90. + 0.5*dellat

         call ez_llll2gd(px,py,xlatgf,xlongf,npts,         xlat0,xlon0,dellat,dellon, 0.0)
      endif

      return
      end
!-----------------------------------------------------------------

