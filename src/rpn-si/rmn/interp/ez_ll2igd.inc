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
!**s/r ez_ll2igd - conversion de coordonnees lat-lon a pts de grille

      subroutine ez_ll2igd(px,py,xlat,xlon,npts,ni,nj,grtyp,grref,      ig1,ig2,ig3,ig4,ax,ay,coordflag)
      implicit none
#include "qqqpar.cdk"

      integer coordflag

      external cigaxg,ez_vxyfll,ez_llll2gd,permut
      
      integer i, npts, ni, nj
      real px(npts),py(npts),xlat(npts),xlon(npts)
      character*1 grtyp, grref
      integer ig1,ig2,ig3,ig4
      real ax(ni), ay(nj)
      real pi,pj,dgrw,d60
      real dlat, dlon, xlat0, xlon0
      real xlat1,xlon1,xlat2,xlon2
      real lonref

      integer indx, indy
      integer ez_cherche
      external ez_cherche
      
      if (grref .eq. 'N') then
         call cigaxg(grref,  PI, PJ, D60, DGRW, ig1, ig2, ig3, ig4)
         call ez_vxyfll(px, py, xlat, xlon, npts, d60,dgrw,pi,pj,1)
      endif

      if (grref .eq. 'S') then
         call cigaxg(grref,  PI, PJ, D60, DGRW, ig1, ig2, ig3, ig4)
         call ez_vxyfll(px, py, xlat, xlon, npts, d60,dgrw,pi,pj,2)
      endif

      if (grref.eq.'L') then
         call cigaxg(grref,xlat0,xlon0,dlat,dlon,ig1,ig2,ig3,ig4)
         if (ax(1).lt.0.0) then
            lonref = -180.0
         else
            lonref = 0.0
         endif
         call ez_llll2gd(px,py,xlat,xlon,npts,xlat0,xlon0,dlat,dlon,lonref)
         do i=1,npts
            px(i) = px(i)-1.0
            py(i) = py(i)-1.0
         enddo
      endif

      if (grref.eq.'E') then
         call cigaxg(grref,xlat1,xlon1,xlat2,xlon2,ig1,ig2,ig3,ig4)
         call ez_gfxyfll(xlon,xlat,px,py,npts,xlat1,xlon1,xlat2,xlon2)
      endif
      
      if (coordflag .eq. relatif) then
         do 10 i=1,npts
            indx = ez_cherche(px(i),ax,ni)
            indy = ez_cherche(py(i),ay,nj)
            
            if (indx .ge. ni) indx = ni - 1
            if (indy .ge. nj) indy = nj - 1
            
            px(i) = real(indx)+(px(i)-ax(indx))/(ax(indx+1)-ax(indx))
            py(i) = real(indy)+(py(i)-ay(indy))/(ay(indy+1)-ay(indy))
 10      continue

      endif

      return
      end
