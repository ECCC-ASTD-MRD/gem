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
      subroutine ez_xtrap(zo, px, py, npts, z, ni,nj,ordint,codxtrap,valxtrap)
      implicit none
      
      integer npts,ni,nj
      real zo(npts),px(npts),py(npts)
      real z(ni,nj)
      integer ordint,codxtrap
      
      integer n, i, j, i1, j1, offl, offr
      real rmin, rmax, tempmin, tempmax, valxtrap
      
#include "ez_qqqxtrp.cdk"

      if (ordint .eq. cubique) then
         offr = 2
         offl = 1
      else
         offr = 0
         offl = 0
      endif

      rmin = z(1,1)
      rmax = z(1,1)

      do j=1,nj
         do i=1,ni
            if(z(i,j).lt.rmin) rmin = z(i,j)
            if(z(i,j).gt.rmax) rmax = z(i,j)
         enddo
      enddo
      
      tempmin = rmin - 0.05*(rmax - rmin)
      tempmax = rmax + 0.05*(rmax - rmin)
      rmin = tempmin
      rmax = tempmax
      
      if (codxtrap.eq.voisin) then
         do n=1,npts
            i = ifix(px(n))
            j = ifix(py(n))
            if (i.lt.1.or.j.lt.1.or.i.gt.ni.or.j.gt.nj) then
               i1 = min(ni, max(1, nint(px(n))))
               j1 = min(nj, max(1, nint(py(n))))
               zo(n) = z(i1,j1)
            endif
         enddo
      endif
      
      if (codxtrap .eq. minimum) then
         do n=1,npts
            zo(n) = rmin
         enddo
      endif
      
      if (codxtrap .eq. maximum) then
         do n=1,npts
            zo(n) = rmax
         enddo
      endif
      
      if (codxtrap .eq. valeur) then
         do n=1,npts
            zo(n) = valxtrap
         enddo
      endif
      
      return
      end
      
