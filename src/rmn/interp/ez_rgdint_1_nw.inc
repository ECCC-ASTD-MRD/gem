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
      subroutine ez_rgdint_1_nw(zo,px,py,npts,z,ni,j1,j2)
      implicit none

      integer npts,ni,j1,j2,i,j,n
      real zo(npts),px(npts),py(npts)
      real z(ni,j1:j2)
      real*8 dx,dy,y2,y3

#include "zlin8.cdk"

      do n=1,npts
         i = min(ni-1,max(1,ifix(px(n))))
         j = min(j2-1,max(j1,ifix(py(n))))
         
         dx = px(n) - i
         dy = py(n) - j
         
         y2=zlin(dble(z(i,j  )),dble(z(i+1,j  )),dx)
         y3=zlin(dble(z(i,j+1)),dble(z(i+1,j+1)),dx)
         
         zo(n)=zlin(y2,y3,dy)
      enddo

      return 
      end
