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
      subroutine ez_irgdint_1_nw(zo,px,py,npts,ax,ay,z,ni,nj)
      implicit none
      integer npts,ni,nj,wrap,limite
      integer i, j,n, iplus1
      real zo(npts),px(npts),py(npts)
      real ax(ni),ay(nj)
      real z(ni,nj)
      real*8 x, y, x1, x2, y1,y2,dx,dy

#include "zlin8.cdk"
      
      do n=1,npts
         i = min(ni-1,max(1,ifix(px(n))))
         j = min(nj-1,max(1,ifix(py(n))))
         
         x1=ax(i)
         x2=ax(i+1)
         
         x = ax(i) + (x2-x1)*(px(n)-i)
         y = ay(j) + (ay(j+1)-ay(j))*(py(n)-j)
         
         dx = (x - x1)/(x2-x1)
         dy = (y - ay(j))/(ay(j+1)-ay(j))
         
         y1 = zlin(dble(z(i,j)),dble(z(i+1,j)),dx)
         y2 = zlin(dble(z(i,j+1)),dble(z(i+1,j+1)),dx)
         zo(n) = zlin(y1,y2,dy)
      enddo
      
      return
      end
      
