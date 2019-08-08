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
      subroutine ez_irgdint_1_w(zo,px,py,npts,ax,ay,z,ni,j1,j2,wrap)
      implicit none
      integer npts,ni,nj,wrap,limite
      integer i, j,j1,j2,n, iplus1
      real zo(npts),px(npts),py(npts)
      real ax(ni),ay(j1:j2)
      real z(ni,j1:j2)
      
      real*8 x, y, x1, x2, y1,y2,dx,dy

#include "zlin8.cdk"
      do n=1,npts
         limite = ni+2-wrap
         i = min(ni-2+wrap,max(1,ifix(px(n))))
         j = min(j2-1,max(j1+1, ifix(py(n))))
         
         if (j.lt.0) then
            j = j-1
         endif
         
         iplus1 = i+1
         
         x1=ax(i)
         if (iplus1.le.ni) then
            x2 = ax(iplus1)
         endif

         if (wrap.gt.0.and.(i.eq.(ni-2+wrap))) then
            iplus1  = mod(limite+i+1,limite)
            x2=ax(2)+ax(ni)
         endif
         
         x = x1 + (x2-x1)*(px(n)-i)
         y = ay(j) + (ay(j+1)-ay(j))*(py(n)-j)
         
         dx = (x - x1)/(x2-x1)
         dy = (y - ay(j))/(ay(j+1)-ay(j))
         
         y1 = zlin(dble(z(i,j)),dble(z(iplus1,j)),dx)
         y2 = zlin(dble(z(i,j+1)),dble(z(iplus1,j+1)),dx)
         zo(n) = zlin(y1,y2,dy)
      enddo
      
      return
      end
      
