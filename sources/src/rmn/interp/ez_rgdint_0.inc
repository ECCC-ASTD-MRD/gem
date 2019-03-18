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
      subroutine ez_rgdint_0(zo,px,py,npts,z,ni,j1,j2)
      implicit none
      
      integer npts,ni,j1,j2,degree,wrap,i,j,n
      real zo(npts),px(npts),py(npts)
      real z(ni,j1:j2)
      
      do n=1,npts
         i = min(ni,max(1,nint(px(n))))
         j = min(j2,max(j1,nint(py(n))))
         
         zo(n)=z(i,j)
      enddo
      
      return
      end
      
