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
      subroutine ez_fillspole(zout, zin, ni, j1, j2, valpole)
      implicit none

      integer i,j,ni,j1,j2
      real zout(ni, 4),zin(ni,j1:j2),valpole
      
      do j=j1,j1+2
         do i=1,ni
            zout(i,j-j1+2) = zin(i,j)
         enddo
      enddo
      
      do i=1,ni
         zout(i,1) = valpole
      enddo
         

      return
      end
      
