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
      subroutine ez_fillaxzx(axdst, ax, ni, i1, i2)
      integer ni, i1, i2
      real axdst(i1:i2), ax(ni)

      integer i

      print *, ni, i1, i2
      do i=1,ni
         axdst(i) = ax(i)
         print *, axdst(i), ax(i)
      enddo
      axdst(0) = ax(ni-1) - 360.0
      axdst(ni)=ax(1)+360.0
      axdst(ni+1)=ax(2)+360.0

      return
      end
