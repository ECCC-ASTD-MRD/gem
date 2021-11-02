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
      subroutine ez_calcpolemodule(polemod, uu, vv, ni, ax, grtyp, grref)
      implicit none

      integer ni
      real polemod,uu(ni),vv(ni), ax(ni), module
      character*1 grtyp, grref

      integer i

      if (grtyp.eq.'Z'.and.grref.eq.'E') then
         polemod = 0.0
         do i=1,ni-1
            module = sqrt(uu(i)*uu(i)+vv(i)*vv(i))
            polemod = polemod + module*(ax(i+1)-ax(i))
         enddo
         polemod = polemod / 360.0
         return
      endif

      polemod = 0.0
      do i=1,ni
         polemod = polemod + sqrt(uu(i)*uu(i)+vv(i)*vv(i))
      enddo
      polemod = polemod / (1.0 * ni)

      return
      end

