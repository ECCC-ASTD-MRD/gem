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

      subroutine ez_calcpoleval(poleval, z, ni, ax, grtyp, grref)
      implicit none
      
      integer ni, wrap
      real poleval,z(ni),ax(ni)
      character*1 grtyp, grref
      
      integer i
      
      if (grtyp.eq.'Z'.and.grref.eq.'E') then
         poleval = 0.0
         do i=1,ni-1
            poleval = poleval + z(i)*(ax(i+1)-ax(i))
         enddo
         if ((ax(ni)-ax(1)) /= 0.0) then
            poleval = poleval / (ax(ni)-ax(1))
         endif
         return
      endif

      poleval = 0.0
      do i=1,ni
         poleval = poleval + z(i)
      enddo
      poleval = poleval / (1.0 * ni)

      return
      end

