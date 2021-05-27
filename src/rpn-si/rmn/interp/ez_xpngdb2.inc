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
      subroutine ez_xpngdb2(zout,zi,ni,nj,j1,j2,hem,symetrie)
      
      implicit none
#include "qqqpar.cdk"

      external permut
      
      integer i1,i2,j1,j2,ni,nj
      integer hem,symetrie
      character*1 grtyp
      integer i,j,njsur2
      
      
      real zout(ni,j1:j2)
      real zi(ni,nj)
      integer sym
      real sign
      integer ier,ii
      
      if (symetrie.eq.0) then
         sign = -1.0
      else
         sign = 1.0
      endif
      
      if (hem .eq. nord) then
         do j=1,nj
            do i=1,ni
               zout(i,j)  = zi(i,j)
            enddo
         enddo

         do j=2,nj
            do i=1,ni
               zout(i,2-j)  = sign * zi(i,j)
            enddo
         enddo
      endif

      if (hem .eq. sud) then
         do j=1,nj
            do i=1,ni
               zout(i,j)  = zi(i,j)
            enddo
         enddo
         
         do j=2,nj
            do i=1,ni
               zout(i,nj+j-1)  = sign * zi(i,nj-j+1)
            enddo
         enddo
      endif
      
         
      return
      end
      
