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
!****************************************************************************
      subroutine ez_polrint(vpolnor,vpolsud,zi,ni,nj,grtyp,grref,hem,vecteur,ax,ay)
      
      implicit none

      integer ni,nj,hem
      real zi(ni,nj)
      real ax(ni),ay(nj)
      character*1 grtyp,grref
      logical vecteur      
      integer n
      real vpolnor, vpolsud, sum

      
#include "qqqpar.cdk"

!
!Definition des variables locales
!

      integer i,j
      integer n1, n2, s1, s2

      if(vecteur) then
         return
      endif
      
      if(grtyp.eq.'L'.or.grtyp.eq.'N'.or.grtyp.eq.'S'.or.grtyp.eq.'!'.or.(grtyp.eq.'Z'.and.grref.ne.'E')) then
         return
      endif
      
      if (grtyp .eq. 'B') then
         vpolnor = zi(1, nj)
         vpolsud = zi(1, 1)
         return
      endif
      
      if(grtyp .eq. 'A' .or. grtyp .eq. 'G'.or. (grtyp.eq.'Z'.and.grref.eq.'E')) then
         sum = 0.0
         do i=1,ni
            sum = sum + zi(i,nj)
         enddo
         vpolnor = sum / (1.0*ni)
         
         sum = 0.0
         do i=1,ni
            sum = sum + zi(i,1)
         enddo
         vpolsud = sum / (1.0*ni)
      endif
      
      return
      end
