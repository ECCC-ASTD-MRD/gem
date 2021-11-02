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
!*********************************************************************
!**s/r ez_glat - calcul des latitudes d'une grille gaussienne
!
!  auteur: Yves Chartier. Mars 1991.
!****
!
      subroutine ez_glat(latroots,groots,nj,hem)
      implicit none

      integer nj, hem,ier

#include "qqqpar.cdk"
#include "pi.cdk"
      
      external dgauss
      real latroots(*), groots(*)
    
      integer j,npoly
      real temp

      if (hem .ne. GLOBAL) then
         npoly = nj * 2
      else
         npoly = nj
      endif
      call dgauss(npoly, groots, global)
      
      do j=1,npoly/2
         temp = groots(j)
         groots(j)=groots(npoly+1-j)
         groots(npoly+1-j)=temp
      enddo
      
      if (hem .ne. nord) then
         do j=1,nj
            latroots(j)= 90. - rdtodg * acos(groots(j))
         enddo
         
      endif
      
      if (hem .eq. NORD) then
         do 20 j=1,nj
            latroots(j)=90.0-rdtodg*acos(groots(j+nj))
 20      continue
      endif

      return
      end
