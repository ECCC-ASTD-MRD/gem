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
!**function cherche
!
	integer function ez_cherche(val, tableau, nbelem)
      implicit none
	real val
	integer nbelem
	real tableau(nbelem)
!
!AUTHOR     Yves Chartier               Jan 1991
!
!REVISION
!     Documentation update              Jul 1994      
!
!LANGUAGE   Fortran 77
!
!OBJECT (cherche)
!     Return the index of the table element equal or closest to "val",
!     using binary search. It is assumed that the table is sorted in ascending
!     order.
!
!FILES
!
!ARGUMENTS
! NAMES          I/O  TYPE  A/S        DESCRIPTION
! val             O    R     S         value to search for
! tableau         I    R     A         table of lookup values
! nbelem          I    I     S         length of table
!
!IMPLICIT
!
!MODULES
!
!*
	integer debut, milieu, fin

	debut = 1
	fin = nbelem
	milieu = (debut+fin)*0.5

10	if (milieu.ne.debut) then
          if (val.le.tableau(milieu)) then
	      fin   = milieu
          else
	      debut = milieu
	   endif
          milieu = (debut+fin)*0.5
          goto 10
	endif

	ez_cherche = milieu

	return
	end

