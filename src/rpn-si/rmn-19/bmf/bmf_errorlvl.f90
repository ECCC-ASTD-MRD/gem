!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2005  Environnement Canada
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
  subroutine bmf_error(errlvl)
!
! SUBROUTINE BMF_error, L. Corbeil
!
! ARGUMENTS
!
! DESCRIPTION
!
! Subroutine that initializes/modifies the error level 
!
  use bmf_mod
  implicit none
  integer errlvl

! 0= Tous
! 1= moins warnings
! 2= moins erreurs
! >5 aucun (de grace, ne faites pas ca!!!)


  bmf_errorlvl=errlvl

  end subroutine bmf_error
