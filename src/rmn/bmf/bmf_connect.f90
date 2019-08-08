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
  function bmf_connect(no_pe) result(ierr)
  use bmf_modsplit
  implicit none
  integer, intent(IN) ::  no_pe
  integer :: ierr

  integer fnom
  external fnom

  if(split_unit(no_pe).eq.0) then
     ierr=FNOM(split_unit(no_pe),split_files(no_pe),'SEQ/UNF',0)
  else
     ierr=0
  endif

  end function bmf_connect
