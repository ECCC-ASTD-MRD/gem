!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2007  Environnement Canada
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
!object
!  A few fortran helpers to deal with logical fortran type in C
!author
!  Stephane Chamberland, aug 2007
!** end of rpn-doc sections


!**s/r fc2_logical2int
subroutine ftn2c_logical2int(dest,src,n)
  implicit none
!arguments
  integer :: n
  integer, dimension(n) :: dest
  logical, dimension(n) :: src
!object
!  handler function call in C to store/retreive fortran logical
!** end of rpn-doc sections
  dest = 0
  where(src) dest = 1
  return
end subroutine ftn2c_logical2int


!**s/r f2c_int2logical
subroutine ftn2c_int2logical(dest,src,n)
  implicit none
!arguments
  integer :: n
  integer, dimension(n) :: src
  logical, dimension(n) :: dest
!object
!  handler function call in C to store/retreive fortran logical
!** end of rpn-doc sections
  dest = src.ne.0
  return
end subroutine ftn2c_int2logical


!**s/r f_logical_move
subroutine ftn2c_logical_move(dest,src,n)
  implicit none
!arguments
  integer :: n
  logical, dimension(n) :: dest,src
!object
!  handler function call in C to store/retreive fortran logical
!** end of rpn-doc sections
  dest=src
  return
end subroutine ftn2c_logical_move
