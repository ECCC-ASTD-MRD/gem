!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
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
!InTf!
!!integer function RPN_COMM_split(mex,nx,nxg,minx,maxx,nxl,nxlmax,halox,nx0,fill) !InTf!
	integer function RPN_COMM_split(mex,nx,nxg,minx,maxx,nxl,nxlmax,&
     &                   halox,nx0,fill)
	use rpn_comm
	implicit none                                             !InTf!
	integer, intent(IN) :: nx, mex, nxg, halox                !InTf!
	integer, intent(OUT) :: minx,maxx,nxl,nxlmax,nx0          !InTf!
	logical, intent(IN) :: fill                               !InTf!
!
!arguments
!  I	nx	number of PEs (1D)
!  I	mex	PE number of this PE (1D)
!  I	nxg	Global dimension (1D) of data
!  O	minx	start index for this tile (1D) (halo accounted for)
!  O	maxx	end index for this tile (1D) (halo accounted for)
!  O	nxl	length of this tile (1D)
!  O	nxlmax	length of longest tile (tile no 1)
!  I	halox	halo length (1D)
!  O	nx0	start of this tile in glogal space (1D)
!  I	fill	if .true. make sure dimansion is odd (deprecated)
!
!*
!
	integer count(nx),depl(nx), nxn, ierr
	integer, external :: RPN_COMM_limit

	RPN_COMM_split = -1

	ierr = RPN_COMM_limit(mex, nx, 1, nxg, nx0,&
     &     nxn,count, depl)
	if(ierr.ne.0) then
	   write(rpn_u,*) 'RPN_COMM_split: invalid distribution, ABORT'
	   return
	endif

	minx = 1 - halox
	nxl = nxn-nx0+1
	nxlmax = count(1)
	maxx = nxlmax + halox
	if(fill) maxx = maxx + 1 - mod(nxlmax,2)

	RPN_COMM_split = 0

	return
	end function RPN_COMM_split                              !InTf!
