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
	subroutine RPN_COMM_qadl(ipe,ni,nj,i0,in,j0,jn)   !InTf!
	use rpn_comm
	implicit none                                     !InTf!
    integer, intent(OUT) :: i0,in,j0,jn               !InTf!
    integer, intent(IN) :: ipe,ni,nj                  !InTf!
    integer :: ilim, jlim, ifact
!	include 'rpn_comm.h'
!	include 'mpif.h'
	ilim = (ni + pe_nx -1) / pe_nx
	jlim = (nj + pe_ny -1) / pe_ny
	ifact = ipe/pe_nx
	j0 = jlim*(ipe/pe_nx) + 1
	jn = j0 + jlim -1
	i0 = ilim*mod(ipe,pe_nx) + 1
	in = i0 + ilim -1
	in=min(in,ni)
	jn=min(jn,nj)
	return
	end subroutine RPN_COMM_qadl   !InTf!
