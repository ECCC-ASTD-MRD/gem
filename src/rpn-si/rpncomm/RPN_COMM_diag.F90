!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/

	subroutine RPN_COMM_diag(Sub,Array,Minx,Maxx,Miny,Maxy,Ni,Nj,Nk,&
     &             R_8,Nr,Tmp_8,Fast)
	use rpn_comm
	implicit none
	external sub
	integer Minx,Maxx,Miny,Maxy,Ni,Nj,Nk,Nr
	logical Fast
	real array(Minx:Maxx,Miny:Maxy,Nk)
	real *8 R_8(Nk,Nr), Tmp_8(Miny:Maxy,Nk,Nr)
!
!	include 'rpn_comm.h'
!	include 'mpif.h'
!
	integer nelem,ierr
	integer status(mpi_status_size)
!
!	pass 1 WEST to EAST
!
	nelem = (Maxy-Miny+1)*Nk*Nr
!
	if(bnd_west) then  ! initialize WEST to EAST collect
	   call Sub(Array,Minx,Maxx,Miny,Maxy,Ni,Nj,Nk,R_8,Nr,Tmp_8,&
     &              1,.true.)
	else  ! get result from WEST neighbor continue collection
	   call MPI_RECV(Tmp_8,nelem,MPI_DOUBLE_PRECISION,&
     &                   pe_id(pe_mex-1,pe_mey),pe_id(pe_mex-1,pe_mey),&
     &                   pe_defcomm,status,ierr)
	   call Sub(Array,Minx,Maxx,Miny,Maxy,Ni,Nj,Nk,R_8,Nr,Tmp_8,&
     &              1,.false.)
	endif
	if(.not. bnd_east) then  ! send result to EAST neighbor
	   call MPI_SSEND(Tmp_8,nelem,MPI_DOUBLE_PRECISION,&
     &                    pe_id(pe_mex+1,pe_mey),pe_me,&
     &                    pe_defcomm,ierr)
	endif
!
!	broadcast partial result to row
!
	call MPI_BCAST(Tmp_8,nelem,MPI_DOUBLE_PRECISION,&
     &                 pe_id(pe_nx-1,pe_mey),pe_myrow,ierr)
!
!	pass 2 NORTH to SOUTH  (EAST column only)
!
	if(bnd_north) then  ! initialize NORTH to SOUTH collect
	  call Sub(Array,Minx,Maxx,Miny,Maxy,Ni,Nj,Nk,R_8,Nr,Tmp_8,&
     &              2,.true.)
	else  ! get result from NORTH neighbor continue collection
	  call MPI_RECV(R_8,Nr*Nk,MPI_DOUBLE_PRECISION,&
     &                   pe_id(pe_mex,pe_mey+1),pe_id(pe_mex,pe_mey+1),&
     &                   pe_defcomm,status,ierr)
	  call Sub(Array,Minx,Maxx,Miny,Maxy,Ni,Nj,Nk,R_8,Nr,Tmp_8,&
     &              2,.false.)
	endif
	if(.not. bnd_south) then  ! send result to SOUTH neighbor
	  call MPI_SSEND(R_8,Nr*Nk,MPI_DOUBLE_PRECISION,&
     &                    pe_id(pe_mex,pe_mey-1),pe_me,&
     &                    pe_defcomm,ierr)
	endif
!
!	broadcast final result to column
!
	call MPI_BCAST(R_8,Nr*Nk,MPI_DOUBLE_PRECISION,&
     &                 pe_id(pe_mex,0),pe_mycol,ierr)
	return
	end
