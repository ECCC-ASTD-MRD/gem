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


!**S/R RPN_COMM_swapns: north/south data exchange
!
!	IN:
!	nwen:	integer, number of words to send to north pe
!	wden:	array of integer, words to send to north pe
!	nwes:	integer, number of words to send to south pe
!	wdes:	array of integer, words to send to south pe
!	nwin:	integer, number of words to receive from N pe	
!	nwis:	integer, number of words to receive from S pe	
!	periody:logical, periodicity in N-S direction
!
!	OUT:
!	wdin:	array of integer, words received from N pe
!	wdis:	array of integer, words received from S pe
!	nwrn:	integer, # of words received from N pe
!	nwrs:	integer, # of words received from S pe
!	status:	integer, error flag
!****

	subroutine RPN_COMM_swapns(nwen,wden,nwes,wdes,&
     &          nwin,nwrn,wdin,nwis,nwrs,wdis,periody,&
     &          status)
	use rpn_comm
	implicit none

!	include 'rpn_comm.h'
!	include 'mpif.h'

	integer nwen,wden(nwen),nwes,wdes(nwes)
	integer nwin,nwrn,wdin(nwin),nwis,nwrs,wdis(nwis),status
	logical periody
	integer i,ierr,request,northpe,southpe
	logical north,south,flag
	integer tag_send_to_north, tag_recv_from_south
	integer tag_send_to_south, tag_recv_from_north
	integer req(4),stat(MPI_STATUS_SIZE,4)

	if(pe_ny.eq.1) then
	   do i=1,min(nwes,nwin)
	      wdin(i)=wdes(i)
	   enddo
	   do i=1,min(nwen,nwis)
	      wdis(i)=wden(i)
	   enddo
	   nwrs=min(nwen,nwis)
	   nwrn=min(nwes,nwin)
	   return
	endif
        north=(bnd_north) .and. (.not.periody)
        south=(bnd_south) .and. (.not.periody)
        northpe=pe_id(pe_mex,pe_mey+1)
        southpe=pe_id(pe_mex,pe_mey-1)

!       generate tag for north and south exchange

	tag_send_to_north = northpe
	tag_recv_from_south = pe_me
	tag_send_to_south = pe_tot + southpe
	tag_recv_from_north = pe_tot + pe_me

	if(.not.north) then
	   call mpi_irecv( wdin, nwin, mpi_integer, northpe,&
     &         tag_recv_from_north,pe_defcomm, req(2),ierr )
	endif

	if(.not.south) then
	   call mpi_irecv( wdis, nwis, mpi_integer, southpe,&
     &         tag_recv_from_south,PE_DEFCOMM, req(4),ierr )
	endif

	if(.not.north) then
	   call mpi_isend( wden, nwen, mpi_integer, northpe,&
     &         tag_send_to_north,PE_DEFCOMM, req(1),ierr )
	endif

	if(.not.south) then
	   call mpi_isend( wdes, nwes, mpi_integer, southpe,&
     &         tag_send_to_south,PE_DEFCOMM, req(3),ierr )
	endif

	if((.not.north).and.(.not.south)) then
	   call mpi_waitall(4,req,stat,ierr)
	   call mpi_get_count(stat(1,2),mpi_integer,nwrn,ierr)
	   call mpi_get_count(stat(1,4),mpi_integer,nwrs,ierr)
	   if(nwrn.gt.nwin) then
	      write(rpn_u,*) 'ABORT PE=',pe_me,': Received ',nwrn,&
     &        ' elements from north in a ',nwin,' elements array!'
	      call MPI_abort(PE_DEFCOMM,1,ierr)
	   endif
	   if(nwrs.gt.nwis) then
	      write(rpn_u,*) 'ABORT PE=',pe_me,': Received ',nwrs,&
     &        ' elements from south in a ',nwis,' elements array!'
	      call MPI_abort(PE_DEFCOMM,1,ierr)
	   endif

	else
	 
!       cannot be south and north at the same time (case pe_ny=1)
	   if(north) then
	     call mpi_waitall(2,req(3),stat,ierr)
	     call mpi_get_count(stat(1,2),mpi_integer,nwrs,ierr)
	   if(nwrs.gt.nwis) then
	      write(rpn_u,*) 'ABORT PE=',pe_me,': Received ',nwrs,&
     &        ' elements from south in a ',nwis,' elements array!'
	      call MPI_abort(PE_DEFCOMM,1,ierr)
	   endif
	   else ! south...
	     call mpi_waitall(2,req(1),stat,ierr)
	     call mpi_get_count(stat(1,2),mpi_integer,nwrn,ierr)
	   if(nwrn.gt.nwin) then
	      write(rpn_u,*) 'ABORT PE=',pe_me,': Received ',nwrn,&
     &        ' elements from north in a ',nwin,' elements array!'
	      call MPI_abort(PE_DEFCOMM,1,ierr)
	   endif
	   endif
	endif
 510	continue
        return
        end

