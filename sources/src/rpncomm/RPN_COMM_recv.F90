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

        SUBROUTINE RPN_COMM_recv(buffer,count,datatype,source,tag,&
     &                           com,status,ierr)
!	Luc Corbeil, 2000-11-21
!	mpi recv

        implicit none
        integer count,comm,ierr,tag,source
        integer status
        integer buffer
        integer datyp
        character(len=*) datatype,com
        integer RPN_COMM_datyp,RPN_COMM_comm
        logical RPN_COMM_grank

        include 'mpif.h'

        integer status2(MPI_STATUS_SIZE)
        datyp=rpn_comm_datyp(datatype)
	comm=rpn_comm_comm(com)

        call mpi_recv(buffer,count,datyp,source,tag,comm,status2,ierr)
        call mpi_get_count( status2, datyp, status, ierr)

	return
	end
