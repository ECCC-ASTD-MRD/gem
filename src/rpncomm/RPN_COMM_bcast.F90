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

        SUBROUTINE RPN_COMM_bcast(buffer,count,datatype,root,com,ierr)
!	Luc Corbeil, 2000-11-20
!	mpi broadcast
        implicit none
        integer count,root,comm,ierr,group
        integer buffer
        integer datyp
        character(len=*) datatype,com
        integer RPN_COMM_datyp,RPN_COMM_comm
        logical RPN_COMM_grank
!*
!        include 'rpn_comm.h'
!        include 'mpif.h'
        comm=rpn_comm_comm(com)
        datyp=rpn_comm_datyp(datatype)

        if(.not.RPN_COMM_grank(com)) return
        call mpi_bcast(buffer,count,datyp,root,comm,ierr)

	return
	end

!InTf!
      SUBROUTINE RPN_COMM_bcastc(buffer,count,datatype,root,com,ierr) !InTf!
!	Luc Corbeil, 2000-11-20
!	mpi broadcast
      implicit none                                                   !InTf!
      integer, intent(IN) :: count,root                               !InTf!
      integer, intent(OUT) :: ierr                                    !InTf!
      character(len=*), intent(INOUT) :: buffer                       !InTf!
      character(len=*), intent(IN) :: datatype,com                    !InTf!
      integer comm,group,datyp
      integer, external :: RPN_COMM_datyp,RPN_COMM_comm
      logical, external :: RPN_COMM_grank
!*
!        include 'rpn_comm.h'
!        include 'mpif.h'
      comm=rpn_comm_comm(com)
      datyp=rpn_comm_datyp(datatype)

      if(.not.RPN_COMM_grank(com)) return
      call mpi_bcast(buffer,count,datyp,root,comm,ierr)

      return
      end SUBROUTINE RPN_COMM_bcastc                                  !InTf!
