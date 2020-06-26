! RPN_COMM - Library of useful routines for C and FORTRAN programming
! Copyright (C) 2020  Division de Recherche en Prevision Numerique
!                     Environnement Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.
!
subroutine RPN_COMM_win_allocate_shared(comm, msize, win, baseptr, err)  !InTf! 
  use ISO_C_BINDING
  implicit none
!!  import :: C_PTR                            !InTf!   
  integer, intent(IN)  :: comm               !InTf!   communicator
  integer, intent(IN)  :: msize              !InTf!   size in integers of shared memory area
  integer, intent(OUT) :: win                !InTf!   window communicator
  type(C_PTR), intent(OUT) :: baseptr        !InTf!   base address of shared memory area
  integer, intent(OUT)  :: err               !InTf!   status RPN_COMM_OK or RPN_COMM_ERROR
  
  include 'mpif.h'
  include 'RPN_COMM_constants.inc'

  integer(KIND=MPI_ADDRESS_KIND) :: wsiz
  integer :: dispunit, ierr, rank

  err = RPN_COMM_ERROR      ! precondition for failure

  call MPI_comm_rank(comm, rank, ierr)
  if(ierr .ne. MPI_SUCCESS) return
  if(rank == 0) then          ! everythng allocated by rank 0
    wsiz = msize
  else
    wsiz = 0
  endif
  dispunit = 4                ! words (integers/floats)
  wsiz = wsiz * dispunit    ! size in Bytes
  call MPI_win_allocate_shared(wsiz, dispunit, MPI_INFO_NULL, comm, baseptr, win, ierr)   ! create window and shared memory area
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_win_shared_query(win, MPI_PROC_NULL, wsiz, dispunit, baseptr, ierr)            ! get local address for shared memory area
  if(ierr .ne. MPI_SUCCESS) return
  
  err = RPN_COMM_OK
  return
end subroutine RPN_COMM_win_allocate_shared  !InTf! 
