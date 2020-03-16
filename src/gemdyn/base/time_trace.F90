! Hopefully useful software for FORTRAN
! Copyright (C) 2019  Division de Recherche en Prevision Numerique
!                     Environnement Canada
!
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
module time_trace_mod
  use ISO_C_BINDING
  implicit none
#define NEED_PRIVATE
#include <time_trace.hf>
#undef NEED_PRIVATE

end module
!
! user callable subroutines that need to be written in Fortran because of interfacing issues
!
subroutine time_trace_barr(t, tag, barrier, comm, barrier_code)  ! insert a new time trace entry (2 entries if barrier is true)
  use ISO_C_BINDING
  use time_trace_mod
  implicit none
  external :: barrier_code
  type(time_context), intent(IN) :: t              ! opaque time context pointer (from time_trace_init)
  integer, intent(IN) :: tag                       ! tag number for this timing point (MUST be >0 and <128M)
  integer, intent(IN) :: comm                      ! MPI communicator (only used if barrier flag is true)
  logical, intent(IN) :: barrier                   ! if true, call MPI_barrier with timing points before and after

  integer :: ierr
  integer(kind=8) :: tt1, tt2

  tt1 = what_time_is_it() 
  if(barrier) then
    call barrier_code(comm, ierr)                  ! barrier call if needed
    tt2 = what_time_is_it()                        ! time entry after barrier (used to measure imbalance)
  else
    tt2 = -1                                       ! no barrier, set to zero
  endif
  call new_time_tag(t, tag, tt1, tt2)              ! insert into timing table

  return
end subroutine time_trace_barr

subroutine time_trace_dump_binary(t, filename, ordinal)   ! dump timings int file filename_nnnnnn.txt (nnnnnn from ordinal)
  use ISO_C_BINDING
  use time_trace_mod
  implicit none
  type(time_context), intent(IN) :: t              ! opaque time context pointer
  character(len=*), intent(IN) :: filename         ! file name prefix (will be trimmed to remove trailing blanks if any)
  integer, intent(IN) :: ordinal                   ! numbered extension to file name (nnnnnn) (normally MPI rank)

  interface
    subroutine c_dump_binary(t, filename, ordinal) bind(C,name='TimeTraceDumpBinary')
      import :: time_context, C_CHAR, C_INT
      type(time_context), intent(IN), value :: t
      character(C_CHAR), dimension(*), intent(IN) :: filename
      integer(C_INT), intent(IN), value :: ordinal
    end subroutine c_dump_binary
  end interface

  call c_dump_binary(t, trim(filename)//achar(0), ordinal) ! call C routine with properly terminated string

  return
end subroutine time_trace_dump_binary

subroutine time_trace_dump_text(t, filename, ordinal)   ! dump timings int file filename_nnnnnn.txt (nnnnnn from ordinal)
  use ISO_C_BINDING
  use time_trace_mod
  implicit none
  type(time_context), intent(IN) :: t              ! opaque time context pointer
  character(len=*), intent(IN) :: filename         ! file name prefix (will be trimmed to remove trailing blanks if any)
  integer, intent(IN) :: ordinal                   ! numbered extension to file name (nnnnnn) (normally MPI rank)

  interface
    subroutine c_dump_text(t, filename, ordinal) bind(C,name='TimeTraceDumpText')
      import :: time_context, C_CHAR, C_INT
      type(time_context), intent(IN), value :: t
      character(C_CHAR), dimension(*), intent(IN) :: filename
      integer(C_INT), intent(IN), value :: ordinal
    end subroutine c_dump_text
  end interface

  call c_dump_text(t, trim(filename)//achar(0), ordinal) ! call C routine with properly terminated string
end subroutine time_trace_dump_text

#if defined SELF_TEST
program test_trace
  use ISO_C_BINDING
  implicit none
#include <time_trace.hf>
#if ! defined(NO_MPI)
  include 'mpif.h'
  integer :: ierr
#else
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: MPI_COMM_NULL = -1
#endif
  integer :: i, tag, rank
  type(time_context) :: t
  type(C_PTR), dimension(10) :: array
  integer(C_INT), dimension(10) :: larray
  external :: MPI_barrier

  rank = 0
#if ! defined(NO_MPI)
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, rank, ierr)
#endif
  call time_trace_init(t)
  print *,'-----------------------------'
  call time_trace(t, 0)
  call time_trace_step(t, 0)
  print *,'============================='
  call time_trace_barr(t, 1, .true., MPI_COMM_WORLD, MPI_barrier)
  do i = 1, 3
    call time_trace_step(t, i)
    print *,'+++++++ start step =',i
    tag = 10*i+1
    call time_trace_barr(t, tag, .true., MPI_COMM_WORLD, MPI_barrier)
    print *,'+++++++ end   step =',i
    tag = 10*i+2
    call time_trace(t, tag)
  enddo
  call time_trace_step(t, 4)
  print *,'============================='
  call time_trace(t, 999999)
  call time_trace_dump_text(t, 'time_list', rank)
  call time_trace_dump_binary(t, 'time_list', rank)
  call time_trace_get_buffers(t, array, larray, 10)
  write(6,'(10I6)')larray
#if ! defined(NO_MPI)
  call MPI_finalize(ierr)
#endif

  stop
end program
#if defined(NO_MPI)
  subroutine MPI_barrier(dummy1, dummy2)
    integer, intent(IN) :: dummy1, dummy2
    if(loc(dummy1) == loc(dummy2)) print *,'OOPS !!'
  end subroutine MPI_barrier
#endif
#endif
