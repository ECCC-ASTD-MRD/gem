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
program test_trace   ! Fortran test program for library (C and Fortran code)
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
  integer :: i, tag, rank, nbeads, nbent, nused, nprc, prc, j
  type(time_context) :: t
  type(C_PTR), dimension(10) :: array
  integer(C_INT), dimension(10) :: larray
  integer(C_INT), dimension(:,:), allocatable :: blob
  integer(C_INT), dimension(:), pointer   :: local
  integer(C_INT), dimension(:,:), allocatable, target :: global
  external :: MPI_barrier

  rank = 0
  nprc = 1
#if ! defined(NO_MPI)
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprc, ierr)
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
  array(1) = C_NULL_PTR
  array(1) = time_trace_get_buffer_data(t, nbeads, nbent, 1)
  if(C_ASSOCIATED(array(1))) then
    print *,'nbeads =',nbeads,' , nbent =',nbent
    call time_trace_single_text(array(1), nbent, 'time_dump'//achar(0), rank + 1)
  else
    print *,'ERROR getting timing data'
    goto 777
  endif
  call c_f_pointer(array(1), local, [nbent])  ! C pointer to Fortran pointer
  if(rank == 0) then
    allocate(global(nbent,max(nprc,2)))
    print *,'global allocated ',nbent,max(nprc,2)
  endif
#if ! defined(NO_MPI)
  call MPI_gather(local, nbent, MPI_INTEGER, global, nbent, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
#else
  global(:,1) = local
  global(:,2) = local
#endif
  if(rank == 0) then
!     print *,'nbeads =',nbeads,' , nbent =',nbent
!     call time_trace_single_text(array(1), nbent, 'time_dump'//achar(0), 1)
    allocate(blob(2+2*nprc+2,100))
    print *,'blob allocated ',2+2*nprc+2,100
    blob = 0
!     nused = time_trace_expand(array(1), nbent, blob, blob(3,1), 6, 14, 0)
    nused = time_trace_expand(C_LOC(global(1,1)), nbent, blob, blob(3,1), 2+2*nprc+2, 14, 0)
    print *,'nused =',nused

    do prc = 2, nprc
!       nused = time_trace_expand(array(1), nbent, blob, blob(5,1), 6, 14, 1)
      nused = time_trace_expand(C_LOC(global(1,prc)), nbent, blob, blob(1+2*prc,1), 2+2*nprc+2, 14, 1)
      print *,'nused =',nused
    enddo

    do i = 1, abs(nused)
      if(blob(2,i) == -1) then   ! step flag, combine 2 unsigned 32 bit integers into a 64 bit integer
        print 100,blob(1,i), blob(2,i), &
                  ((i8_from_2_i4(blob(2*j+1,i), blob(2*j+2,i)), 0) , j=1,nprc)
      else
        print 101,blob(1:2+nprc*2,i)
      endif
    enddo
  endif
777 continue
#if ! defined(NO_MPI)
  call MPI_finalize(ierr)
#endif

  stop
100 format(2I10,5(I18,I2))
101 format(12I10)
end program
#if defined(NO_MPI)
  subroutine MPI_barrier(dummy1, dummy2)
    integer, intent(IN) :: dummy1, dummy2
    if(loc(dummy1) == loc(dummy2)) print *,'OOPS !!'
  end subroutine MPI_barrier
#endif
#endif
