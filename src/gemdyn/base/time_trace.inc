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
type, BIND(C) :: time_context  ! wrappper to give a type to a C pointer
  type(C_PTR) :: t
end type

interface

  subroutine time_trace_init(t) bind(C,name='TimeTraceInit') ! create and initialize a new time trace context
    import :: time_context
    implicit none
    type(time_context), intent(OUT) :: t             ! opaque time context pointer (passed to other routines)
  end subroutine time_trace_init

  subroutine time_trace_dump_text(t, filename, ordinal)   ! dump timings int file filename_nnnnnn.txt (nnnnnn from ordinal)
    import :: time_context
    implicit none
    type(time_context), intent(IN) :: t              ! opaque time context pointer (from time_trace_init)
    character(len=*), intent(IN) :: filename         ! file name prefix (will be trimmed to remove trailing blanks if any)
    integer, intent(IN) :: ordinal                   ! numbered extension to file name (nnnnnn) (normally MPI rank)
  end subroutine time_trace_dump_text

  subroutine time_trace_dump_binary(t, filename, ordinal)   ! dump timings int file filename_nnnnnn.txt (nnnnnn from ordinal)
    import :: time_context
    implicit none
    type(time_context), intent(IN) :: t              ! opaque time context pointer
    character(len=*), intent(IN) :: filename         ! file name prefix (will be trimmed to remove trailing blanks if any)
    integer, intent(IN) :: ordinal                   ! numbered extension to file name (nnnnnn) (normally MPI rank)
  end subroutine time_trace_dump_binary

  subroutine time_trace_step(t, n) bind(C,name='TimeTraceStep')  ! set step value for subsequent calls to time_trace
    import :: time_context, C_INT
    implicit none
    type(time_context), intent(IN), value :: t              ! opaque time context pointer (from time_trace_init)
    integer(C_INT), intent(IN), value :: n                  ! time step number
  end subroutine time_trace_step

  subroutine time_trace_barr(t, tag, barrier, comm, barrier_code)  ! insert a new time trace entry (2 entries if barrier is true)
    import :: time_context
    type(time_context), intent(IN) :: t              ! opaque time context pointer (from time_trace_init)
    integer, intent(IN) :: tag                       ! tag number for this timing point (MUST be >0 and <128M)
    logical, intent(IN) :: barrier                   ! if true, call MPI_barrier with timing points before and after
    integer, intent(IN) :: comm                      ! MPI communicator (only used if barrier flag is true)
    external :: barrier_code
  end subroutine time_trace_barr

  subroutine time_trace(t, tag)  bind(C,name='TimeTrace') ! insert a new time trace entry (no barrier)
    import :: time_context, C_INT
    implicit none
    type(time_context), intent(IN), value :: t              ! opaque time context pointer (from time_trace_init)
    integer(C_INT), intent(IN), value :: tag                ! tag number for this timing point (MUST be >0 and <128M)
  end subroutine time_trace

  subroutine time_trace_get_buffers(t, array, larray, n) bind(C,name='TimeTraceGetBuffers')
    import :: C_PTR , C_INT, time_context
    implicit none
    type(time_context), intent(IN), value :: t              ! opaque time context pointer (passed to other routines)
    type(C_PTR), dimension(n), intent(OUT) :: array         ! to receive pointers to buffers
    integer(C_INT), dimension(n), intent(OUT) :: larray     ! to receive lengths of buffers
    integer(C_INT), intent(IN), value :: n                  ! size of array and larray 
  end subroutine time_trace_get_buffers

end interface
