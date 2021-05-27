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
!****iP* RPN_COMM/io_pe  management routines for  IO PE sets (internal)
! FUNCTION
!   Create, free, manage sets of processes uses to perform IO (or other) operations
!   in the RPN_COMM library (version 4.5.16 and above)
! AUTHOR
!   M.Valin Recherche en Prevision Numerique 2015
! DESCRIPTION
!
!   internal module RPN_COMM_io_pe_tables (output arguments in UPPER CASE)
!   (not directly callable by user, called by user callable routines)
!   (contains the static ioset table and control variables)
!   - .t /.f.   = is_valid_io_setno(setno)
!   - ordinal   = is_io_pe(setno)
!   - comm      = io_pe_comm(setno)
!   - list(:,2) = io_pe_coord(setno)
!   - ivalue    = io_pe_call(setno,callback,argv)
!   - npes      = io_pe_size(setno)
!   - ngrps     = io_pe_groups(setno)
!   - setno     = free_ioset(setno)
!   - setno     = create_ioset(npes,method)
!******
!****P* rpn_comm/io_pe  (IO PE sets management package)
! FUNCTION
!   Create, free, manage sets of processes used to perform IO (or other) operations
!   in the RPN_COMM library (version 4.5.16 and above)
! AUTHOR
!   M.Valin Recherche en Prevision Numerique 2015
! DESCRIPTION
!
!  user callable routines / functions (output/inout arguments in UPPER CASE)
!                 CREATE / FREE
!   - setno      = RPN_COMM_create_io_set(npes,method)
!   - call         RPN_COMM_make_io_pe_list(X,Y,npes,penx,peny,method)
!   - setno      = RPN_COMM_free_io_set(setno)
!                  INFORMATION
!   - .t /.f.    = RPN_COMM_is_valid_io_setno(setno)
!   - ordinal_s  = RPN_COMM_is_io_pe(setno)
!   - ordinal_g  = RPN_COMM_io_pe_gridid(setno,n)
!   - comm       = RPN_COMM_io_pe_comm(setno)
!   - npes       = RPN_COMM_io_pe_size(setno)
!   - grps       = RPN_COMM_io_pe_groups(setno)
!   - list(:)    = RPN_COMM_io_pe_idlist(setno)
!   - coord(:,:) = RPN_COMM_io_pe_coord(setno)
!                    OPERATE
!   - ivalue     = RPN_COMM_io_pe_callback(setno,callback,argv)
!   - call         RPN_COMM_io_pe_bcast(BUFFER,count,datatype,root,setno,IERR)
!                   VALIDATION
!   - 0/-1       = RPN_COMM_check_ioset(setno, x ,y, npes, penx, peny, peme, diag)
!   - logical    = RPN_COMM_io_pe_valid_set(X,Y,npes,penx,peny,diag,method)
!
!     setno      set number, assigned at IO PE set creation time by RPN_COMM_create_io_set
!     method     PE spread method, integer, must be 0 for now
!     peme       ordinal in "GRID" of this PE
!     diag       if .true. and peme == 0, print diagnostic messages
!     npes       number of PEs in IO set
!     grps       number of "groups" ni thi IO PE set
!     comm       MPI communicator for this IO PE set (integer)
!     ordinal_s  ordinal in IO PE set communicator  (origin 0)
!     ordinal_g  ordinal in "GRID"  (origin 0)
!     list       one dimensional integer array, list(I) is the ordinal in the "GRID"
!                communicator of the Ith IO PE in set. (size of list is npes)  (origin 0)
!     coord      two dimensional integer array (npes,2) containing the X and Y coordinates
!                in the "GRID" of the IO PEs in the set.   (origin 0)
!                X coordinates in coord(:,1), Y coordinates in coord(:,2)
!     argv       "blind" argument of type(C_PTR) that will be passed to callback
!     ivalue     integer return value of callback function
!     callback   the name of a user supplied function, passed to RPN_COMM_io_pe_callback.
!                RPN_COMM_io_pe_callback MUST be called on ALL members of the IO PE set
!                but may also be called by non members of the IO set. the callback function
!                will only be called on the members of the IO set. the return value of 
!                RPN_COMM_io_pe_callback will be 0 on non members and the return value of 
!                callback on the members of the IO set. the call to callback is done as
!                ivalue = callback(argv,setno,comm,pe_me,x,y,npes)
!                pe_me is the ordinal of the PE in the IO PE set communicator (NOT "GRID")
! Notes
!     include 'RPN_COMM.inc'
!     is necessary to get the correct interface description
!
!     the normal sequence of operations is 
!     1 - RPN_COMM_create_io_set
!     2 - get info, do operations on set "setno"
!     3 - RPN_COMM_free_io_set
!
!     for now a maximum of 16 IO PE sets may be defined at any point in time
!     (when an IO set has been freed, it does not count anymore)
!
!     RPN_COMM_make_io_pe_list and RPN_COMM_check_ioset are really internal routines,
!     users should rather call RPN_COMM_io_pe_valid_set that is intended for public use
!     and better documented
!
!     FOR MODEL DEVELOPERS:
!
!     The IO PEs distribution is done in a hopefully optimal way to maximize IO performance
!     when using a large number of nodes.
!     For output purposes, between 1 and 2 IO PEs per node is usually considered optimal.
!     For input purposes one might use more PEs per node if the input operation involves
!     heavvy computations (like horizontal interpolation).
!     If the number of PEs in the X direction is pex and the number of PEs in the Y direction
!     is pey, min(pex,pey)*min(pex,pey) is likely to be the maximum number of PEs in a set.
!     The RPN_COMM_io_pe_valid_set subroutine (totally standalone) can be used to see if
!     the requested IO PE set is possible given pex and pey.
!******
#define IN_RPN_COMM_io_pes
#if defined(SELF_TEST)
#define STAND_ALONE
#endif
!
! STAND_ALONE mode is really used for debugging
! in normal mode, these routines rely on internal module rpn_comm
!
module RPN_COMM_io_pe_tables   ! this module may also be used by data distribution routines
  use ISO_C_BINDING
#if ! defined(STAND_ALONE)
  use rpn_comm
#else
#define RPN_COMM_OK 0
#define RPN_COMM_ERROR -1
#endif
#include "RPN_COMM_interfaces_int.inc"
  type :: RPN_COMM_io_set
    integer, dimension(:), pointer :: x    ! x coordinate in grid of IO PEs in this set
    integer, dimension(:), pointer :: y    ! y coordinate in grid of IO PEs in this set
    integer :: ioset                       ! ID for this IO PE set
    integer :: comm                        ! MPI communicator  for this IO PE set     if it belongs to the set (MPI_COMM_NULL otherwise)
    integer :: me                          ! ordinal of this PE in above communicator if it belongs to the set (-1 otherwise)
    integer :: megrid                      ! ordinal of this PE in grid               if it belongs to the set (-1 otherwise)
    integer :: npe                         ! number of IO PEs in this set
    integer :: ngroups                     ! number of PE groups in this set (size of a group may not exceed MIN(pe_nx,pe_ny) )
    integer :: groupsize                   ! size of groups (last group may be shorter)
    integer :: reserved                    ! unused for now
  end type
  integer, save :: iosets=0
  integer, PARAMETER :: MAXIOSETS=16
  type(RPN_COMM_io_set), dimension(MAXIOSETS), save :: io_set
contains
  function is_valid_io_setno(setno) result(valid)  ! is this a valid IO set ?
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    logical :: valid

    valid = .false.
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    valid = (io_set(setno)%ioset == setno)      ! set still valid ?
  end function is_valid_io_setno
!
  function is_io_pe(setno) result(ordinal)  ! is this pe part of IO PE set setno. if yes return rank in set, else return -1
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: ordinal
    integer :: i
#if defined(STAND_ALONE)
    integer :: pe_mex, pe_mey, pe_me, ierr
    integer, external :: RPN_COMM_mype
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
#endif
!
    ordinal = -1                                ! preset for failure
    if( .not. is_valid_io_setno(setno) ) return   ! not a valid IO set
!
    do i = 1 , io_set(setno)%npe                ! loop over IO PE set population
      if(pe_mex == io_set(setno)%x(i) .and. pe_mey == io_set(setno)%y(i)) then  ! I am an IO pe
        ordinal =  io_set(setno)%me             ! return my ordinal in IO PE communicator
        return
      endif
    enddo
    
  end function is_io_pe
!
  function io_pe_comm(setno) result(communicator)  ! is this pe part of IO PE set setno. if yes return set communicator, else return MPI_COMM_NULL
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: communicator
    integer :: i
#if defined(STAND_ALONE)
    integer :: pe_mex, pe_mey, pe_me, ierr
    integer, external :: RPN_COMM_mype
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
#endif

    communicator = MPI_COMM_NULL
    if( .not. is_valid_io_setno(setno) ) return   ! invalid IO set
!
    do i = 1 , io_set(setno)%npe
      if(pe_mex == io_set(setno)%x(i) .and. pe_mey == io_set(setno)%y(i)) then  ! I am an IO pe
        communicator =  io_set(setno)%comm
        return
      endif
    enddo
    
  end function io_pe_comm
!
! the caller of this function will be responsible for freeing the returned array memory
!
  function io_pe_coord(setno) result(list)  ! return pointer to x y coordinates in grid of PEs if set is valid, else return a null pointer
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer, dimension(:,:), pointer :: list
    integer :: npes

    nullify(list)
    if( .not. is_valid_io_setno(setno) ) return   ! invalid IO set
!
    npes = io_set(setno)%npe
    allocate(list(npes,2))
    list(1:npes,1) = io_set(setno)%x            ! list(:,1) x coordinates of IO PEs
    list(1:npes,2) = io_set(setno)%y            ! list(:,2) y coordinates of IO PEs
  end function io_pe_coord
!
  function io_pe_call(setno,callback,argv) result(status)
    implicit none
    integer, intent(IN) :: setno      ! set number as returned by create_ioset
    type(C_PTR) :: argv               ! pointer describing arguments to callback
    integer, external :: callback     ! user routine to be called
    integer :: status
    integer, dimension(:), pointer :: argvf

    status = -1
    if( .not. is_valid_io_setno(setno) ) return   ! invalid IO set
!
    status = 0
    if(io_set(setno)%me < 0) return   ! not a member, nothing to do
!
    status = callback(argv,setno,io_set(setno)%comm,io_set(setno)%me,   &
                      io_set(setno)%x,io_set(setno)%y,io_set(setno)%npe)
  end function io_pe_call
!
  function io_pe_size(setno) result(population)  ! return set population if set is valid, else return -1
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: population

    population = -1
    if( .not. is_valid_io_setno(setno) ) return   ! invalid IO set
!
    population = io_set(setno)%npe
  end function io_pe_size
!
  function io_pe_groups(setno) result(ngroups)  ! return number of groups if set is valid, else return -1
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: ngroups

    ngroups = -1
    if( .not. is_valid_io_setno(setno) ) return   ! invalid IO set
!
    ngroups = io_set(setno)%ngroups
  end function io_pe_groups
!
  function free_ioset(setno) result(freed)   ! free table space associated with IO PE set setno
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: freed
    integer :: ierr
!
    freed = -1
    if( .not. is_valid_io_setno(setno) ) return   ! invalid IO set
!
    deallocate(io_set(setno)%x)
    nullify(io_set(setno)%x)
    deallocate(io_set(setno)%y)
    nullify(io_set(setno)%y)
    io_set(setno)%ioset = -1
    io_set(setno)%npe = -1
    if(io_set(setno)%comm .ne. MPI_COMM_NULL) call mpi_comm_free(io_set(setno)%comm,ierr)
    io_set(setno)%comm = MPI_COMM_NULL
    io_set(setno)%me = -1
    io_set(setno)%ngroups = 0
    io_set(setno)%groupsize = 0
    freed = setno
  end function free_ioset
!
  function create_ioset(npes,method) result(setno)  ! pe_indomm AKA "GRID" is assumed as a communicator
    implicit none
    integer, intent(IN) :: npes
    integer, intent(IN) :: method
    integer :: setno
    integer :: ierr, my_color, i, j, newset, ordinal_in_set

#if defined(STAND_ALONE)
    include 'mpif.h'
    integer :: pe_mex, pe_mey, pe_indomm, pe_me, pe_nx, pe_ny
    integer, external :: RPN_COMM_comm, RPN_COMM_mype
    pe_indomm = RPN_COMM_comm('GRID')          ! grid communicator
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
    call mpi_allreduce(pe_mex,pe_nx,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)  ! roundabout way to get pe_nx
    pe_nx = pe_nx + 1
    call mpi_allreduce(pe_mey,pe_ny,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)  ! roundabout way to get pe_ny
    pe_ny = pe_ny + 1
#endif
!
    setno = -1                       ! precondition for failure
    if(iosets >= MAXIOSETS) return   ! OOPS, table is full
    newset = iosets + 1
    do i = 1 , iosets
      if(io_set(i)%ioset == -1) newset = i   ! recycle freed set
    enddo
    iosets = max(iosets,newset)
!
    allocate(io_set(newset)%x(npes))
    allocate(io_set(newset)%y(npes))
    io_set(newset)%ioset = -1
    io_set(newset)%npe = -1
    io_set(newset)%comm = MPI_COMM_NULL
    io_set(newset)%me = -1
    io_set(newset)%megrid = -1
    io_set(newset)%ngroups = 0
    io_set(newset)%groupsize = 0
    call RPN_COMM_make_io_pe_list(io_set(newset)%x, io_set(newset)%y, npes, pe_nx, pe_ny, method)
    if(io_set(newset)%x(1) .ne. -1) then
      if( RPN_COMM_check_ioset(newset,io_set(newset)%x ,io_set(newset)%y, npes, pe_nx, pe_ny, pe_me, .true.) == -1 ) io_set(newset)%x(1) = -1
    endif
    if(io_set(newset)%x(1) == -1) then        ! miserable failure, list of IO pes could not be created
      deallocate(io_set(newset)%x,io_set(newset)%y)
      nullify(io_set(newset)%x)
      nullify(io_set(newset)%y)
!      setno = -1    ! not needed, it is already -1
      return
    endif
    io_set(newset)%ioset = newset
    io_set(newset)%npe = npes
    io_set(newset)%groupsize = min(pe_nx, pe_ny, npes)   ! size of groups for this IO PE set (last group may be shorter)
    io_set(newset)%ngroups = (npes+io_set(newset)%groupsize-1)/(io_set(newset)%groupsize)   ! number of groups in this IO PE set
    my_color = MPI_UNDEFINED
    ordinal_in_set = -1
    do i = 1 , npes
      if(pe_mex == io_set(newset)%x(i) .and. pe_mey == io_set(newset)%y(i)) then  ! I am an IO pe
        my_color = 1
        ordinal_in_set = i
        exit
      endif
    enddo
    call MPI_COMM_SPLIT(pe_indomm,my_color,ordinal_in_set,io_set(newset)%comm,ierr)   ! communicator for this set
    if(io_set(newset)%comm .ne. MPI_COMM_NULL) then                          ! get rank in communicator if part of set
      call mpi_comm_rank(io_set(newset)%comm,io_set(newset)%me,ierr)
      io_set(newset)%megrid = pe_me                                          ! rank in grid
    endif
    setno = newset
#if defined(SELF_TEST)
    if(pe_me == 0) then
      print *,'DEBUG: creating IO set ',newset,' npes=',npes
      print *,'DEBUG: npe, groups, me',io_set(newset)%npe,io_set(newset)%ngroups,io_set(newset)%me
    endif
#endif
  end function create_ioset
!
end module RPN_COMM_io_pe_tables
!
!=========================   embedded self test   ===============================
!
#if defined(STAND_ALONE)
  subroutine RPN_COMM_io_pe_test(pe_nx,pe_ny,pe_me,pe_mex,pe_mey)
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
!  integer, external :: RPN_COMM_create_io_set, RPN_COMM_free_io_set, RPN_COMM_io_pe_callback
!  integer, external :: RPN_COMM_is_io_pe, RPN_COMM_io_pe_size
!  interface
!    function RPN_COMM_io_pe_coord(setno) result(list)
!      implicit none
!      integer, intent(IN) :: setno
!      integer, dimension(:,:), pointer :: list
!    end function RPN_COMM_io_pe_coord
!  end interface
  integer, intent(IN) :: pe_nx,pe_ny,pe_me,pe_mex,pe_mey
#else
  subroutine RPN_COMM_io_pe_test()
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces.inc>
#endif
!****f* rpn_comm/example example code
! EXAMPLE
  integer, external :: RPN_COMM_io_pe_test_callback   ! the callback function returns an integer value
  integer setno,nio,me_io,setno2,setno3,status
  integer, dimension(1), target :: argv               ! the argument list for the callback function
  integer, dimension(:,:), pointer :: iopelist        ! will receive a coordinate list
  integer, dimension(1) :: tbcst                      ! this will be broadcast in the test
!
! IGNORE
  if(pe_me == 0)  print *,'DEBUG: pe_nx,pe_ny',pe_nx,pe_ny
  nio = min(pe_nx,pe_ny)
  print 100,'RPN_COMM_io_pe test program, pe_nx,pe_ny,pe_me,pe_mex,pe_mey,nio=',pe_nx,pe_ny,pe_me,pe_mex,pe_mey,nio
! EXAMPLE
  setno = RPN_COMM_create_io_set(nio+2,0)                   ! create a set of IO PEs containing nio+2 members
  me_io = RPN_COMM_is_io_pe(setno)                          ! is this PE a member of this IO PE set ?
  if(me_io .ne. -1) then
    print *,"I am a proud IO pe !"                          ! YES it is
  else
    print *,"I am a lazy  NON-IO pe !"                      ! NO it is not
  endif
  print *,"set number, size of set='",setno,RPN_COMM_io_pe_size(setno)         ! get size of this IO PE set
  setno2 = RPN_COMM_create_io_set(nio,0)                                       ! crete another set containing nio PEs
  print *,"set number, size of set='",setno2,RPN_COMM_io_pe_size(setno2)       ! get size of this other IO PE set (should be nio+2)
  setno = RPN_COMM_free_io_set(setno)                                          ! delete IO set
  print *,'DEBUG: freed IO set ',setno
  print *,"set number, size of set='",setno,RPN_COMM_io_pe_size(setno)         ! this should return -1
  setno = RPN_COMM_create_io_set(nio,0)                                        ! re create IO PE set, this time with nio PEs
  print *,"set number, size of set='",setno,RPN_COMM_io_pe_size(setno)         ! this should return nio now
  setno3 = RPN_COMM_create_io_set(nio-1,0)                                     ! crete another set containing nio-1 PEs
  print *,"set number, size of set='",setno3,RPN_COMM_io_pe_size(setno3)       ! this should return nio-1
  argv(1) = pe_me                                                              ! argument array for callback function
  status = RPN_COMM_io_pe_callback(setno3,RPN_COMM_io_pe_test_callback,C_LOC(argv(1)))    ! call to callback function(see example below)
  print *,"after callback, status,argv=",status,argv(1)                        ! status 0 from non members of IO set, return of callback function on members
  iopelist => RPN_COMM_io_pe_coord(setno3)                                     ! get grid coordinates of PEs in IO set setno3
  print *,"PE list x=",iopelist(:,1)                                           ! row 1, X coordinates in "GRID"
  print *,"PE list y=",iopelist(:,2)                                           ! row 2, Y coordinates in "GRID"
  tbcst = pe_me
  if(RPN_COMM_is_io_pe(setno3) .ne. -1)then  ! part of the set only, bcst pe_me of highest rank in set
    call RPN_COMM_io_pe_bcast(tbcst,1,'MPI_INTEGER',RPN_COMM_io_pe_size(setno3)-1,setno3,status) 
!   root for broadcast is RPN_COMM_io_pe_size(setno3)-1, get "GRID" ordinal of PE with highest rank (last PE) in IO set
!   root for broadcast 0 would have returned the "GRID" ordinal of first (lowest rank) PE in set
    print *,'tbcst after broadcast',tbcst                                      ! and the answer is ...
  endif
!******
100 format(A,10I5)
  return
end subroutine
!===============================================================================
! beginning of USER CALLABLE routines/functions
!===============================================================================
!
!    status = callback(argv,setno,io_set(setno)%comm,io_set(setno)%me,   &
!                      io_set(setno)%x,io_set(setno)%y,io_set(setno)%npe)
!****f* rpn_comm/example_callback example of callback function
! EXAMPLE
function RPN_COMM_io_pe_test_callback(argv,setno,comm,me,x,y,npe) result(status)  ! demo callback function
! this function will only be called on members of the IO PE set (setno)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: setno,comm,me,npe      ! set number, MPI communicator of this IO PE set, this PE's ordinal in set
  integer, dimension(npe), intent(IN) :: x,y
  type(C_PTR), intent(IN) :: argv               ! "blind pointer" to argument list
!
  integer :: status
  integer, dimension(:), pointer :: argvf        ! what the "blind pointer" will be mapped into
!
  call C_F_POINTER(argv,argvf,[1])                                ! transform C pointer into Fortran pointer
  print *,"IO PE CALLBACK set, me, argument = ",setno,me,argvf(1) ! an integer array of dimension 1 was expected
  print *,"CALLBACK x=",x                                         ! list of X "GRID" coordinates
  print *,"CALLBACK y=",y                                         ! list of Y "GRID" coordinates
  argvf(1) = -(100 +argvf(1))                                     ! we modify the argument list (optional)
  status = -123                                                   ! status to be returned to caller on this PE
  return
end
!******
!
!****f* rpn_comm/RPN_COMM_create_io_set
! SYNOPSIS
function RPN_COMM_create_io_set(npes,method) result(setno)  !InTf!  ! create a set of IO PEs
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                                             !InTf!
! ARGUMENTS
  integer, intent(IN)  :: npes      !InTf!  ! number of IO PEs desired in set ( set size will be min(npes,pe_nx,pe_ny) )
  integer, intent(IN)  :: method    !InTf!  ! method 0 is the only one supported for the time being
  integer :: setno                  !InTf!  ! set number created ( -1 if creation failed )
!******

  setno = create_ioset(npes,method)
  return
end function RPN_COMM_create_io_set  !InTf!  
!
!****f* rpn_comm/RPN_COMM_free_io_set
! SYNOPSIS
function RPN_COMM_free_io_set(setno) result(freed)  !InTf!  ! delete a set of IO PEs created by RPN_COMM_create_io_set
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: freed                    !InTf!  ! setno if operation succeeded, -1 if it failed
!******
!
  freed = free_ioset(setno)
end function RPN_COMM_free_io_set     !InTf!  
!
!****f* rpn_comm/RPN_COMM_is_valid_io_setno
! SYNOPSIS
function RPN_COMM_is_valid_io_setno(setno) result(valid)   !InTf!  ! is this IO set number valid ?
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                      !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno       !InTf!  ! set number as returned by RPN_COMM_create_io_set
  logical :: valid                   !InTf!  ! .true. if set is valid, .false. otherwise
!******
  valid = is_valid_io_setno(setno)
end function RPN_COMM_is_valid_io_setno   !InTf!
!
!****f* rpn_comm/RPN_COMM_is_io_pe
! SYNOPSIS
function RPN_COMM_is_io_pe(setno) result(ordinal)   !InTf!  ! is this PE part of IO set setno ? 
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                      !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno       !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: ordinal                 !InTf!  ! ordinal in IO PE set if a member, -1 if not
!******
  ordinal = is_io_pe(setno)
end function RPN_COMM_is_io_pe        !InTf!  
!
!****f* rpn_comm/RPN_COMM_io_pe_gridid
! SYNOPSIS
function RPN_COMM_io_pe_gridid(setno,n) result(ordinal)   !InTf!  ! ordinal in "grid" of PE n of IO set setno
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                      !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno       !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer, intent(IN) :: n           !InTf!  ! ordinal of PE for which ordinal in grid is sought (ORIGIN 0)
  integer :: ordinal                 !InTf!  ! ordinal in "grid" if n not larger than IO set size, -1 otherwise
!******
!
  ordinal = -1
  if( .not. is_valid_io_setno(setno) )  return  ! invalid IO set
  if(n<0 .or. n >= io_set(setno)%npe) return  ! n >= population of IO set setno or negative
#if ! defined(SELF_TEST)
  ordinal = pe_id(io_set(setno)%x(n+1) , io_set(setno)%y(n+1))
#endif
end function RPN_COMM_io_pe_gridid        !InTf!  
!
!****f* rpn_comm/RPN_COMM_io_pe_comm
! SYNOPSIS
function RPN_COMM_io_pe_comm(setno) result(communicator)  !InTf!   ! get communicator of IO set setno
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: communicator             !InTf!  ! communicator for IO PE set if a member, null communicator otherwise
!******
  communicator = io_pe_comm(setno)
end function RPN_COMM_io_pe_comm  !InTf!  
!
!****f* rpn_comm/RPN_COMM_io_pe_size
! SYNOPSIS
function RPN_COMM_io_pe_size(setno) result(population)   !InTf!  ! get size of IO set setno
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: population               !InTf!  ! population of IO PE set, -1 if set is not valid
!******
  population = io_pe_size(setno)
end function RPN_COMM_io_pe_size      !InTf!  
!
!****f* rpn_comm/RPN_COMM_io_pe_groups
! SYNOPSIS
function RPN_COMM_io_pe_groups(setno) result(ngroups)   !InTf!  ! get number of groups in IO set setno
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: ngroups                  !InTf!  ! ngroups of IO PE set, -1 if set is not valid
!******
  ngroups = io_pe_groups(setno)
end function RPN_COMM_io_pe_groups    !InTf!  
!
! the caller of this function will be responsible for freeing the returned array memory
!
!****f* rpn_comm/RPN_COMM_io_pe_idlist
! SYNOPSIS
function RPN_COMM_io_pe_idlist(setno) result(idlist)   !InTf!  ! get list of PE ordinals in IO set setno
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                               !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno                !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer, dimension(:), pointer :: idlist      !InTf!  ! list of IO PE set, null pointer if set is not valid
!******

  integer :: setsize, i

  setsize = io_pe_size(setno)
  if(setsize > 0) then
    allocate(idlist(setsize))
    do i = 1, setsize
#if defined(SELF_TEST)
      idlist(i) = 0
#else
      idlist(i) = pe_id(io_set(setno)%x(i),io_set(setno)%y(i))    ! list(n) is ordinal in grid of IO PE n
#endif
    enddo
  else
    idlist => NULL()
  endif
end function RPN_COMM_io_pe_idlist                !InTf!  
!
! the caller of this function will be responsible for freeing the returned array memory
!
!****f* rpn_comm/RPN_COMM_io_pe_coord
! SYNOPSIS
function RPN_COMM_io_pe_coord(setno) result(list)   !InTf!  ! get x and y coordinates in grid of PEs in IO set setno
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none                               !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno                !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer, dimension(:,:), pointer :: list    !InTf!  ! list of IO PE set, null pointer if set is not valid
!******
  list => io_pe_coord(setno)                    ! list(:,1) x coordinates, list(:,2) y coordinates of IO PEs
end function RPN_COMM_io_pe_coord              !InTf!  
!
!****f* rpn_comm/RPN_COMM_io_pe_callback
! SYNOPSIS
function RPN_COMM_io_pe_callback(setno,callback,argv) result(status) !InTf!
! Notes
! callback to be executed only on members of an IO PE set (usually called by all PEs on grid)
! if this PE is not a member of the IO set, nothing happens, status = 0
! if the IO set is invalid, nothing happens, status = -1
! if this PE is a member of the set, callback is called and its return value is put into status
!
! the calling sequence of callback is as follows
! status = callback(argv,setno,comm,ordinal,x,y,npe)
! type(C_PTR) :: argv 
! integer :: setno         ! IO PE set number
! integer :: comm          ! communicator used by this IO PE set
! integer :: ordinal       ! ordinal of this PE in the set
! integer, dimension(npe) :: x  ! x grid coordinates for IO PEs in this set
! integer, dimension(npe) :: y  ! y grid coordinates for IO PEs in this set
! integer :: npe           ! number of PEs in this set
!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
!!  import :: C_PTR      !InTf!
  implicit none                       !InTf!
! ARGUMENTS
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by create_ioset
  type(C_PTR) :: argv                 !InTf!  ! "pickled" arguments to callback
  integer, external :: callback       !InTf!  ! user routine to be called
  integer :: status                   !InTf!
!******
!
  status = io_pe_call(setno,callback,argv)
!
end function RPN_COMM_io_pe_callback  !InTf!  
!
! broadcast within an IO PE set
! this routine MUST be called by ALL PEs in the IO PE set
! but MAY be called by other PEs (the call will be ignored by non members of the set
!
!****f* rpn_comm/RPN_COMM_io_pe_bcast
! SYNOPSIS
subroutine RPN_COMM_io_pe_bcast(buffer,count,datatype,root,setno,ierr) ! cannot easily publish interface because of buffer
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none
! ARGUMENTS
  integer, intent(IN) :: setno                ! set number as returned by RPN_COMM_create_io_set
  integer, intent(IN) :: count                ! number of elements of type datatype
  integer, intent(INOUT), dimension(*) :: buffer   ! data to be broadcast
  character (len=*), intent(IN) :: datatype   ! datatype RPN_COMM style
  integer, intent(IN) :: root                 ! root PE in set for broadcast (typically 0)
  integer, intent(OUT) :: ierr                ! status upon completion
!******
  integer :: dtype
!
  ierr = RPN_COMM_ERROR
  if( is_io_pe(setno) < 0) return             ! bad set or not a member of the set
  call MPI_bcast(buffer,count,RPN_COMM_datyp(datatype),root,io_pe_comm(setno),ierr)
end subroutine RPN_COMM_io_pe_bcast
!
! is this IO PE configuration possible ?
!
!****f* rpn_comm/RPN_COMM_io_pe_valid_set
! SYNOPSIS
function RPN_COMM_io_pe_valid_set(x,y,npes,penx,peny,diag,method) result(status)       !InTf!
!  check if this IO PE configuration is possible
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_io_pe_tables
  implicit none
! ARGUMENTS
    integer, intent(IN) :: npes                !InTf!   ! number of PEs in set
    integer, dimension(npes), intent(OUT) :: x !InTf!   ! x coordinates of PEs in set
    integer, dimension(npes), intent(OUT) :: y !InTf!   ! y coordinates of PEs in set
    integer, intent(IN) :: penx, peny          !InTf!   ! number of PEs along x and y in PE grid
    integer, intent(IN) :: method              !InTf!   ! fill method
    logical, intent(IN) :: diag                !InTf!   ! if .true. print IO PE map and diagnostics
    integer :: status                          !InTf!   ! RPN_COMM_OK or RPN_COMM_ERROR
! Notes
!   some error messages are printed in case of error
!   diagnostics include a map of IO PEs
!******

    integer :: setno

    setno = -1
    status = RPN_COMM_ERROR
    call RPN_COMM_make_io_pe_list(x,y,npes,penx,peny,method)
    if(x(1) == -1) return   ! miserable failure at creation
    if( RPN_COMM_check_ioset(0, x ,y, npes, penx, peny, 0, diag) .ne. 0 ) return
    status = RPN_COMM_OK
end function RPN_COMM_io_pe_valid_set !InTf!
!
!=========================    END of user callable functions    ===============================
!
#if defined SELF_TEST
program self_test
  implicit none
  include 'mpif.h'
  integer ierr,pe_me,pe_mex,pe_mey,npes,pe_nx,pe_ny,pe_indomm
  integer, external :: RPN_COMM_mype, RPN_COMM_comm

  call mpi_init(ierr)
  pe_nx = 1
  pe_ny = 1
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  pe_indomm = RPN_COMM_comm('GRID')
  ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  call mpi_allreduce(pe_mex,pe_nx,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)
  pe_nx = pe_nx + 1
  call mpi_allreduce(pe_mey,pe_ny,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)
  pe_ny = pe_ny + 1
  call RPN_COMM_io_pe_test(pe_nx,pe_ny,pe_me,pe_mex,pe_mey)
  call mpi_finalize(ierr)
  stop
  end
!
integer function RPN_COMM_comm(what)
  implicit none
  include 'mpif.h'
  character (len=*) :: what
  RPN_COMM_comm = MPI_COMM_WORLD
  return
end
!
integer function RPN_COMM_datyp(what)
  implicit none
  include 'mpif.h'
  character (len=*) :: what
  RPN_COMM_datyp = MPI_UNDEFINED
  if(trim(what) .ne. 'MPI_INTEGER') return
  RPN_COMM_datyp = MPI_INTEGER
end function RPN_COMM_datyp
!
integer function RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  implicit none
  include 'mpif.h'
  integer, intent(OUT) :: pe_me,pe_mex,pe_mey
  integer :: i
  integer :: ierr,npes

  call mpi_comm_rank(MPI_COMM_WORLD,pe_me,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  RPN_COMM_mype = 0
  pe_mex = npes - 1
  pe_mey = 0
  do i=7,2,-1
    if(mod(npes,i) == 0 .and. npes .ne. i)then
      pe_mey = (pe_me) / i
      pe_mex = mod(pe_me, i)
      exit
    endif
  enddo
end function RPN_COMM_mype
#endif
