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
#define IN_RPN_COMM_window
#if defined(WITH_DOC)
!****P* rpn_comm/windows  (simplified/restricted version of MPI one sided communications package)
! DESCRIPTION
!   This is a simplified and restricted version of MPI-2 one sided communications.
!
!   When creating the "communication window", some attributes are 
!   determined for said window and will not change during its life.
!   1 - the MPI communicator
!   2 - the MPI data type of the elements contained in the window
!   3 - the size of the window (in number of elements)
!   4 - an array large enough to contain these elements
!       (this array will either be supplied by the caller or allocated internally.
!
!   all further operations refer to the window by its "identifier" 
!   Fortran type : type(rpncomm_window)  (include 'RPN_COMM_types.inc')
!
!   the window creation is a COLLECTIVE operation, all members of the communicator group
!   must call RPN_COMM_i_win_create
!
!   "exposing" a window and terminating the "exposure" of a window are also a COLLECTIVE operation
!
!   remote get/put operations, i.e. sending read/write requests targeting the communication
!   window belonging to a remote PE may only happen when the window is "exposed"
!   remote operations are NOT ALLOWED when the window is "not exposed"
!
!   two modes of one sided communication are available.
!   - active mode
!      a) for the whole communicator group  (fence)
!      b) for a subset of the communicator group (scales better when said group is large)
!         (start/complete/post/wait)
!   - passive mode
!   this mode is selected when calling the routine that "exposes" the communication window
!
!   local get/put operations, i.e. reading/writing from/into the window
!   belonging the local PE may only happen when the window is "not exposed"
!   local operations are NOT ALLOWED whe the window is "exposed" as the result of
!   such an operation would be "undefined"
!
!   a typical sequence of operations would be
!   1 - create a one sided communication window (COLLECTIVE operation)
!    repeat as needed
!    {
!      2a - modify the local copy of the window (if and as needed)
!      2b - "expose" the window and select active or passive mode (COLLECTIVE operation)
!      2c - perform get/put operations targeting remote PEs (as needed on each PE)
!      2d - "end exposing" the window (COLLECTIVE operation)
!      2e - get modifications from the local copy of the window (if and as needed)
!    }
!   3 - free the one sided communication window (COLLECTIVE operation)
!
!   window creation/destruction : RPN_COMM_i_win_create, RPN_COMM_i_win_free
!   window status queries       : RPN_COMM_i_valid_win, RPN_COMM_i_win_check
!   window info queries         : RPN_COMM_i_win_get_ptr, RPN_COMM_i_win_get_size
!   window "exposition"         : RPN_COMM_i_win_open, RPN_COMM_i_win_close
!   window get operations       : RPN_COMM_i_win_get_r, RPN_COMM_i_win_get_l
!   window put operations       : RPN_COMM_i_win_put_r, RPN_COMM_i_win_put_l
!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
!******
#endif
!****iP* rpn_comm/windows  one sided communication window management module
! SYNOPSIS
module RPN_COMM_windows
!===============================================================================
! one sided communication window management code 
! (used internally by user callable routines)
!===============================================================================
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'
  include 'RPN_COMM_types_int.inc'
  include 'RPN_COMM_constants.inc'
!******
  integer, parameter :: RPN_COMM_MAX_WINDOWS = 64
  integer, save :: integer_size = 0                                       ! will contain extent of MPI_INTEGER
  type(rpncomm_windef), dimension(:), pointer, save :: win_tab => NULL()  ! rpn comm window table
  logical, save :: debug_mode = .true.
  type(rpncomm_operator), save :: op = NULL_rpncomm_operator
 contains
!****if* RPN_COMM_windows/init_win_tab
! SYNOPSIS
  subroutine init_win_tab
!===============================================================================
! allocate and initialize internal single sided communication windows table
!===============================================================================
!******
    implicit none
    integer :: i
    integer(kind=MPI_ADDRESS_KIND) :: lb

    if(ASSOCIATED(win_tab)) return                  ! already initialized, nothing to do
    allocate(win_tab(RPN_COMM_MAX_WINDOWS))         ! allocate for up to RPN_COMM_MAX_WINDOWS windows
    do i = 1 , RPN_COMM_MAX_WINDOWS
      win_tab(i) = NULL_rpncomm_windef              ! null entry (see RPN_COMM_types_int.inc)
      win_tab(i)%com = MPI_COMM_NULL                ! with MPI null communicator
      win_tab(i)%typ = MPI_DATATYPE_NULL            ! and MPI null datatype
      win_tab(i)%win = MPI_WIN_NULL                 ! MPI null window
    enddo
    call MPI_type_get_extent(MPI_INTEGER, lb, integer_size, i) ! get size(extent) of MPI_INTEGER
  end subroutine init_win_tab

!****if* RPN_COMM_windows/win_ptr
! SYNOPSIS
  function valid_win_entry(win_ptr) result(is_valid)
!===============================================================================
! check if win_ptr (C pointer pointing into win_tab) points to a valid entry
!===============================================================================
!******
    implicit none
    type(C_PTR), intent(IN), value :: win_ptr        ! pointer to entry in win_tab
    logical :: is_valid

    type(rpncomm_windef), pointer :: win_entry
    type(C_PTR) :: temp

    is_valid = .false.
    if(.not. c_associated(win_ptr)) return           ! win_ptr C pointer is not valid
    if(.not. associated(win_tab)) return             ! win_tab does not exist yet

    call c_f_pointer(win_ptr,win_entry)              ! make Fortran pointer to win_entry from win_ptr

    if(win_entry%indx <=0 .or. win_entry%indx>RPN_COMM_MAX_WINDOWS) return  ! invalid index into window table
    temp = c_loc(win_tab(win_entry%indx))            ! address pointed to by win_entry index component
    if(.not. c_associated(temp,win_ptr)) return      ! component index should point to win_entry itself
    if(.not. c_associated(win_entry%base)) return    ! no array associated with entry
    if(win_entry%com == MPI_COMM_NULL) return        ! no valid communicator
    if(win_entry%typ == MPI_DATATYPE_NULL) return    ! no valid datatype
    if(win_entry%win == MPI_WIN_NULL) return         ! no valid window
    if(win_entry%siz <= 0) return                    ! invalid size (MUST be > 0)

    is_valid = .true.                                ! no more reason to think entry is not valid

  end function valid_win_entry

!****if* RPN_COMM_windows/create_win_entry
! SYNOPSIS
  subroutine create_win_entry(base,typ,siz,comm,indx,ierr)
!===============================================================================
! create a new one sided communication window
!
! if argument base is a NULL pointer (C_NULL_PTR), allocate an internal array
! using the recommended MPI_alloc_mem MPI routine
!===============================================================================
!******
    implicit none
    type(C_PTR), intent(IN), value :: base ! base address of array exposed through window (may be C_NULL_PTR)
    integer, intent(IN) :: typ             ! MPI data type (as defined in mpif.h)
    integer, intent(IN) :: siz             ! number of elements in array associated to this window
    integer, intent(IN) :: comm            ! MPI communicator for this one sided window
    integer, intent(OUT) :: indx           ! index into wintab of created window (used for consistency test)
    integer, intent(OUT) :: ierr           ! return status (RPN_COMM_ERROR or RPN_COMM_OK)

    integer :: i, extent, ierror, group
    integer(kind=MPI_ADDRESS_KIND) win_size  ! may be wider than a default integer
    integer(kind=MPI_ADDRESS_KIND) :: lb

    if(.not. associated(win_tab)) call init_win_tab  ! create and initialize window table if necessary
    ierr = RPN_COMM_ERROR                  ! preset for failure
    indx = -1

    call MPI_comm_group(comm,group,ierror)
    if(ierror .ne. MPI_SUCCESS) return     ! error getting group associated with communicator (bad communicator ?)

    call MPI_type_get_extent(typ, lb, extent, ierror)  ! determine size associated with MPI datatype
    if(ierror .ne. MPI_SUCCESS) return         !  invalid type ? other error ?
    if( mod(extent,integer_size) .ne. 0 ) return    ! extent of data type must be a multiple of integer size

    win_size = extent * siz    ! possibly wider integer needed by MPI_win_create

    do i = 1 , RPN_COMM_MAX_WINDOWS        ! loop over window table entries (we are looking for a free entry)
      if(.not. c_associated(win_tab(i)%base) ) then  ! this entry is free
        if(c_associated(base)) then    ! if user has already allocated space
          if(debug_mode) print *,'DEBUG: user supplied space, size=',win_size
          win_tab(i)%base = base
          win_tab(i)%is_user = .true.  ! user supplied space
        else                           ! otherwise allocate needed space with MPI_alloc_mem
          if(debug_mode) print *,'DEBUG: allocating with MPI_alloc_mem, size=',win_size
          call MPI_alloc_mem(win_size, MPI_INFO_NULL, win_tab(i)%base,ierr)
          if(ierr .ne. MPI_SUCCESS) return       ! could not allocate memory, OUCH !!
          win_tab(i)%is_user = .false.  ! internally supplied space
        endif
        call c_f_pointer(win_tab(i)%base,win_tab(i)%remote,[extent*siz])  ! make Fortran pointer to array associated with window

        call MPI_win_create(win_tab(i)%remote, win_size, extent, MPI_INFO_NULL, comm, win_tab(i)%win, ierror) ! create window
        if(ierror .ne. MPI_SUCCESS) then
          print *,'ERROR: window creation unsuccessful'
          return
        endif
        indx = i
        win_tab(i)%is_open = .false.             ! window is not "exposed"
        win_tab(i)%opr = MPI_OP_NULL             ! MPI accumulate operator associated to window
        win_tab(i)%typ = typ                     ! MPI data type associated to window
        win_tab(i)%ext = extent/integer_size     ! size (extent) of MPI data type in MPI_INTEGER units
        win_tab(i)%indx = indx                   ! entry index points to itself
        win_tab(i)%siz = siz                     ! number of elements in window
        win_tab(i)%com = comm                    ! communicator associated with window
        win_tab(i)%grp = group                   ! group associated with said communicator
        win_tab(i)%active_mode = .true.          ! use active mode by default
        win_tab(i)%s_group = MPI_GROUP_NULL      ! no group by default for remote targets
        win_tab(i)%r_group = MPI_GROUP_NULL      ! no group by default for exposure to remote accesses
        ierr = RPN_COMM_OK                       ! all is well
        return
      endif
    enddo
    return  ! if we fall through the loop, the table was full, we will return RPN_COMM_ERROR
  end subroutine create_win_entry

!****if* RPN_COMM_windows/win_valid
! SYNOPSIS
  function win_valid(window) result(is_valid)
!===============================================================================
! check if window is indeed a valid one
!===============================================================================
!******
  implicit none
  type(rpncomm_window), intent(IN) :: window
  logical :: is_valid
  type(C_PTR) :: temp

  is_valid = .false.
  if(.not. associated(win_tab)) call init_win_tab  ! create and initialize window table if necessary

  if( ieor(window%t1,RPN_COMM_MAGIC) .ne. window%t2 )   return  ! inconsistent tags
  if(window%t2 <= 0 .or. window%t2>RPN_COMM_MAX_WINDOWS) return  ! invalid index

  temp = c_loc(win_tab(window%t2))                 ! address of relevant win_tab entry
  if( .not. c_associated(window%p,temp)) return    ! window%p must point to win_tab(window%t2)

  is_valid = valid_win_entry(window%p)             ! check that window description itself in table is valid
  return

end function win_valid

end module RPN_COMM_windows

!===============================================================================
! test code for one sided communication window package
!===============================================================================
subroutine RPN_COMM_i_win_test(nparams,params)
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
  include 'mpif.h'
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params

  integer, parameter :: DATA_SIZE = 1024
  type(rpncomm_window) :: window, window2
  type(rpncomm_datatype) :: dtype
  type(rpncomm_communicator) :: com
  type(rpncomm_operator) :: oper
  integer :: siz
  type(C_PTR) :: array, array2, cptr, cptr2
  integer, dimension(:), pointer :: fptr, fptr2
  integer, dimension(DATA_SIZE), target :: my_data
  integer :: ierr, i, nerrors, nval, offset, offset_r, expected, got
  integer :: me, me_x, me_y, status, wsiz, npes, target_pe, from_pe
  integer, dimension(100), target :: local, local2
  integer, dimension(:), pointer :: errors

!  debug_mode = .true.
  siz = DATA_SIZE
  status = RPN_COMM_mype(Me,Me_x,Me_y)
  call RPN_COMM_size( RPN_COMM_GRID, npes ,ierr )
  print *,'TEST INFO: this is PE',me+1,' of',npes
  allocate(errors(npes))
  call RPN_COMM_i_comm(RPN_COMM_GRID,com)
  call RPN_COMM_i_datyp(RPN_COMM_INTEGER,dtype)

  do i = 1 , siz
    my_data(i) = -i
  enddo
! ================================== create windows ======================================
  array = c_loc(my_data(1))
  call RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: created window, PE=',me
  else
    print *,'TEST ERROR: cannot create window, PE=',me
    return
  endif

  array2 = C_NULL_PTR
  call RPN_COMM_i_win_create(window2,dtype,siz*2,com,array2,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: created window2, PE=',me
  else
    print *,'TEST ERROR: cannot create window2, PE=',me
    return
  endif
! ================================== get windows properties ======================================
  wsiz = RPN_COMM_i_win_get_size(window,ierr)  ! get size of data array associated to window
  if(me == 0) then
    print *,'TEST INFO: size associated to window =',wsiz
  endif
  cptr = RPN_COMM_i_win_get_ptr(window,ierr)   ! get C pointer to local window data
  call c_f_pointer(cptr,fptr,[wsiz])           ! convert to Fortran pointer
  nerrors = 0
  do i = 1 , wsiz                              ! check values in array (consistency check)
    if(fptr(i) .ne. my_data(i)) nerrors = nerrors + 1
  enddo
  print *,'TEST INFO: nerrors for window =',nerrors
  if(nerrors > 0) return

  wsiz = RPN_COMM_i_win_get_size(window2,ierr)  ! get size of data array associated to window
  if(me == 0) then
    print *,'TEST INFO: size associated to window2 =',wsiz
  endif
  cptr2 = RPN_COMM_i_win_get_ptr(window2,ierr)   ! get C pointer to local window data
  call c_f_pointer(cptr2,fptr2,[wsiz])           ! convert to Fortran pointer
  fptr2 = -1                                     ! initialize contents of array associated with windows2
! ===================================== remote PUT test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
  print *,'TEST INFO: target PE is',target_pe,' accessing PE is',from_pe
! ===================================== active mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: PUT active fence mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      print *,'i=',i,' got',my_data(i),' expected',local(i-offset)
      if( my_data(i) .ne. local(i-offset))  nerrors = nerrors + 1
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' error'
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: PUT active fence mode test end'
!  endif
! ===================================== active group mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 35                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: PUT active group mode test start'
!  endif

  call RPN_COMM_i_win_group(window,[target_pe],[from_pe],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: group creation OK'
  else
    print *,'TEST ERROR: group creation failure'
  endif
  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [12,14,16,18,20] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)

  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      print *,'i=',i,' got',my_data(i),' expected',local(i-offset)
      if( my_data(i) .ne. local(i-offset))  nerrors = nerrors + 1
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' error'
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from PUT',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST ERROR:',nerrors,' unexpected values found from PUT'
!   endif

  call RPN_COMM_i_win_group(window,[-1],[-1],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: window group deletion OK'
  else
    print *,'TEST ERROR: window group deletion failure'
  endif

!  if(me == 0) then
    print *,'TEST INFO: PUT active group mode test end'
!  endif
! ================================== passive mode ======================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 17                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: PUT passive mode test start'
!  endif

  target_pe = mod(me+1,npes)

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [1,3,5,7,9] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      got = my_data(i)
      expected = local(i-offset)
      if( got .ne. expected)  then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' ERROR',nerrors
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from PUT',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST ERROR:',nerrors,' unexpected values found from PUT'
!   endif

!  if(me == 0) then
    print *,'TEST INFO: PUT passive mode test end'
!  endif
! ===================================== remote GET test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
! ===================================== active mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET active fence mode test start'
!  endif

!  local = -(me+1)
  local(1:5) = [2,4,6,8,10]
  do i = 1,nval
    my_data(i+offset) = local(i) + 100*me
  enddo
  local2 = -1
  print *,'TEST INFO: my_data',my_data(offset+1:offset+nval)

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  do i = 1 , nval
      expected = my_data(i+offset)-100*me+100*target_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST INFO:',nerrors,' local unexpected values found from GET'
!   endif

!   if(me == 0) then
    print *,'TEST INFO: GET active fence mode test end'
!   endif
! ===================================== active group mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET active group mode test start'
!  endif

!  local = -(me+1)
  local(1:5) = [12,14,16,18,20]
  do i = 1,nval
    my_data(i+offset) = local(i) + 100*me
  enddo
  local2 = -1
  print *,'TEST INFO: my_data',my_data(offset+1:offset+nval)

  call RPN_COMM_i_win_group(window,[target_pe],[from_pe],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: group creation OK'
  else
    print *,'TEST ERROR: group creation failure'
  endif
  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  do i = 1 , nval
      expected = my_data(i+offset)-100*me+100*target_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST INFO:',nerrors,' local unexpected values found from GET'
!   endif

  call RPN_COMM_i_win_group(window,[-1],[-1],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: window group deletion OK'
  else
    print *,'TEST ERROR: window group deletion failure'
  endif

!   if(me == 0) then
    print *,'TEST INFO: GET active group mode test end'
!   endif
! ===================================== passive mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET passive mode test start'
!  endif

!  local = -(me+1)
  local(1:5) = [1,3,5,7,9]
  do i = 1,nval
    my_data(i+offset) = local(i) + 100*me
  enddo
  local2 = -1
  print *,'TEST INFO: my_data',my_data(offset+1:offset+nval)

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  do i = 1 , nval
      expected = my_data(i+offset)-100*me+100*target_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST INFO:',nerrors,' local unexpected values found from GET'
!   endif

!   if(me == 0) then
    print *,'TEST INFO: GET passive mode test end'
!   endif
! ===================================== remote GET/PUT test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
  print *,'TEST INFO: target PE is',target_pe,' accessing PE is',from_pe
! ===================================== active mode ===================================
  do i = 1 , siz
    my_data(i) = -i - 100*me
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
  offset_r = 20
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET/PUT active fence mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),from_pe,offset_r,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),from_pe,offset_r+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      print *,'i=',i,' got',my_data(i),' expected',local(i-offset)
      if( my_data(i) .ne. local(i-offset))  nerrors = nerrors + 1
    else
      if(my_data(i) .ne. -i - 100*me) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i - 100*me,' error'
      endif
    endif
  enddo
  do i = 1 , nval
      expected = my_data(i+offset_r) + 100*me - 100*from_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET/PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: GET/PUT active fence mode test end'
!  endif
! ===================================== group mode ===================================
  do i = 1 , siz
    my_data(i) = -i - 100*me
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
  offset_r = 20
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET/PUT active group mode test start'
!  endif

  call RPN_COMM_i_win_group(window,[target_pe,from_pe],[target_pe,from_pe],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: group creation OK'
  else
    print *,'TEST ERROR: group creation failure'
  endif
  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),from_pe,offset_r,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),from_pe,offset_r+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      print *,'i=',i,' got',my_data(i),' expected',local(i-offset)
      if( my_data(i) .ne. local(i-offset))  nerrors = nerrors + 1
    else
      if(my_data(i) .ne. -i - 100*me) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i - 100*me,' error'
      endif
    endif
  enddo
  do i = 1 , nval
      expected = my_data(i+offset_r) + 100*me - 100*from_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET/PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

  call RPN_COMM_i_win_group(window,[-1],[-1],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: window group deletion OK'
  else
    print *,'TEST ERROR: window group deletion failure'
  endif
!  if(me == 0) then
    print *,'TEST INFO: GET/PUT active group mode test end'
!  endif
! ===================================== passive mode ===================================
  do i = 1 , siz
    my_data(i) = -i - 100*me
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
  offset_r = 20
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET/PUT passive mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),from_pe,offset_r,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),from_pe,offset_r+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      print *,'i=',i,' got',my_data(i),' expected',local(i-offset)
      if( my_data(i) .ne. local(i-offset))  nerrors = nerrors + 1
    else
      if(my_data(i) .ne. -i - 100*me) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i - 100*me,' error'
      endif
    endif
  enddo
  do i = 1 , nval
      expected = my_data(i+offset_r) + 100*me - 100*from_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET/PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: GET/PUT passive mode test end'
!  endif
! ===================================== remote ACC test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
  print *,'TEST INFO: target PE is',target_pe,' accessing PE is',from_pe
! ================================== set windows properties ======================================
  call RPN_COMM_i_oper(RPN_COMM_MAX,oper)
  call RPN_COMM_i_win_oper(window,oper,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: associated operator to window, PE=',me
  else
    print *,'TEST ERROR: cannot associate operator to window, PE=',me
    goto 777
  endif
! ===================================== active mode ===================================
! initial values for local window are negative, the accumulate operator used is MPI_MAX
! this is then equivalent to a straight PUT
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: ACC active fence mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_acc_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),oper,ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_acc_r(window,c_loc(local(i)),target_pe,offset+i-1,1,oper,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      print *,'i=',i,' got',my_data(i),' expected',local(i-offset)
      if( my_data(i) .ne. local(i-offset))  nerrors = nerrors + 1
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' error'
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from ACC',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from ACC'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: ACC active fence mode test end'
!  endif
! ================================== set windows properties ======================================
  print *,'=============================================================================='
  call RPN_COMM_i_oper(RPN_COMM_BXOR,oper)
  print *,'TEST INFO: operator=',oper%t2
  call RPN_COMM_i_win_oper(window,oper,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: associated operator to window, PE=',me
  else
    print *,'TEST ERROR: cannot associate operator to window, PE=',me
    goto 777
  endif
! ===================================== passive mode ===================================
! initial values for local window are negative, the accumulate operator used is MPI_MAX
! this is then equivalent to a straight PUT
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'TEST INFO: ACC passive mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    return
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_acc_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),oper,ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_acc_r(window,c_loc(local(i)),target_pe,offset+i-1,1,oper,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    return
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    got = my_data(i)
    if(i>offset .and. i<=offset+nval) then
!      expected = local(i-offset)
      expected = -i
      expected = ieor(local(i-offset),expected)
      if( got .ne. expected)  then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,'  ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,'  OK'
      endif
    else
      expected = -i
      if(got .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      endif
    endif
  enddo
  if(nerrors .ne. 0) print *,'TEST INFO:',nerrors,' unexpected local values found from ACC'
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from ACC',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from ACC'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: ACC passive mode test end'
!  endif
! ================================== free windows ======================================
777 continue
  call RPN_COMM_i_win_free(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: freed window, PE=',me
  else
    print *,'TEST ERROR: cannot free window, PE=',me
    return
  endif

  call RPN_COMM_i_win_free(window2,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: freed window2, PE=',me
  else
    print *,'TEST ERROR: cannot free window2, PE=',me
    return
  endif
  return
end subroutine RPN_COMM_i_win_test
!===============================================================================
! beginning of USER CALLABLE routines/functions
!===============================================================================
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_group control membership of get and put groups for this PE and this window
! SYNOPSIS
subroutine RPN_COMM_i_win_group(window,pes_to,pes_from,ierr)  !InTf!
!===============================================================================
! for this one sided communication window,
! provide a list of PEs into which remote operations (get/put/acc) will be performed
! provide a list of PEs from which remote operations (get/put/acc) will be accepted
!===============================================================================
!
! window   (IN)      rpn_comm window (see RPN_COMM_types.inc)
! pes_to   (IN)      array containing the list of pe_s who's window will be 
!                    the target of remote operations from this PE (get/put/acc)
! pes_from (IN)      array containing the list of pe_s that will be accessing
!                    via remote operations this PE's window (get/put/acc)
! ierr     (OUT)     will be set to RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: rpncomm_window                                       !InTf!
! ARGUMENTS
  type(rpncomm_window), intent(IN) :: window                       !InTf!
  integer, dimension(:), intent(IN) :: pes_to                      !InTf!
  integer, dimension(:), intent(IN) :: pes_from                    !InTf!
  integer, intent(OUT) :: ierr                                     !InTf!
!******

  integer :: ierr2, index
  type(rpncomm_windef), pointer :: win_entry
  integer :: n_to, n_from

  ierr = RPN_COMM_ERROR
  if( .not. win_valid(window) )  return  ! bad window reference
  n_to = size(pes_to)
  n_from = size(pes_from)
  if(n_to <= 0 .or. n_from <= 0) return  ! bad list of PEs

!  indx = window%t2                   ! window table entry for this window is win_tab(indx)
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(win_entry%is_open) return           ! no change allowed if window is "exposed"
  if(win_entry%s_group .ne. MPI_GROUP_NULL) then   ! free group, it will be re-creataed
    call MPI_group_free(win_entry%s_group,ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif
  if(win_entry%r_group .ne. MPI_GROUP_NULL) then   ! free group, it will be re-creataed
    call MPI_group_free(win_entry%r_group,ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif

  if(pes_to(1)==-1 .and. pes_from(1)==-1) then     ! cancel groups are desired 
    win_entry%s_group = MPI_GROUP_NULL
    win_entry%r_group = MPI_GROUP_NULL
    ierr = RPN_COMM_OK
    return
  endif

  call MPI_group_incl(win_entry%grp, n_to, pes_to, win_entry%s_group, ierr2)
  if(ierr2 .ne. MPI_SUCCESS) return                ! group include failed

  call MPI_group_incl(win_entry%grp, n_from, pes_from, win_entry%r_group, ierr2)
  if(ierr2 .ne. MPI_SUCCESS) return                ! group include failed

  ierr = RPN_COMM_OK
  return
end subroutine RPN_COMM_i_win_group                                   !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_oper
! SYNOPSIS
subroutine RPN_COMM_i_win_oper(window,oper,ierr)  !InTf!
!===============================================================================
! for this one sided communication window,
! associate an operator to accumulate operations
!===============================================================================
!
! window   (IN)      rpn_comm window (see RPN_COMM_types.inc)
! oper     (IN)      rpncomm_operator (see RPN_COMM_types.inc)
! ierr     (OUT)     will be set to RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
#include <RPN_COMM_interfaces_int.inc>
!!  import :: rpncomm_window                                       !InTf!
!!  import :: rpncomm_operator                                     !InTf!
! ARGUMENTS
  type(rpncomm_window), intent(IN) :: window                       !InTf!
  type(rpncomm_operator), intent(IN) :: oper                       !InTf!
  integer, intent(OUT) :: ierr                                     !InTf!
!******

  integer :: ierr2, index
  type(rpncomm_windef), pointer :: win_entry
  integer :: n_to, n_from

  ierr = RPN_COMM_ERROR
  if(.not. RPN_COMM_i_valid_oper(oper) ) then
    print *,'ERROR: invalid operator in RPN_COMM_i_win_oper'
    return  ! bad operator
  endif
  if( .not. win_valid(window) )  then
    print *,'ERROR: invalid window in RPN_COMM_i_win_oper'
    return  ! bad window reference
  endif

  call c_f_pointer( window%p , win_entry )   ! pointer to win_tab entry (rpncomm_windef)
  win_entry%opr = oper%t2

  ierr = RPN_COMM_OK
  return
end subroutine RPN_COMM_i_win_oper                                   !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_create create a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)  !InTf!
!===============================================================================
! create a one sided communication window (user exposed interface)
!===============================================================================
!
! window (OUT)     rpn_comm window type returned to caller (see RPN_COMM_types.inc)
! dtype  (IN)      rpn_comm datatype descriptor (see RPN_COMM_types.inc)
! siz    (IN)      number of elements of type dtype in window
! com    (IN)      rpn_comm communicator used for window (see RPN_COMM_types.inc)
! array  (IN)      C pointer to array associated with window
!                  if defined (not equal to C_NULL_PTR), this user array will be used
!                  if not defined (equal to C_NULL_PTR), an internal array will be allocated and used
! ierr   (OUT)     RPN_COMM_OK or RPN_COMM_ERROR will be returned
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window, rpncomm_datatype, rpncomm_communicator  !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(OUT) :: window                         !InTf!
  type(rpncomm_datatype), intent(IN) :: dtype                         !InTf!
  integer, intent(IN) :: siz                                          !InTf!
  type(rpncomm_communicator), intent(IN) :: com                       !InTf!
  type(C_PTR), intent(IN), value :: array                             !InTf!
! IGNORE
!!! VOID$ array                                                         !InTf!
!******
  integer :: indx

  ierr = RPN_COMM_ERROR
  window = NULL_rpncomm_window

  call create_win_entry(array,dtype%t2,siz,com%t2,indx,ierr)
  if(ierr .ne. RPN_COMM_OK) return
!print *,'DEBUG: window created indx =',indx
!print *,'DEBUG: window created win =',win_tab(indx)%win
  window%p = c_loc(win_tab(indx))        ! point to entry in window table
  window%t1 = ieor(indx,RPN_COMM_MAGIC)  ! xor with magic token
  window%t2 = indx                       ! index into table

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_create                                  !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_free free a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_free(window,ierr)                           !InTf!
!===============================================================================
! delete a previously created one sided communication window (see RPN_COMM_i_win_create)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(INOUT) :: window                       !InTf!
!******

  integer :: indx, ierr1, ierr2

  ierr = RPN_COMM_ERROR
  if(.not. win_valid(window) ) return
  indx = window%t2

!print *,'DEBUG: window is valid, index =',indx
!print *,'DEBUG: window is valid, win   =',win_tab(indx)%win
  call MPI_win_free(win_tab(indx)%win,ierr2)       ! free window
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: error freeing window'
  endif
  win_tab(indx)%win = MPI_WIN_NULL                 ! MPI null window

  if(.not. win_tab(indx)%is_user) then   ! internal storage allocation, release it
    print *,'INFO: attempting to free internal window storage'
    call MPI_free_mem(win_tab(indx)%remote,ierr1)
    if(ierr1 .ne. MPI_SUCCESS) then
      print *,'ERROR: error freeing window storage'
    endif
  else
    ierr1 = MPI_SUCCESS
  endif
  win_tab(indx) = NULL_rpncomm_windef              ! blank entry
  win_tab(indx)%com = MPI_COMM_NULL                ! MPI null communicator
  win_tab(indx)%typ = MPI_DATATYPE_NULL            ! MPI null datatype

  window = NULL_rpncomm_window
  if(ierr1 == MPI_SUCCESS .and. ierr2 == MPI_SUCCESS) then
    ierr = RPN_COMM_OK
!    print *,'DEBUG: successfully freed window'
  endif

  return
end subroutine RPN_COMM_i_win_free                                    !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_open open a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_open(window,active,ierr)                           !InTf!
!===============================================================================
! "expose" a one sided communication window (see RPN_COMM_i_win_create)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! active(IN)      .true. : use active mode, .false.: use passive mode
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical, intent(IN) :: active                                       !InTf!
!******

  integer :: ierr1, ierr2, indx
  logical ::is_open
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  indx = window%t2                   ! window table entry for this window
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if(is_open .or. ierr2 .ne. RPN_COMM_OK) return    ! ERROR: window is already open (exposed) or not valid

  win_tab(indx)%active_mode = active
  if(win_tab(indx)%active_mode) then           ! active mode
    if(win_tab(indx)%s_group == MPI_GROUP_NULL .and. win_tab(indx)%r_group == MPI_GROUP_NULL) then  ! fence mode
      ierr1 = MPI_SUCCESS
!      call MPI_win_fence(MPI_MODE_NOPRECEDE,win_tab(indx)%win,ierr2)
      call MPI_win_fence(0,win_tab(indx)%win,ierr2)
      if(debug_mode) print *,'DEBUG: opening window, fence mode',win_tab(indx)%win,win_tab(indx)%grp
!      call MPI_win_start(win_tab(indx)%grp,0,win_tab(indx)%win,ierr1)    ! start RMA access epoch on win to members of s_group
!      call MPI_win_post (win_tab(indx)%grp,0,win_tab(indx)%win,ierr2)    ! start RMA exposure epoch on win from members of r_group
    else                                                                                         ! start/complete/post/wait mode
      call MPI_win_start(win_tab(indx)%s_group,0,win_tab(indx)%win,ierr1)    ! start RMA access epoch on win to members of s_group
      call MPI_win_post (win_tab(indx)%r_group,0,win_tab(indx)%win,ierr2)    ! start RMA exposure epoch on win from members of r_group
      if(debug_mode) print *,'DEBUG: opening window, groups mode',win_tab(indx)%win,win_tab(indx)%s_group,win_tab(indx)%r_group
    endif
  else                                         ! nothing to do in passive mode
    if(debug_mode) print *,'DEBUG: opening window, passive mode',win_tab(indx)%win
    ierr1 = MPI_SUCCESS
    ierr2 = MPI_SUCCESS
  endif

  if(ierr1 .eq. MPI_SUCCESS .and. ierr2 .eq. MPI_SUCCESS) then
    ierr = RPN_COMM_OK
    win_tab(indx)%is_open = .true.   ! set window open flag
  endif
  return

end subroutine RPN_COMM_i_win_open                                    !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_close close a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_close(window,ierr)                          !InTf!
!===============================================================================
! stop "exposing" a one sided communication window (see RPN_COMM_i_win_create)
! the result of all remotely performed get/put operations may now be used
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
!******

  integer :: ierr1, ierr2, indx
  logical :: is_not_open
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  indx = window%t2                   ! window table entry for this window
  is_not_open = .not. RPN_COMM_i_win_check(window,ierr2)
  if(is_not_open .or. ierr2 .ne. RPN_COMM_OK) return    ! ERROR: window is not open (exposed) or not valid

  if(win_tab(indx)%active_mode) then           ! active mode
    if( win_tab(indx)%s_group == MPI_GROUP_NULL .and. win_tab(indx)%r_group == MPI_GROUP_NULL ) then  ! fence mode
      ierr1 = MPI_SUCCESS
!      call MPI_win_fence(MPI_MODE_NOSUCCEED,win_tab(indx)%win,ierr2)
      call MPI_win_fence(0,win_tab(indx)%win,ierr2)
      if(debug_mode) print *,'DEBUG: fence mode closing window',win_tab(indx)%win
!      call MPI_win_complete(win_tab(indx)%win, ierr1)    ! Complete RMA access epoch on win started MPI_WIN_START
!      call MPI_win_wait    (win_tab(indx)%win, ierr2)    ! Complete RMA exposure epoch started by MPI_WIN_POST on win
    else                                                                                            ! start/complete/post/wait mode
      call MPI_win_complete(win_tab(indx)%win, ierr1)    ! Complete RMA access epoch on win started MPI_WIN_START
      call MPI_win_wait    (win_tab(indx)%win, ierr2)    ! Complete RMA exposure epoch started by MPI_WIN_POST on win
      if(debug_mode) print *,'DEBUG: group mode closing window',win_tab(indx)%win
    endif
  else                                         ! nothing to do in passive mode
    if(debug_mode) print *,'DEBUG: passive mode closing window',win_tab(indx)%win
    ierr1 = MPI_SUCCESS
!    ierr2 = MPI_SUCCESS
    call MPI_barrier(win_tab(indx)%com,ierr2)   ! AHEM!! seems to be needed for openmpi 1.6.5
  endif

  if(ierr1 .eq. MPI_SUCCESS .and. ierr2 .eq. MPI_SUCCESS) then
    ierr = RPN_COMM_OK
    win_tab(indx)%is_open = .false.   ! unset window open flag
  endif
  return

end subroutine RPN_COMM_i_win_close                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_valid check if a one sided communication window is valid
! SYNOPSIS
function RPN_COMM_i_valid_win(window,ierr) result(is_valid)           !InTf!
!===============================================================================
! find if a one sided communication window description is valid
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : .true. (window description is valid) or .false. (not valid)
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_valid                                                 !InTf!
!******
  type(C_PTR) :: temp

  ierr = RPN_COMM_ERROR
  is_valid = win_valid(window)
  if(.not. is_valid ) return             ! check that window description is valid
  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_valid_win                                     !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_check check if a one sided communication window is "exposed"
! SYNOPSIS
function RPN_COMM_i_win_check(window,ierr) result(is_open)            !InTf!
!===============================================================================
! check if a one sided communication window (see RPN_COMM_i_win_create) is "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : .true. (window "exposed") or .false. (window not "exposed")
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_open                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  is_open = .false.
  if(.not. win_valid(window) ) return

  indx = window%t2                   ! window table entry for this window
  is_open = win_tab(indx)%is_open    ! get open (exposed) flag from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_check                                     !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_ptr get data pointer associated to a one sided communication window
! SYNOPSIS
function RPN_COMM_i_win_get_ptr(window,ierr) result(ptr)                 !InTf!
!===============================================================================
! get a one sided communication window (see RPN_COMM_i_win_create) data pointer
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : C compatible (type(C_PTR)) pointer to the data array associated with window
!                  in case of error, C_NULL_PTR is returned (null pointer)
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR) :: ptr                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  ptr = C_NULL_PTR
  if(.not. win_valid(window) ) return

  indx = window%t2            ! window table entry for this window
  ptr = win_tab(indx)%base    ! get data pointer from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_get_ptr                                      !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_size get data pointer associated to a one sided communication window
! SYNOPSIS
function RPN_COMM_i_win_get_size(window,ierr) result(siz)                 !InTf!
!===============================================================================
! get a one sided communication window (see RPN_COMM_i_win_create) data pointer
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : C compatible (type(C_PTR)) pointer to the data array associated with window
!                  in case of error, C_NULL_PTR is returned (null pointer)
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  integer :: siz                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  siz = -1
  if(.not. win_valid(window) ) return

  indx = window%t2            ! window table entry for this window
  siz = win_tab(indx)%siz    ! get data pointer from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_get_size                                      !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_put_r write into a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_put_r(window,larray,targetpe,offset,nelem,ierr) !InTf!
!===============================================================================
! one sided communication remote put (write) into one sided communication window
! from a local array
! it is an error to attempt a "remote" put when the window is not "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (source of put)
! target (IN)     ordinal in window communicator of remote PE
! offset (IN)     displacement (origin 0) into remote PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: targetpe                                     !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local
  type(rpncomm_windef), pointer :: win_entry
  integer(kind=MPI_ADDRESS_KIND) :: offset_8

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_pointer( larray , local, [nelem* win_entry%ext] )   ! pointer to local array

  if(.not. win_entry%active_mode) then
    call mpi_win_lock(MPI_LOCK_EXCLUSIVE, targetpe, 0, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to lock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: locked window',win_entry%win,' on PE',targetpe
    endif
  endif
  offset_8 = offset
  if(debug_mode) print *,'DEBUG: PUT to',targetpe,' at',offset_8+1,' :',local(1:nelem)
  call MPI_put(local,       nelem,        win_entry%typ,   targetpe,   offset_8,    nelem,        win_entry%typ,   win_entry%win,ierr2)
!              ORIGIN_ADDR, ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK,TARGET_DISP, TARGET_COUNT, TARGET_DATATYPE, WIN,          IERROR)
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: remote PUT failed in window',win_entry%win
    return
  endif
  if(.not. win_entry%active_mode) then
    call mpi_win_unlock(targetpe, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to unlock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: unlocked window',win_entry%win,' on PE',targetpe
    endif
  endif

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_put_r                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_acc_r accumulate into a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_acc_r(window,larray,targetpe,offset,nelem,oper,ierr) !InTf!
!===============================================================================
! one sided communication remote accumulate into one sided communication window
! from a local array
! it is an error to attempt a "remote" accumulate when the window is not "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (source of put)
! target (IN)     ordinal in window communicator of remote PE
! offset (IN)     displacement (origin 0) into remote PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! oper   (IN)     rpn comm operator (see RPN_COMM_types.inc and RPN_COMM_i_oper)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
!!  import :: rpncomm_operator                                        !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(rpncomm_operator), intent(IN) :: oper                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: targetpe                                     !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local
  type(rpncomm_windef), pointer :: win_entry
  integer(kind=MPI_ADDRESS_KIND) :: offset_8

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_pointer( larray , local, [nelem* win_entry%ext] )   ! pointer to local array

  if(.not. win_entry%active_mode) then
    call mpi_win_lock(MPI_LOCK_EXCLUSIVE, targetpe, 0, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to lock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: locked window',win_entry%win,' on PE',targetpe
    endif
  endif
  offset_8 = offset
  if(debug_mode) print *,'DEBUG: ACC to',targetpe,' at',offset_8+1,oper%t2,' :',local(1:nelem)
  if(win_entry%opr == MPI_REPLACE) then                  ! replace operation, use put instead of accumulate
    call MPI_put       (local,       nelem,        win_entry%typ,   targetpe,   offset_8,    nelem,        win_entry%typ,           win_entry%win, ierr2)
  else
    call MPI_accumulate(local,       nelem,        win_entry%typ,   targetpe,   offset_8,    nelem,        win_entry%typ, oper%t2,  win_entry%win, ierr2)
!              ORIGIN_ADDR, ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK,TARGET_DISP, TARGET_COUNT, TARGET_DATATYPE, WIN,          IERROR)
  endif
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: remote ACC failed in window',win_entry%win
    return
  endif
  if(.not. win_entry%active_mode) then
    call mpi_win_unlock(targetpe, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to unlock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: unlocked window',win_entry%win,' on PE',targetpe
    endif
  endif

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_acc_r                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_put_l write into a local one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_put_l(window,larray,offset,nelem,ierr)      !InTf!
!===============================================================================
! one sided communication local put (write) into one sided communication window
! from a local array
! it is an error to attempt a "local" put when the window is "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (source of put)
! offset (IN)     displacement (origin 0) into this PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open
  integer :: ierr2, i, indx, extent
  type(rpncomm_windef), pointer :: win_entry
  integer, dimension(:), pointer :: local
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  extent = win_entry%ext
  call c_f_pointer( larray , local, [nelem*extent] )                   ! pointer to local array
  do i = 1, nelem*extent
    win_entry%remote(i+offset*extent) = local(i)
  enddo

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_put_l                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_get_r read a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_get_r(window,larray,target,offset,nelem,ierr) !InTf!
!===============================================================================
! one sided communication remote get (read) from one sided communication window
! into a local array
! it is an error to attempt a "remote" get when the window is not "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (destination of get)
! target (IN)     ordinal in window communicator of remote PE
! offset (IN)     displacement (origin 0) into remote PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: target                                       !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local
  type(rpncomm_windef), pointer :: win_entry
  integer(kind=MPI_ADDRESS_KIND) offset_8

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_pointer( larray , local, [nelem*win_entry%ext] )   ! pointer to local array

  if(.not. win_entry%active_mode) then
    call mpi_win_lock(MPI_LOCK_EXCLUSIVE, target, 0, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif

  offset_8 = offset
  call MPI_get(local,nelem,win_entry%typ,target,offset_8,nelem,win_entry%typ,win_entry%win,ierr2) 
  if(ierr2 .ne. MPI_SUCCESS) return

  if(.not. win_entry%active_mode) then
    call mpi_win_unlock(target, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_get_r                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_l read a local one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_get_l(window,larray,offset,nelem,ierr)      !InTf!
!===============================================================================
! one sided communication local get (read) from one sided communication window
! into a local array
! it is an error to attempt a "local" get when the window is "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (destination of get)
! offset (IN)     displacement (origin 0) into this PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open
  integer :: ierr2, i, indx, extent
  type(rpncomm_windef), pointer :: win_entry
  integer, dimension(:), pointer :: local
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  extent = win_entry%ext
  call c_f_pointer( larray , local, [nelem*extent] )                   ! pointer to local array
  do i = 1, nelem*extent
     local(i) = win_entry%remote(i+offset*extent)
  enddo

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_get_l                                   !InTf!
