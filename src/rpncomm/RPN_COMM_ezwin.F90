#define IN_RPN_COMM_ezwin
module RPN_COMM_ezwin_mod
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM_types.inc'
  type :: rpn_comm_ezwin
    type(C_PTR) :: p
    integer     :: ix
    integer     :: win
    integer     :: sz
  end type
  integer, parameter :: RPN_COMM_MAGIC2=int(Z'0AFEFADE')
  integer, parameter :: MAX_EZWINDOWS=64
  type(rpn_comm_ezwin), dimension(MAX_EZWINDOWS), save :: ezwtab
  integer, save :: ent=0
  contains
  function is_invalid(mywin) result(status)
    implicit none
    type(rpncomm_window), intent(IN) :: mywin
    logical :: status

    status = .true.
    if(.not. C_ASSOCIATED(mywin%p)) return
    if( mywin%t1 .ne. ieor(mywin%t2,RPN_COMM_MAGIC2)) return
    status = .false.
  end function is_invalid
end module
subroutine RPN_COMM_ezwin_create(sz,comm,mywin,ierr)       !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window                           !InTf!
  integer, intent(IN) :: sz,comm                           !InTf!
  integer, intent(OUT) :: ierr                             !InTf!
  type(rpncomm_window), intent(OUT) :: mywin               !InTf!
  integer(KIND=MPI_ADDRESS_KIND) :: sz8, address

  mywin%p  = C_NULL_PTR    ! precondition to error condition
  mywin%t1 = 0
  mywin%t2 = 0

  if(ent >= MAX_EZWINDOWS) return    ! table is full
  ent = ent + 1    ! bump table pointer

  sz8 = sz * 4     ! 32 bit integers to bytes
  ! displacement size is 4 bytes (32 bit integers)
  call MPI_WIN_ALLOCATE(sz8, 4, MPI_INFO_NULL, comm, ezwtab(ent)%p,ezwtab(ent)%win,ierr)
  if(ierr .ne. MPI_SUCCESS) then    ! error creating window
    ent = ent - 1                   ! suppress entry in table
    return
  endif
!   print *,"INFO: created window of size",sz8," bytes"

  mywin%p = ezwtab(ent)%p
  ezwtab(ent)%ix = ent
  ezwtab(ent)%sz = sz
  mywin%t1 = ieor(ent,RPN_COMM_MAGIC2)
  mywin%t2 = ent
end subroutine RPN_COMM_ezwin_create                       !InTf!

function RPN_COMM_ezwin_id(mywin) result(winid)            !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window                           !InTf!
  type(rpncomm_window), intent(IN) :: mywin                !InTf!
  integer :: winid                                         !InTf!

  winid = MPI_WIN_NULL
  if(is_invalid(mywin)) return
  winid = ezwtab(mywin%t2)%win
end function RPN_COMM_ezwin_id                             !InTf!

function RPN_COMM_ezwin_ptr(mywin) result(winptr)          !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window,    C_PTR                 !InTf!
  type(rpncomm_window), intent(IN) :: mywin                !InTf!
  type(C_PTR) :: winptr                                    !InTf!

  if(is_invalid(mywin))then
    winptr = C_NULL_PTR
  else
    winptr = ezwtab(mywin%t2)%p
  endif
end function RPN_COMM_ezwin_ptr                            !InTf!

function RPN_COMM_ezwin_size(mywin) result(sz)             !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window                           !InTf!
  type(rpncomm_window), intent(IN) :: mywin                !InTf!
  integer :: sz                                            !InTf!

  sz = -1
  if(is_invalid(mywin)) return
  sz = ezwtab(mywin%t2)%sz
end function RPN_COMM_ezwin_size                           !InTf!

subroutine RPN_COMM_ezwin_fetch_add(mywin, add, fetch, rank, offset, ierr)       !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window                           !InTf!
  type(rpncomm_window), intent(IN) :: mywin                !InTf!
  integer, intent(IN) :: add, rank, offset                 !InTf!
  integer, intent(OUT) :: fetch, ierr                      !InTf!

  integer(KIND=MPI_ADDRESS_KIND) target_disp
  integer :: assert

  ierr = MPI_ERR_OTHER   ! error code in case of premature return
  if(is_invalid(mywin)) return              ! invalid window object
  if(offset >= ezwtab(mywin%t2)%sz) return  ! OUCH, offset points outside of the window
  target_disp = offset
  assert = 0
  call MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, assert, ezwtab(mywin%t2)%win, ierr)
  call MPI_FETCH_AND_OP(add, fetch, MPI_INTEGER, rank, target_disp, MPI_SUM, ezwtab(mywin%t2)%win, ierr)
  call MPI_Win_unlock(rank, ezwtab(mywin%t2)%win, ierr)
end subroutine RPN_COMM_ezwin_fetch_add                    !InTf!

subroutine RPN_COMM_ezwin_get(mywin,pz,nw,rank,offset,ierr)       !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window,    C_PTR                 !InTf!
  type(rpncomm_window), intent(IN) :: mywin                !InTf!
  type(C_PTR), intent(IN), value :: pz                     !InTf!
  integer, intent(IN) :: nw, rank, offset                  !InTf!
  integer, intent(OUT) :: ierr                             !InTf!

  integer(KIND=MPI_ADDRESS_KIND) target_disp
  integer :: assert
  integer, dimension(:), pointer :: z

  ierr = MPI_ERR_OTHER   ! error code in case of premature return
  if(is_invalid(mywin)) return                 ! invalid window object
  if(offset+nw >= ezwtab(mywin%t2)%sz) return  ! OUCH, offset+nw points outside of the window
  target_disp = offset
  assert = 0
  call c_f_pointer(pz,z,[nw])
  call MPI_Win_lock(MPI_LOCK_SHARED, rank, assert, ezwtab(mywin%t2)%win, ierr)
  call MPI_get(z, nw, MPI_INTEGER, rank, target_disp, nw, MPI_INTEGER, ezwtab(mywin%t2)%win, ierr)
  call MPI_Win_unlock(rank, ezwtab(mywin%t2)%win, ierr)
end subroutine RPN_COMM_ezwin_get                          !InTf!

subroutine RPN_COMM_ezwin_put(mywin,pz,nw,rank,offset,ierr)      !InTf!
  use RPN_COMM_ezwin_mod
  implicit none
  include 'mpif.h'
!!  import :: rpncomm_window,    C_PTR                 !InTf!
  type(rpncomm_window), intent(IN) :: mywin                !InTf!
  type(C_PTR), intent(IN), value :: pz                     !InTf!
  integer, intent(IN) :: nw, rank, offset                  !InTf!
  integer, intent(OUT) :: ierr                             !InTf!

  integer(KIND=MPI_ADDRESS_KIND) target_disp
  integer :: assert
  integer, dimension(:), pointer :: z

  ierr = MPI_ERR_OTHER   ! error code in case of premature return
  if(is_invalid(mywin)) return                 ! invalid window object
  if(offset+nw >= ezwtab(mywin%t2)%sz) return  ! OUCH, offset+nw points outside of the window
  target_disp = offset
  assert = 0
  call c_f_pointer(pz,z,[nw])
! print 101,"PUT",z(1:nw)," to",rank,target_disp
! 101 format(A,4I6,A,2I6)
  call MPI_Win_lock(MPI_LOCK_SHARED, rank, assert, ezwtab(mywin%t2)%win, ierr)
  call MPI_put(z, nw, MPI_INTEGER, rank, target_disp, nw, MPI_INTEGER, ezwtab(mywin%t2)%win, ierr)
  call MPI_Win_unlock(rank, ezwtab(mywin%t2)%win, ierr)
end subroutine RPN_COMM_ezwin_put                          !InTf!

#if defined(SELF_TEST)
program test_ezwin
  use ISO_C_BINDING
  include 'RPN_COMM.inc'
  include 'mpif.h'

interface
  include 'new.inc'
end interface

  integer :: comm0, error, szreq, szans, myrank, mysize, i, j, k
  type(rpncomm_window) :: req_win, ans_win
  integer :: me_x, me_y, fetch, partner
  integer, dimension(-1:4,-1:4) :: petab
  type(C_PTR) :: ptrreq, ptrans
  integer, dimension(:), pointer :: freq, fans
  integer, dimension(10), target :: reqbuf, ansbuf
  integer :: ansoffset
  integer, dimension(100) :: expected
  real, dimension(100) :: expectedf
  integer, dimension(2,50) :: req
  real, target :: tmp
  type :: packet
    SEQUENCE
    integer :: n
    integer :: me
    integer :: off
    real :: val
  end type
  type(packet), target :: pak
  type :: area
    SEQUENCE
    integer :: offset
    type(packet), dimension(10) :: pak
  end type
  type(area), pointer :: pakp

  call mpi_init(error)
  comm0 = MPI_COMM_WORLD
  call mpi_comm_rank(comm0,myrank,error)
  call mpi_comm_size(comm0,mysize,error)
  petab = -1
  me_y = myrank/4
  me_x = mod(myrank,4)
  k = 0
  do j = 0,3
  do i = 0,3
    petab(i,j)=k
    k = k + 1
  enddo
  enddo
  if(me_x == 0 .and. me_y == 0)then
    do j = 3,0,-1
      print 100,petab(0:3,j)
    enddo
  endif
100 format(40I4)

  szreq = 1024*1024*3
  szans = 1024*1024
  call RPN_COMM_ezwin_create(szreq,comm0,req_win,error)
  ptrreq = RPN_COMM_ezwin_ptr(req_win)
  call c_f_pointer(ptrreq,freq,[szreq/4])
  freq(1:szreq/4) = 0
  freq(1) = 1
  call c_f_pointer(ptrreq,pakp)

  call RPN_COMM_ezwin_create(szans,comm0,ans_win,error)
  ptrans = RPN_COMM_ezwin_ptr(ans_win)
  call c_f_pointer(ptrans,fans,[szans/4])
  fans(1:szans/4) = 0
  expected = -1

  ansoffset = 0
  do j = -1,1
  do i = -1,1
    if(i .ne.0 .or. j .ne. 0) then
      partner = petab(i+me_x,j+me_y)
      if(partner .ne. -1) then
        call RPN_COMM_ezwin_fetch_add(req_win, 4, fetch, partner, 0, error)
!         print *,"request to",partner," offset =",fetch
        req(1,1+ansoffset) = partner
        req(2,1+ansoffset) = fetch
        reqbuf(1:4) = [1,myrank,ansoffset,100*myrank]
        pak = packet(1,myrank,ansoffset,1.0*myrank)
        expected(1+ansoffset) = reqbuf(4) + partner ! expected answer from partner is element 4 of request + partner rank
        expectedf(1+ansoffset) = pak%val + .01*partner
!         call RPN_COMM_ezwin_put(req_win,c_loc(reqbuf(1)),4,partner,fetch,error)
        call RPN_COMM_ezwin_put(req_win,c_loc(pak),4,partner,fetch,error)
        ansoffset = ansoffset + 1
      endif
    endif
  enddo
  enddo
  call mpi_barrier(comm0,ierr)
  print 103,"target/offset pairs: ",req(:,1:ansoffset)
103 format(A,20(2I3,3X))
  print *,"total requests size =",freq(1)-1
!   print 101,freq(2:freq(1))
  print 105,pakp%pak(1:(freq(1)-1)/4)
101 format(10(4I5,3x))
105 format(10(3I5,F6.2,3x))

  do i = 2,freq(1)-1,4  ! answer to requests
    ansbuf(1) = myrank + freq(i+3)    ! answer = me + element 4 of request
    tmp = myrank*.01 + transfer(freq(i+3),tmp)
!     call RPN_COMM_ezwin_put(ans_win,c_loc(ansbuf(1)),1,freq(i+1),freq(i+2),error)
    call RPN_COMM_ezwin_put(ans_win,c_loc(tmp),1,freq(i+1),freq(i+2),error)
  enddo
  call mpi_barrier(comm0,ierr)
!   print 102,"expected",expected(1:ansoffset)
  print 104,"expected",transfer(expectedf(1:ansoffset),tmp,ansoffset)
  print 104,"answers ",transfer(fans(1:ansoffset),tmp,ansoffset)
102 format(A,40I6)
104 format(A,40F6.2)
  call mpi_finalize(error)
  stop
end program
#endif
