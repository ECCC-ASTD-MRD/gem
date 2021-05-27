subroutine rpn_comm_test_007
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM.inc'
!  include 'RPN_COMM_types.inc'
  integer, PARAMETER :: npts = 170000
  integer :: npes, myrank, nprint
  integer, dimension(npts) :: pe, index
  real, dimension(npts,1) :: source
  real, dimension(npts,5) :: source2
  real, dimension(npts) :: rpe
  integer :: i, status, ierr, j, totpts
  type(rpncomm_context) :: context
  real, dimension(:,:), pointer, save :: dest => NULL()
  real, dimension(:,:), pointer, save :: dest2 => NULL()
  real*8 :: T1,T2

  nprint = min(npts,15)
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
  do i=1,npts
    source(i,:)=i
    source2(i,:)=i
  enddo
  call random_number(rpe)
  pe = nint( rpe*(npes-1) )
  pe( 1:npts:170) = -1
  if(myrank==0) then
    T1=MPI_Wtime()
    call shuffle(pe)
    T2=MPI_Wtime()-T1
    print *,'time for shuffle of',npts,' elements =',nint(T2*1000000),' usec'
  endif
!  T1=MPI_Wtime()
!  call I_mrgrnk(pe,index,npts)
!  T2=MPI_Wtime()-T1
!  print *,'time for sorting of',npts,' elements =',nint(T2*1000000),' usec'
!  goto 9
!  pe(5:7) = -1
  if(myrank==0) write(0,1)'PE=',pe(1:npts:max(1,npts/20))
1 format(A,20I7)
  status = RPN_COMM_spread_context(context,'GRID',0,pe,npts)
!  goto 111
  status = RPN_COMM_spread(context,source,npts,size(source,2),dest)
  if(associated(dest)) then
    write(0,1)'dest is now associated and has dimensions',shape(dest)
  endif
  if(status >= 0)  write(0,*)'dest  received',status,' data point(s)'
  if(status <  0)  write(0,*)'ERROR while dest was receiving data'
  call mpi_reduce(status,totpts,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(myrank==0) write(0,1)'total number of received points=',totpts
  write(0,1)'--------------------------------------------'
!  status = RPN_COMM_spread(context,source,npts,size(source,2),dest)
111 continue
  write(0,1)'--------------------------------------------'
!    do j=npts,1,-1
!    write(0,*)'SOURCE data:',nint(source2(j,:))
!    enddo
  status = RPN_COMM_spread(context,source2,npts,size(source2,2),dest2)
  if(associated(dest2)) then
    write(0,1)'dest2 is now associated and has dimensions',shape(dest2)
  endif
  if(status >= 0)  write(0,*)'dest2 received',status,' data point(s)'
  if(status <  0)  write(0,*)'ERROR while dest2 was receiving data'
  call mpi_reduce(status,totpts,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(myrank==0) write(0,1)'total number of received points=',totpts
  
  write(0,*)'--------------------------------------------'
!  status = RPN_COMM_spread(context,source,npts,size(source2,2),dest2)
  write(0,*)'============================================='
9 call mpi_finalize()
  return

  contains

  subroutine Shuffle(a)
    integer, intent(inout) :: a(:)
    integer :: i, randpos, temp
    real :: r

    do i = size(a), 2, -1
      call random_number(r)
      randpos = int(r * i) + 1
      temp = a(randpos)
      a(randpos) = a(i)
      a(i) = temp
    end do
  end subroutine Shuffle

!  integer function RPN_COMM_comm(str)
!  include 'mpif.h'
!  character (len=*) :: str
!  RPN_COMM_comm = MPI_COMM_WORLD
!  return
!  end function RPN_COMM_comm
  
  Subroutine I_mrgrnk (XDONT, IRNGT, NVAL)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In)  :: NVAL
! __________________________________________________________
      Integer :: XVALA, XVALB
!
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
!      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
  End Subroutine I_mrgrnk
end subroutine rpn_comm_test_007
