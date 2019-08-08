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
#define IN_RPN_COMM_spread
!InTf!
function RPN_COMM_spread(contxt,source,npts,ndata,dest)  result(status)   !InTf!
  use ISO_C_BINDING                                                       !InTf!
!!  import :: rpncomm_context    !InTf!
  implicit none                                                           !InTf!
#include "RPN_COMM_int.inc"
  include 'RPN_COMM_spread.inc'
  include 'mpif.h'

  type(rpncomm_context), intent(IN) :: contxt          ! object obtained from RPN_COMM_spread_context                !InTf!
  integer, intent(IN) :: npts, ndata                   ! dimensions of source array (used only on root PE)           !InTf!
  real, dimension(npts,ndata), intent(IN) :: source    ! source array (used only on root PE)                         !InTf!
  real, dimension(:,:), pointer, intent(INOUT) :: dest ! pointer to output data                                      !InTf!
  integer :: status                                    ! < 0 : error,  >=0 : number of valid npoints in dest array   !InTf!
!
! NOTES:
!   contxt  is a typed object obtained from function RPN_COMM_spread_context
!           that describes how the data will be spread from the "root" PE (who has the data) and the "client" PEs
!           who will receive some (or none) of that data
!           the application code is responsible for preserving this pointer (save attribute strongly recommended)
!   source  is the data to be "spread". there are npts samples (of ndata values each) to be spread between the PEs
!   npts    number of data samples (must be consistent with metadata from "contxt")
!   ndata   number of values in each sample (must be consistent with dimensions of dest)
!   dest    is a Fortran pointer to a real 2D array supplied by the user. If said pointer is "NULL", it will be allocated
!           with the proper dimensions by RPN_COMM_spread. If said pointer is already associated upon entry, its dimensions
!           will be checked to make sure that the first dimension is ndata and the secon dimension corresponds to the
!           number of samples that this node will receive (information found in the metadata pointed to by "contxt")
!   status  contains the number of valid data samples received. a negative value indocates an error. some PEs may receive no
!           data from the "spread" operation in which case a status value of zero will be returned

  type(spread_context), pointer :: context_entry
  integer, dimension(:), pointer :: offset, order, nlocal
  real, dimension(:,:), allocatable :: source2
  integer :: i, j, ierr, npoints, npe, shiftmsk
  logical :: debug

  debug = .false.
  call c_f_pointer( contxt%p, context_entry)   ! convert C pointer passed by user into a Fortran pointer
  shiftmsk = ishft(1,context_entry%shiftcnt) - 1

  npoints = max(1,context_entry%npoints)      ! 0 length might disturb scatterv, using length of one (1)
  if(context_entry%ntotal .ne. npts .and. context_entry%rootpe == context_entry%myrank) then  ! number of data samples is not consistent
    goto 111
  endif
  if( .not. associated(dest) ) then
    if(debug) write(0,1)'allocating dest with dimensions',ndata,npoints
    allocate(dest(ndata,npoints))
  else
    if( .not. all(shape(dest) == (/ndata,npoints/) ) ) then  ! check that the dimensions of dest are consistent with ndata and "contxt" metadata
      write(0,*)'SIZE ERROR: (RPN_COMM_spread) (dest) expected',shape(dest),' got',(/ndata,npoints/)
      goto 111
    endif
  endif
  call c_f_pointer( context_entry%offset, offset, (/context_entry%n_pes/)   )   ! restore Fortran pointer to offset table
  call c_f_pointer( context_entry%nlocal, nlocal, (/context_entry%n_pes/)   )   ! restore Fortran pointer to length table
  call c_f_pointer( context_entry%order , order , (/context_entry%ntotal/) )    ! restore Fortran pointer to sort index table
  if(debug) then
    if(context_entry%rootpe == context_entry%myrank) npe=size(nlocal)
    if(context_entry%rootpe == context_entry%myrank) write(0,1)'order=',iand(order(1:npoints:npoints/5),shiftmsk)
    if(context_entry%rootpe == context_entry%myrank) write(0,1)'npe   =',npe
    if(context_entry%rootpe == context_entry%myrank) write(0,1)'nlocal=',nlocal(1:npe:max(1,npe/5))
    if(context_entry%rootpe == context_entry%myrank) write(0,1)'offset=',offset(1:npe:max(1,npe/5))
    if(context_entry%rootpe == context_entry%myrank) write(0,1)'shift =',context_entry%shiftcnt
    if(context_entry%rootpe == context_entry%myrank) write(0,1)'context_entry%ntotal=',context_entry%ntotal
    write(0,1)'context_entry%npoints=',context_entry%npoints
  endif

  if(context_entry%rootpe == context_entry%myrank) then  ! root PE has some reordering to do before sending data
    allocate(source2(ndata,npts))
    if(debug) write(0,1)'source2 shape:',shape(source2)
    do i = offset(1)+1, npts   ! offset(1) contains the number of unused points
    do j = 1, ndata            ! do not bother copuing  them into temporary array source2
      source2(j,i) = source(iand(order(i),shiftmsk),j)    ! reorder and transpose data
    enddo
    enddo
  else
    allocate(source2(1,1))  ! dummy allocate to ensure that pointer is valid
  endif
  call mpi_scatterv(source2, max(1,nlocal)*ndata, offset*ndata , MPI_REAL,     &   ! "spread" data samples across PEs as per "contxt"
                    dest   , npoints*ndata                     , MPI_REAL,     &
                    context_entry%rootpe, context_entry%comm, ierr)
  if(debug) then
    do j=context_entry%npoints,1,-(context_entry%npoints/5)
      write(0,1)'j, received data:',j,nint(dest(1:ndata,j))
    enddo
  endif
  deallocate(source2)             ! deallocate temporary array
  status = context_entry%npoints  ! number of valid data samples (>=0)
  return

111 continue
  status = -1

1 format(A,20I7)
  return

end function RPN_COMM_spread                                                                         !InTf!
!InTf!
function RPN_COMM_spread_context(contxt,com,rootpe,pe,npts) result(status)                          !InTf!
  use ISO_C_BINDING                                                                                 !InTf!
!!  import :: rpncomm_context                                                                       !InTf!
  implicit none                                                                                     !InTf!
#include "RPN_COMM_int.inc"
  include 'RPN_COMM_spread.inc'
  include 'mpif.h'

  type(rpncomm_context), intent(OUT) :: contxt     ! C pointer to metadata describing "spread" operation         !InTf!
  character (len=*), intent(IN) :: com             ! RPN_COMM communicator                                       !InTf!
  integer, intent(IN) :: npts                      ! number of data points                                       !InTf!
  integer, intent(IN) :: rootpe                    ! root PE for the spread operation                            !InTf!
  integer, dimension(npts), intent(IN) :: pe       ! destination table, data point i will be sent to PE pe(i)    !InTf!
  integer :: status                                ! 0 if successful, non zero otherwise                         !InTf!
!
! NOTES:
!   contxt     C pointer (see ISO_C_BINBING documentation) to metadata describing the "spread" operation to be performed
!              the application code is expected to keep this information (save attribute strongly recommended)
!              this C pointer is expected to be passed to subsequent RPN_COMM_spread operations
!              each spread pattern (pe,rootpe,com,npts) needs its own "contxt"
!   com        RPN_COMM style communicator for the "spread" operation (e.g. "GRID")
!   rootpe     ordinal of PE that "spreads" the data to the client PEs in communicator com
!   npts       number of data samples to be "spread"
!   pe         pe(i) contains the ordinal in communicator com of the PE that is to receive sample(i)
!              pe(i) == -1 indicates that the sample will not be received by any PE (presumably out of the grid)

  type(spread_context), pointer :: context_entry   ! metadata
  integer, dimension(:), pointer :: offset, order, nlocal
!  integer, external :: RPN_COMM_comm
  integer :: i, max_pe, n_pes, myrank, ierr, comm, shiftcnt

  comm = RPN_COMM_comm(com)               ! get MPI communicator
  call mpi_comm_size( comm, n_pes ,ierr )  ! get number of PEs in communicator
  max_pe = n_pes -1
  call mpi_comm_rank( comm, myrank ,ierr ) ! my rank in communicator

  allocate(context_entry)
  contxt%p = c_loc(context_entry)   ! blind pointer returned to user (to be passed later to RPN_COMM_spread)
  context_entry%comm    = comm                    ! MPI communicator
  context_entry%myrank  = myrank                  ! my rank in MPI communicator
  context_entry%rootpe  = rootpe                  ! root PE rank in MPI communicator

  if(myrank == rootpe) then        ! root PE needs full tables
    context_entry%ntotal  = npts                    ! total number of points
    allocate(order(npts))
    context_entry%n_pes   = n_pes                   ! number of PE s
    allocate(offset(n_pes))
    allocate(nlocal(n_pes))
  else                             ! non root PE allocates minimum dimansion dummy tables
    context_entry%ntotal  = 1                    ! fake total number of points
    allocate(order(1))
    context_entry%n_pes   = 1                   ! fake number of PE s
    allocate(offset(1))
    allocate(nlocal(1))
  endif

  context_entry%npoints = 0                       ! will be adjusted when value is known
  context_entry%order   = c_loc(order(1))         ! pointer to reordering table
  context_entry%offset  = c_loc(offset(1))        ! pointer to offset table
  context_entry%nlocal  = c_loc(nlocal(1))        ! pointer to PE population table
  context_entry%shiftcnt= 0                       ! point number field width
  shiftcnt = 0
  do i = 16, 30                                   ! find if 32 bits can accomodate both point number and PE number
    if(npts < ishft(1,i) .and. max_pe < ishft(1,32-i)) then   ! both fit
      context_entry%shiftcnt= i                   ! shiftcount OK
      shiftcnt = i
      exit
    endif
  enddo
  if(myrank == rootpe) write(0,*)'INFO: (RPN_COMM_spread_context), shift count =',shiftcnt
  if(shiftcnt == 0) go to 3                       ! no way to fit both PE number and point number in 32 bits

  nlocal = 0
  offset(1)=0
  order(1)=99999
  if(myrank == rootpe) then        ! root PE computes full tables, other PEs do nothing
    call I_mrgrnk(pe,order,npts)       ! index sort by target PE number, result in array "order"
    do i = 1, npts
      if(pe(i) >= 0) then              ! point is in grid
        order(i) = ior( ishft(pe(i),shiftcnt) , order(i) ) ! add PE number into index array  (max 512000 points)
        nlocal(pe(i)+1) = nlocal(pe(i)+1) + 1              ! increase count for destination PE i
      else
        offset(1) = offset(1) + 1        ! point not in grid, will sort at b eginning of index, bump first offset
        if(pe(i) > max_pe) then          ! OOPS, target PE number too large
          nlocal = -1                    ! a count of -1 will be sent to all PEs to indicate error
          status = -1
          go to 2
        endif
      endif
    enddo
    do i = 2, n_pes
      offset(i) = offset(i-1) + nlocal(i-1)   ! build offset table for mpi_scatterv from nlocal table
    enddo
    offset = min(npts-1,offset)  ! in case send count for last pes is 0, let's not overindex 
  endif

2 call mpi_scatter(nlocal,1,MPI_INTEGER,context_entry%npoints,1,MPI_INTEGER,rootpe,comm,ierr) ! send local count to each PE
  if(ierr .ne. MPI_SUCCESS) status = -2
  if(context_entry%npoints == -1 .or. ierr .ne. MPI_SUCCESS) go to 3    ! OOPS, target PE number was too large or MPI error

  status = 0                   ! all is well, return
  return

3 continue                     ! upon error, deallocate everything, set error status
  deallocate(order)
  deallocate(offset)
  deallocate(nlocal)
  deallocate(context_entry)
  contxt%p = C_NULL_PTR           ! and return a NULL contxt pointer
  status = -1
  return
  
1 format(A,20I5)
  contains

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

end function RPN_COMM_spread_context                                                    !InTf!
