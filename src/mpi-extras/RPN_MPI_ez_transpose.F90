!/! RPN_MPI - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
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
!****P* rpn_mpi/Transpose
! DESCRIPTION
! simplified/restricted version of MPI transpose package
!
! proposed layout for GEM solver with FFTW (array shapes)
!
!
!                         ===  with xz and xy transposes  ===
!
!                         dimension(lnix,      lnjy,      gnk)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                 FWD transpose X                     inverse transpose X
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lnix,  lnjy,  lnkx , npex)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
! forward FFT(x) <------> dimension(gni,      lnjy,      lnkx)  <------> inverse FFT(x)
! normalisation           |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lniy,  lnjy,  lnkx,  npey)  [ all i indices on same PE ]
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                 FWD transpose Y                     inverse transpose Y
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lniy,  lnjy,  lnkx,  npey)  [ all j indices on same PE ]
!                                            ^
!                                            |
!                                            v
!                                 tridiagonal solver along j
!                         
!
!                 ===  with xz transpose and ripple solver along y ===
!
!                         dimension(lnix,      lnjy,      gnk)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                 FWD transpose X                     inverse transpose X
!                         |                                  ^
!                         |                                  |
!                         v                                  |
!                         dimension(lnix,  lnjy,  lnkx , npex)
!                         |                                  ^
!                         |                                  |
!                         v                                  |
! forward FFT(x) <------> dimension(gni,  lnjy,  nnuma,  lnkx)  <------> inverse FFT(x)
!                       [ SHARED MEMORY on NUMA node] (nnuma PEs)
!                                            ^
!                                            |
!                                            v
!                     ripple distributed tridiagonal solver along j
!
! example of usage (code taken from TEST_104.F90 in rpn_mpi package)
!
! #include 'RPN_MPI.hf'
!
!   integer(kind=8), dimension(:,:,:), allocatable   :: z, z2
!   integer(kind=8), dimension(:,:,:,:), allocatable :: zt, zy, zty
!   integer :: row_comm, col_comm        
! ! eventually : type(RPN_MPI_Comm) :: row_comm, col_comm
! 
!   call RPN_MPI_transpose_setup(gnk, lnk, row_comm, col_comm, ier)
!
!   call RPN_MPI_transpose_xz(LoC(z), LoC(zt), .true., lnimax*2, lnjmax, gnk, lnk, row_comm, ier)
!
!   call RPN_MPI_ez_transpose_xz(LoC(z), LoC(zt), .true., lnimax*2, lnjmax, lnk, ier)
!
!   call RPN_MPI_transpose_xy(LoC(zy), LoC(zty), .false., lnimaxy*2, lnjmax, lnk, col_comm, ier)
!
!   call RPN_MPI_ez_transpose_xz(LoC(z2), LoC(zt), .false., lnimax*2, lnjmax, lnk, ier)
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
!******
module RPN_MPI_transpose_mod
  use rpn_mpi_mpif
  implicit none
#include <RPN_MPI_mpi_symbols.hf>
  integer, dimension(:), pointer, save :: sendcount => NULL()     ! tables for alltoallv
  integer, dimension(:), pointer, save :: senddispl => NULL()
  integer, dimension(:), pointer, save :: recvcount => NULL()
  integer, dimension(:), pointer, save :: recvdispl => NULL()
  integer, dimension(:), pointer, save :: countnk   => NULL()
  integer, save :: nk   = 0                                       ! number of levels
  integer, save :: npex = 0                                       ! number of PEs in a row
  integer, save :: npey = 0                                       ! number of PEs in a column
  integer, save :: rowcom = MPI_COMM_NULL                         ! row communicator
  integer, save :: colcom = MPI_COMM_NULL                         ! column communicator
  integer, save :: rowrank= 999999
  integer, save :: colrank= 999999
  integer, save :: moxz = 1                                       ! default modulo for staggered alltoall is 1
  integer, save :: moxy = 1                                       ! default modulo for staggered alltoall is 1
  integer, parameter :: MAXTMG = 8
  real(kind=8), dimension(4,MAXTMG), save :: times                ! collect timings for last MAXTMG to fwd/inverse xz and xy transposes
  integer, dimension(4), save :: counts = [0, 0, 0, 0]            ! counters for fwd/inverse xz and xy transposes
 contains
  subroutine RPN_MPI_transpose_alloc(m)
    implicit none
    integer, intent(IN) :: m

    if(associated(countnk))   deallocate(countnk)        ! deallocate if already allocated
    if(associated(sendcount)) deallocate(sendcount)
    if(associated(recvcount)) deallocate(recvcount)
    if(associated(senddispl)) deallocate(senddispl)
    if(associated(recvdispl)) deallocate(recvdispl)
    allocate(countnk(m), sendcount(m), recvcount(m), senddispl(m), recvdispl(m))
    return
  end subroutine RPN_MPI_transpose_alloc
end module RPN_MPI_transpose_mod

subroutine RPN_MPI_print_transpose_times ! print collected transpose timings 
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
  write(0,1) 'fwd xz :',times(1,:)
  write(0,1) 'inv xz :',times(2,:)
  write(0,1) 'fwd xy :',times(3,:)
  write(0,1) 'inv xy :',times(4,:)
1 format(A,8F10.6)
end subroutine RPN_MPI_print_transpose_times

!****f* rpn_mpi/RPN_MPI_transpose_setup
! DESCRIPTION
! setup routine for EZ transposes
!
! setup routine that needs to be called before subsequent calls to
! RPN_MPI_ez_transpose_xz and RPN_MPI_ez_transpose_xy
!
! gnk       : total number of levels
! lnkx      : local number of levels in transposed array zt (the SUM of lnkx MUST be equal to gnk)
! row_comm  : row communicator for transpose_xz
! col_comm  : column communicator for transpose_xy
! ierr      : error flag, same as MPI
!
! SEE ALSO
!  Transpose
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! SYNOPSIS
 subroutine RPN_MPI_transpose_setup(gnk, lnkx, row_comm, col_comm, ierr) !InTf!
! IGNORE
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
! ARGUMENTS
  integer, intent(IN) :: row_comm, col_comm                              !InTf!
  integer, intent(IN) :: gnk, lnkx                                       !InTf!
  integer, intent(OUT) :: ierr                                           !InTf!
!******

  integer :: m

  nk     = gnk
  rowcom = row_comm
  call MPI_Comm_size(rowcom, m, ierr)                ! get size of row coumunicator
  call MPI_Comm_rank(rowcom, rowrank, ierr)          ! rank in row
  if(npex < m) call RPN_MPI_transpose_alloc(m)       ! allocate arrays if not large enough
  npex   = m
  call MPI_Allgather(lnkx, 1, MPI_INTEGER, countnk, 1, MPI_INTEGER, rowcom, ierr) ! get lnk counts
  if(sum(countnk(1:npex)) .ne. gnk .or. &
     any(countnk(1:npex) < 0)) then                  ! ERROR, missing or extra levels
    ierr = MPI_ERROR
    write(0,*) 'ERROR: the sum of lnk(:) does not match gnk'
    return
  endif

  colcom = col_comm
  call MPI_Comm_size(colcom, npey, ierr)             ! get size of column coumunicator
  call MPI_Comm_rank(colcom, colrank, ierr)          ! rank in column

  return
 end subroutine RPN_MPI_transpose_setup                                  !InTf!

!****f* rpn_mpi/RPN_MPI_ez_transpose_xz
! DESCRIPTION
! ezversion of xz transpose (row communicator is implicit)
! see RPN_MPI_transpose_xz for detailed description of arguments
!
! SEE ALSO
!  RPN_MPI_transpose_xz Transpose
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! SYNOPSIS
 subroutine RPN_MPI_ez_transpose_xz(z, zt, forward, lnix, lnjy, lnkx, ierr) !InTf!
! IGNORE
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
! ARGUMENTS
!! ! RPN_MPI_Loc is essentially the wrapped address of some array
  type(RPN_MPI_Loc), intent(IN), value :: z, zt                 !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lnix, lnjy, lnkx                       !InTf!
  integer, intent(OUT) :: ierr                                  !InTf!
! IGNORE
! little white lie in interface, z, zt are advertised as addresses passed by value
!   integer, dimension(lnix,lnjy,*), intent(INOUT)      :: z      ! NO HALO in arrays 
!   integer, dimension(lnix,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npex
!******
  interface
    subroutine RPN_MPI_transpose_xz(z, zt, forward, lnix, lnjy, gnk, lnkx, row_comm, ierr)
      import :: RPN_MPI_Loc
      type(RPN_MPI_Loc), intent(IN), value :: z, zt
      logical, intent(IN) :: forward
      integer, intent(IN) :: lnix, lnjy, lnkx, gnk
      integer, intent(IN) :: row_comm
      integer, intent(OUT) :: ierr 
    end subroutine RPN_MPI_transpose_xz
    subroutine RPN_MPI_transpose_xz_s(z, zt, forward, lnix, lnjy, gnk, lnkx, row_comm, ierr)
      import :: RPN_MPI_Loc
      type(RPN_MPI_Loc), intent(IN), value :: z, zt
      logical, intent(IN) :: forward
      integer, intent(IN) :: lnix, lnjy, lnkx, gnk
      integer, intent(IN) :: row_comm
      integer, intent(OUT) :: ierr 
    end subroutine RPN_MPI_transpose_xz_s
  end interface

  ierr = MPI_ERROR
  if(npex == 0 .or. rowcom == MPI_COMM_NULL .or. nk ==0) return ! not initialized properly

  call RPN_MPI_transpose_xz(z, zt, forward, lnix, lnjy, nk, lnkx, rowcom, ierr)

 end subroutine RPN_MPI_ez_transpose_xz                         !InTf!

!****f* rpn_mpi/RPN_MPI_transpose_xz
! DESCRIPTION
! xz transpose
!
! transpose z along x
!   the z dimension gets distributed over the row, 
!   the x dimension slices get gathered along the row
!
! in the original array a given PE has 
!   all data( dimension(,,gnk) ) along z (all k indices, gnk values)
!   one slice along y (length lnjy)
!   one slice along x (length lnix)
!
! in the transposed array, a given PE has
!   one slice along z (length lnkx)
!   one slice along y (length lnjy)
!   all slices along x (all i indices, gni values, npex slices)
!
! forward transpose (forward == .true.)
! z         : original array, dimension(lnix,lnjy,gnk)
! zt        : transposed arrray, dimension(lnix,lnjy,lnkx,npex)
! lnix      : number of local points along x (assumed to be IDENTICAL on ALL PEs in row)
! lnjy      : number of local points along y (assumed to be IDENTICAL on ALL PEs in row)
! gnk       : total number of levels
! lnkx      : local number of levels in transposed array zt (the SUM of lnkx MUST be equal to gnk)
! row_comm  : row communicator for transpose
! ierr      : error flag, same as MPI
!
! SEE ALSO
!  Transpose
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! SYNOPSIS
 subroutine RPN_MPI_transpose_xz(z0, zt0, forward, lnix, lnjy, gnk, lnkx, row_comm, ierr) !InTf!
! IGNORE
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
! ARGUMENTS
!! ! RPN_MPI_Loc is the wrapped address of some array
  type(RPN_MPI_Loc), intent(IN), value :: z0, zt0               !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lnix, lnjy, gnk, lnkx                  !InTf!
  integer, intent(IN) :: row_comm                               !InTf!
  integer, intent(OUT) :: ierr                                  !InTf!
!   integer, dimension(lnix,lnjy,gnk), intent(INOUT)    :: z         ! NO HALO in arrays 
!   integer, dimension(lnix,lnjy,lnkx,npex), intent(INOUT) :: zt     ! last dimension is npex
!******
  integer :: ix, slot, n, m
  integer, dimension(:), pointer :: z, zt

  call C_F_POINTER(z0%p,  z,  [lnix*lnjy*gnk])
  call C_F_POINTER(zt0%p, zt, [lnix*lnjy*lnkx*npex])
  ierr = MPI_ERROR
  if(npex == 0 .or. row_comm .ne. rowcom .or. nk .ne. gnk) then  ! initialize or reinitialize ?
    nk = gnk
    rowcom = row_comm
    call MPI_Comm_size(rowcom, m, ierr)                ! get size of row coumunicator
    call MPI_Comm_rank(rowcom, rowrank, ierr)          ! rank in row
    if(npex < m) call RPN_MPI_transpose_alloc(m)       ! allocate arrays if not large enough
    npex = m                                           ! new array dimension
    call MPI_Allgather(lnkx, 1, MPI_INTEGER, countnk, 1, MPI_INTEGER, rowcom, ierr) ! get lnk counts
    if(sum(countnk(1:npex)) .ne. gnk) goto 222         ! ERROR, missing or extra levels
  endif

  if(forward) then
    slot = 1
  else
    slot = 2
  endif
  ix = mod(counts(slot), MAXTMG)+1
  counts(slot) = ix

  times(slot,ix) = MPI_Wtime()
  if(colcom .ne. MPI_COMM_NULL) moxz = 1               ! can be set to 2/4/... deactivated for now
  if(forward) then                                     ! z -> zt
    sendcount(1:npex) = countnk(1:npex) * lnix * lnjy  ! may be zero for some PEs
    recvcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    senddispl(1) = 0
    recvdispl(1) = 0
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
    call MPI_Alltoallv(z , sendcount, senddispl, MPI_INTEGER,               &
                       zt, recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
  else                                                 ! zt -> z
    sendcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    recvcount(1:npex) = countnk(1:npex) * lnix * lnjy  ! may be zero for some PEs
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
    call MPI_Alltoallv(zt, sendcount, senddispl, MPI_INTEGER,               &
                       z , recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
  endif
  times(slot,ix) = MPI_Wtime() - times(slot,ix)

1 return
222 continue
  write(0,*) 'ERROR: the sum of lnk() does not match gnk'
  goto 1
 end subroutine RPN_MPI_transpose_xz !InTf!

! version with "rolling wave" alltoall
 subroutine RPN_MPI_transpose_xz_s(z0, zt0, forward, lnix, lnjy, gnk, lnkx, row_comm, ierr) !InTf!
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
!! ! RPN_MPI_Loc is essentially the wrapped address of some array
  type(RPN_MPI_Loc), intent(IN), value :: z0, zt0               !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lnix, lnjy, gnk, lnkx                  !InTf!
  integer, intent(IN) :: row_comm                               !InTf!
  integer, intent(OUT) :: ierr                                  !InTf!
! IGNORE
! little white lie in interface, z, zt are advertised as addresses passed by value
!   integer, dimension(lnix,lnjy,gnk), intent(INOUT)    :: z      ! NO HALO in arrays 
!   integer, dimension(lnix,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npex
!******

  integer :: m, ix, slot, n
  integer, dimension(:), pointer :: z, zt

  call C_F_POINTER(z0%p,  z,  [lnix*lnjy*gnk])
  call C_F_POINTER(zt0%p, zt, [lnix*lnjy*lnkx*npex])
  ierr = MPI_ERROR
  if(npex == 0 .or. row_comm .ne. rowcom .or. nk .ne. gnk) then  ! initialize or reinitialize ?
    nk = gnk
    rowcom = row_comm
    call MPI_Comm_size(rowcom, m, ierr)                ! get size of row coumunicator
    call MPI_Comm_rank(rowcom, rowrank, ierr)          ! rank in row
    if(npex < m) call RPN_MPI_transpose_alloc(m)       ! allocate arrays if not large enough
    npex = m                                           ! new array dimension
    call MPI_Allgather(lnkx, 1, MPI_INTEGER, countnk, 1, MPI_INTEGER, rowcom, ierr) ! get lnk counts
    if(sum(countnk(1:npex)) .ne. gnk) goto 222         ! ERROR, missing or extra levels
  endif

  if(forward) then
    slot = 1
  else
    slot = 2
  endif
  ix = mod(counts(slot), MAXTMG)+1
  counts(slot) = ix

  times(slot,ix) = MPI_Wtime()
  if(colcom .ne. MPI_COMM_NULL) moxz = 1               ! can be set to 2/4/... deactivated for now
  if(forward) then                                     ! z -> zt
    sendcount(1:npex) = countnk(1:npex) * lnix * lnjy  ! may be zero for some PEs
    recvcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    senddispl(1) = 0
    recvdispl(1) = 0
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
    do m = 0 , moxz - 1                                ! perform alltoall by goups of rows to reduce traffic
      if(mod(colrank,moxz) == m) then
        call MPI_Alltoallv(z , sendcount, senddispl, MPI_INTEGER,               &
                           zt, recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
      endif
      if(moxz > 1) call MPI_Barrier(colcom, ierr)
    enddo
  else                                                 ! zt -> z
    sendcount(1:npex) = lnkx * lnix * lnjy             ! may be zero for this PE
    recvcount(1:npex) = countnk(1:npex) * lnix * lnjy  ! may be zero for some PEs
    do n = 2, npex
      senddispl(n) = sendcount(n-1) + senddispl(n-1)
      recvdispl(n) = recvcount(n-1) + recvdispl(n-1)
    enddo
    do m = 0 , moxz - 1
      if(mod(colrank,moxz) == m) then
        call MPI_Alltoallv(zt, sendcount, senddispl, MPI_INTEGER,               &
                           z , recvcount, recvdispl, MPI_INTEGER, rowcom, ierr)
      endif
      if(moxz > 1) call MPI_Barrier(colcom, ierr)
    enddo
  endif
  times(slot,ix) = MPI_Wtime() - times(slot,ix)

1 return
222 continue
  write(0,*) 'ERROR: the sum of lnk() does not match gnk'
  goto 1
 end subroutine RPN_MPI_transpose_xz_s !InTf!

!****f* rpn_mpi/RPN_MPI_ez_transpose_xy
! DESCRIPTION
! ez version of xy transpose (column communicator is implicit)
!
!  see RPN_MPI_transpose_xy for detailed description of arguments
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! SEE ALSO
!  RPN_MPI_transpose_xy Transpose
! SYNOPSIS
 subroutine RPN_MPI_ez_transpose_xy(z, zt, forward, lniy, lnjy, lnkx, ierr) !InTf!
! IGNORE
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
! ARGUMENTS
!! ! RPN_MPI_Loc is essentially the wrapped address of some array
  type(RPN_MPI_Loc), intent(IN), value :: z, zt                 !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lniy, lnjy, lnkx                       !InTf!
  integer, intent(OUT) :: ierr                                  !InTf!
! IGNORE
! little white lie in interface, z, zt are advertised as addresses passed by value
!   integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: z      ! NO HALO in arrays 
!   integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npey
!******
  interface
    subroutine RPN_MPI_transpose_xy(z, zt, forward, lniy, lnjy, lnkx, col_comm, ierr)
      import :: RPN_MPI_Loc
      type(RPN_MPI_Loc), intent(IN), value :: z, zt
      logical, intent(IN) :: forward
      integer, intent(IN) :: lniy, lnjy, lnkx
      integer, intent(IN) :: col_comm
      integer, intent(OUT) :: ierr 
    end subroutine RPN_MPI_transpose_xy
    subroutine RPN_MPI_transpose_xy_s(z, zt, forward, lniy, lnjy, lnkx, col_comm, ierr)
      import :: RPN_MPI_Loc
      type(RPN_MPI_Loc), intent(IN), value :: z, zt
      logical, intent(IN) :: forward
      integer, intent(IN) :: lniy, lnjy, lnkx
      integer, intent(IN) :: col_comm
      integer, intent(OUT) :: ierr 
    end subroutine RPN_MPI_transpose_xy_s
  end interface

  call RPN_MPI_transpose_xy(z, zt, forward, lniy, lnjy, lnkx, colcom, ierr)
 end subroutine RPN_MPI_ez_transpose_xy                         !InTf!

!****f* rpn_mpi/RPN_MPI_transpose_xy
! DESCRIPTION
! xy transpose
!
! transpose x along y, 
!   the original and transposed arrays have the same dimensions
!   the x dimension slices get distributed over the column, 
!   the y dimension slices get gathered along the column
!
! in the original array a given PE has
!   all slices along x (all i indices, gni values, npex slices)
!   one slice along y (length lnjy)
!   one slice along z (length lnkx)
!
! in the transposed array, a given PE has
!   one slice along x (length lniy)
!   all slices along y (all j indices, gnj values, npey slices)
!   one slice along z (length lnkx)
!
! caveat:
!   the output of RPN_MPI_transpose_xz is not a suitable input to RPN_MPI_transpose_xy
!   unless npex == npey (not likely)
!   a reshaping (lni_X,lnj_y,lnk_x,npe_X) -> (lni_Y,lnj_y,lnk_x,npe_Y) is needed
!   the output of RPN_MPI_transpose_xy is not a suitable input to RPN_MPI_transpose_xz
!   a reshaping (lni_Y,lnj_y,lnk_x,npe_Y) -> (lni_X,lnj_y,lnk_x,npe_X) is needed
!
! forward transpose (forward == .true.)
! z         : original array, dimension(lniy,lnjy,lnkx,npey)
! zt        : transposed arrray, dimension(lniy,lnjy,lnkx,npey)
! lniy      : number of local points along x (assumed to be IDENTICAL on ALL PEs in column)
! lnjy      : number of local points along y (assumed to be IDENTICAL on ALL PEs in column)
! lnkx      : local number of levels in transposed array zt (the SUM of lnkx MUST be equal to gnk)
! col_comm  : communicator for this transpose
! ierr      : error flag, same as MPI
!
! SEE ALSO
!  RPN_MPI_ez_transpose_xy Transpose
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! SYNOPSIS
 subroutine RPN_MPI_transpose_xy(z0, zt0, forward, lniy, lnjy, lnkx, col_comm, ierr) !InTf!
! IGNORE
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
! ARGUMENTS
!! ! RPN_MPI_Loc is essentially the wrapped address of some array
  type(RPN_MPI_Loc), intent(IN), value :: z0, zt0               !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lniy, lnjy, lnkx                       !InTf!
  integer, intent(IN) :: col_comm                               !InTf!
  integer, intent(OUT) :: ierr                                  !InTf!
! IGNORE
! little white lie in interface, z, zt are advertised as addresses passed by value
!   integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: z      ! NO HALO in arrays 
!   integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npey
!******

  integer :: ix, slot
  integer, dimension(:), pointer :: z, zt

  call C_F_POINTER(z0%p,  z,  [lniy*lnjy*lnkx])
  call C_F_POINTER(zt0%p, zt, [lniy*lnjy*lnkx])
  ierr = MPI_ERROR
  if(npey == 0 .or. col_comm .ne. colcom) then   ! update/initialize internal tables if necessary
    colcom = col_comm
    call MPI_Comm_size(colcom, npey, ierr)
    call MPI_Comm_rank(colcom, colrank, ierr)
  endif

  if(rowcom .ne. MPI_COMM_NULL) moxy = 1         ! can be set to 2/4/... deactivated for now
  if(forward) then
    slot = 3
  else
    slot = 4
  endif
  ix = mod(counts(slot), MAXTMG)+1
  counts(slot) = ix

  times(slot,ix) = MPI_Wtime()
  if(forward) then                                     ! z -> zt
    call MPI_Alltoall(z , lniy*lnjy*lnkx, MPI_INTEGER,               &
                      zt, lniy*lnjy*lnkx, MPI_INTEGER, colcom, ierr)
  else                                                 ! zt -> z
    call MPI_Alltoall(zt, lniy*lnjy*lnkx, MPI_INTEGER,               &
                      z , lniy*lnjy*lnkx, MPI_INTEGER, colcom, ierr)
  endif
  times(slot,ix) = MPI_Wtime() - times(slot,ix)

 return
 end subroutine RPN_MPI_transpose_xy !InTf!

! version with "rolling wave" alltoall
 subroutine RPN_MPI_transpose_xy_s(z0, zt0, forward, lniy, lnjy, lnkx, col_comm, ierr) !InTf!
  use ISO_C_BINDING
  use RPN_MPI_transpose_mod
  implicit none
!! import :: RPN_MPI_Loc                                        !InTf!
!! ! RPN_MPI_Loc is essentially the wrapped address of some array
  type(RPN_MPI_Loc), intent(IN), value :: z0, zt0               !InTf!
  logical, intent(IN) :: forward                                !InTf!
  integer, intent(IN) :: lniy, lnjy, lnkx                       !InTf!
  integer, intent(IN) :: col_comm                               !InTf!
  integer, intent(OUT) :: ierr                                  !InTf!
! little white lie in interface, z, zt are advertised as addresses passed by value
!   integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: z      ! NO HALO in arrays 
!   integer, dimension(lniy,lnjy,lnkx,*), intent(INOUT) :: zt     ! last dimension is npey

  integer :: ix, slot, m
  integer, dimension(:), pointer :: z, zt

  call C_F_POINTER(z0%p,  z,  [lniy*lnjy*lnkx])
  call C_F_POINTER(zt0%p, zt, [lniy*lnjy*lnkx])
  ierr = MPI_ERROR
  if(npey == 0 .or. col_comm .ne. colcom) then   ! update/initialize internal tables if necessary
    colcom = col_comm
    call MPI_Comm_size(colcom, npey, ierr)
    call MPI_Comm_rank(colcom, colrank, ierr)
  endif

  if(rowcom .ne. MPI_COMM_NULL) moxy = 1         ! can be set to 2/4/... deactivated for now
  if(forward) then
    slot = 3
  else
    slot = 4
  endif
  ix = mod(counts(slot), MAXTMG)+1
  counts(slot) = ix

  times(slot,ix) = MPI_Wtime()
  do m = 0, moxy-1                                ! perform alltoall by goups of columns to reduce traffic
    if(mod(rowrank,moxy) == m) then
      if(forward) then                                     ! z -> zt
        call MPI_Alltoall(z , lniy*lnjy*lnkx, MPI_INTEGER,               &
                          zt, lniy*lnjy*lnkx, MPI_INTEGER, colcom, ierr)
      else                                                 ! zt -> z
        call MPI_Alltoall(zt, lniy*lnjy*lnkx, MPI_INTEGER,               &
                          z , lniy*lnjy*lnkx, MPI_INTEGER, colcom, ierr)
      endif
    endif
    if(moxy > 1) call MPI_Barrier(rowcom, ierr)
  enddo
  times(slot,ix) = MPI_Wtime() - times(slot,ix)

  return
 end subroutine RPN_MPI_transpose_xy_s !InTf!



