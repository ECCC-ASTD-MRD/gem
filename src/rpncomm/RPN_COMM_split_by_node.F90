! RPN_COMM - Library of useful routines for C and FORTRAN programming
! Copyright (C) 2019  Division de Recherche en Prevision Numerique
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

module split_by_node
  save
  integer, parameter :: MAX_CACHE=16
  integer, dimension(MAX_CACHE) :: cold      ! original communicator
  integer, dimension(MAX_CACHE) :: cnew      ! communicator for PEs on same node
  integer, dimension(MAX_CACHE) :: commio    ! dommunicator for rank on same node peers
  integer, save :: ncached = 0               ! number of cached entries

end module split_by_node
!****P* rpn_comm/commnicators  (communicator management package)
!******
!InTf!
!****f* rpn_comm/RPN_COMM_split_by_node  split a communicator on a host basis
! SYNOPSIS
subroutine RPN_COMM_split_by_node(origcomm, nodecomm, peercomm, noderank, peerrank, isiz, err)   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! IGNORE
  use ISO_C_BINDING
  use split_by_node
  implicit none
! ARGUMENTS
  integer, intent(IN)  :: origcomm  ! MPI communicator to split on a host basis        !InTf!
  integer, intent(OUT) :: nodecomm  ! new communicator to be used py PEs on same host  !InTf!
  integer, intent(OUT) :: peercomm  ! communicator for node peers                      !InTf!
  integer, intent(OUT) :: noderank  ! rank in new communicator                         !InTf!
  integer, intent(OUT) :: peerrank  ! rank in node peers                               !InTf!
  integer, intent(OUT) :: isiz      ! size of new communicator                         !InTf!
  integer, intent(OUT) :: err       ! error code                                       !InTf!
!******
  include 'mpif.h'
  interface
    function gethostid() result(id) BIND(C,name='gethostid')
      import :: C_LONG
      integer(C_LONG) :: id
    end function gethostid
  end interface
  include 'RPN_COMM_constants.inc'

!   integer :: myhost, myhost0, myhost1, tmpcomm, i, rank
  integer :: i, rank, ierr

  err = RPN_COMM_ERROR      ! precondition for failure
  noderank = -1
  peerrank = -1
  isiz = 0
  nodecomm = MPI_COMM_NULL
  peercomm = MPI_COMM_NULL

  do i = 1, ncached
    if(cold(i) == origcomm) then  ! cached entry found
      nodecomm = cnew(i)
      peercomm = commio(i)
    endif
  enddo

  if(nodecomm == MPI_COMM_NULL) then          ! nothing useful found in cache
    call mpi_comm_rank(origcomm, rank, ierr)
    if(ierr .ne. MPI_SUCCESS) return
    call mpi_comm_split_type(origcomm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, nodecomm, ierr)
    if(ierr .ne. MPI_SUCCESS) return
!     myhost  = gethostid()                    ! host id
!     myhost0 = iand(myhost , Z'7FFFFFFF')     ! lower 31 bits
!     myhost1 = iand( ishft(myhost, -31) , 1 ) ! upper bit
! 
!     call MPI_Comm_split(origcomm , myhost0, rank, tmpcomm, err)    ! split origcomm using the lower 31 bits of host id , weight=rank in origcomm
!     if(err .ne. MPI_SUCCESS) return
!     call MPI_Comm_split(tmpcomm ,myhost1, rank, nodecomm, err)     ! re split using the upper bit of host id , weight=rank in origcomm
!     if(err .ne. MPI_SUCCESS) return
!     call MPI_Comm_free(tmpcomm, err)                               ! tmpcomm no longer needed
    if(ierr .ne. MPI_SUCCESS) return
    call MPI_Comm_rank(nodecomm, noderank, ierr);                   ! rank of this PE on this SMP node
    if(ierr .ne. MPI_SUCCESS) return
    call MPI_Comm_split(origcomm, noderank, rank, peercomm, ierr)   ! split origcomm using noderank as the color to make node peers communicator
    if(ierr .ne. MPI_SUCCESS) return
    call MPI_Comm_rank(peercomm, peerrank, ierr);                   ! rank of this PE in node peers communicator
    if(ierr .ne. MPI_SUCCESS) return

    if(ncached < MAX_CACHE) then                ! add to cache if cache is not full
      ncached = ncached + 1
      cold(ncached)   = origcomm
      cnew(ncached)   = nodecomm
      commio(ncached) = peercomm
    endif
  endif

  call MPI_Comm_rank(nodecomm, noderank,ierr);                       ! rank of this PE on this SMP node (belonging to origcomm)
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_Comm_rank(peercomm, peerrank,ierr);                       ! rank of this PE in the peers communicator
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_Comm_size(nodecomm, isiz, ierr);                          ! number of PEs on this SMP node (belonging to origcomm)
  if(ierr .ne. MPI_SUCCESS) return

  err = RPN_COMM_OK
  return
  
end subroutine RPN_COMM_split_by_node  !InTf!
!****f* rpn_comm/RPN_COMM_split_by_socket  split a communicator on a socket basis
! SYNOPSIS
subroutine RPN_COMM_split_by_socket(origcomm, nodecomm, sockcomm, peercomm, noderank, sockrank, peerrank, isiz, err)   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2020
! IGNORE
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_constants.inc'
! ARGUMENTS
  integer, intent(IN)  :: origcomm  ! MPI communicator to split on a socket basis        !InTf!
  integer, intent(OUT) :: nodecomm  ! new communicator to be used py PEs on same node    !InTf!
  integer, intent(OUT) :: sockcomm  ! new communicator to be used py PEs on same socket  !InTf!
  integer, intent(OUT) :: peercomm  ! new communicator for socket peers                  !InTf!
  integer, intent(OUT) :: noderank  ! rank in node communicator                          !InTf!
  integer, intent(OUT) :: sockrank  ! rank in socket communicator                        !InTf!
  integer, intent(OUT) :: peerrank  ! rank in socket peers                               !InTf!
  integer, intent(OUT) :: isiz      ! size of socket communicator                        !InTf!
  integer, intent(OUT) :: err       ! error code                                         !InTf!
!******
  integer :: socket, ierr, rank
  integer :: lcpus, sockets_per_node, nnuma, ht, map

  call RPN_COMM_split_by_node(origcomm, nodecomm, peercomm, noderank, peerrank, isiz, err) ! split by node
  if(ERR .ne. RPN_COMM_OK) return

  err = RPN_COMM_ERROR      ! precondition for failure
  call get_logical_cpu_configuration(lcpus, sockets_per_node, nnuma, ht, map)
  socket = noderank / (lcpus/sockets_per_node)          ! rank_on_node / cpus_per_socket
!   if(noderank < isiz / 2 ) then
!     socket = 0
!   else
!     socket = 1
!   endif
  call mpi_comm_split(nodecomm, socket, noderank, sockcomm, ierr)  ! re split by socket
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_Comm_rank(sockcomm, sockrank,ierr);                       ! rank in socket communicator
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_Comm_size(sockcomm, isiz, ierr);                          ! number of PEs on this SMP node (belonging to origcomm)
  if(ierr .ne. MPI_SUCCESS) return

  call MPI_Comm_rank(origcomm, rank,ierr);                           ! rank of this PE in the original communicator
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_Comm_split(origcomm, sockrank, rank, peercomm, ierr)      ! create socket peers communicator
  if(ierr .ne. MPI_SUCCESS) return
  call MPI_Comm_rank(peercomm, peerrank,ierr);                       ! rank of this PE in the socket peers communicator
  if(ierr .ne. MPI_SUCCESS) return

  err = RPN_COMM_OK
  return
end subroutine RPN_COMM_split_by_socket  !InTf!

