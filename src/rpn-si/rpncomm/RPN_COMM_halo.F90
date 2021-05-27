module rpn_comm_halos
implicit none
include 'mpif.h'
save

integer, parameter :: HALO_NOT_FOUND=0
integer, parameter :: HALO_ERROR=-1

private :: init_halotab
public  :: find_halo_entry
private :: create_halo_entry

logical, private :: initialized=.false.
integer, private :: next_halo=0

type :: halo
  sequence
  integer :: maxni, maxnj, maxnk
  integer :: halox, haloy
  integer :: comm
  integer, dimension(4) :: ew_request
  integer, dimension(4) :: ns_request
  integer, dimension(:), pointer :: to_w, from_w, to_e, from_e, to_n, from_n, to_s, from_s
end type

integer, private, parameter :: MAX_HALOS=256
type(halo), private, dimension(:), pointer :: halotab

contains
  subroutine init_halotab
  implicit none
  integer :: i
  allocate(halotab(MAX_HALOS))
  do i=1,MAX_HALOS
    halotab(i)%comm = MPI_COMM_NULL
  enddo
  initialized = .true.
  next_halo = 1
  end subroutine init_halotab
  
  function find_halo_entry(ni,nj,nk,halox,haloy,comm) result(index)
  implicit none
  integer :: index
  integer, intent(IN) :: ni,nj,nk,halox,haloy,comm
  integer :: i

  index = HALO_NOT_FOUND
  do i=1,next_halo-1
    if(comm /= halotab(i)%comm) cycle
    if(halox /= halotab(i)%halox) cycle
    if(haloy /= halotab(i)%haloy) cycle
    if(ni*nk > halotab(i)%maxni*halotab(i)%maxnk) cycle
    if(nj*nk > halotab(i)%maxnj*halotab(i)%maxnk) cycle
    if(ni*nk < 0.9*halotab(i)%maxni*halotab(i)%maxnk) cycle
    if(nj*nk < 0.9*halotab(i)%maxnj*halotab(i)%maxnk) cycle
    index = i
    return
  enddo
  end function find_halo_entry
  
  function create_halo_entry(ni,nj,nk,halox,haloy,comm,east,west,north,south) result(index)
  implicit none
  integer :: index
  integer, intent(IN) :: ni,nj,nk,halox,haloy,comm
  integer, intent(IN) :: east,west,north,south    ! neighbors

  index = HALO_ERROR
  if(next_halo >= MAX_HALOS) return
! save halo description
  halotab(next_halo)%halox = halox
  halotab(next_halo)%haloy = haloy
  halotab(next_halo)%maxni = ni
  halotab(next_halo)%maxnj = nj
  halotab(next_halo)%maxnk = nk
  halotab(next_halo)%comm = comm
! allocate send-receive buffers
  allocate(halotab(next_halo)%to_w(nj*nk*halox))
  allocate(halotab(next_halo)%to_e(nj*nk*halox))
  allocate(halotab(next_halo)%from_w(nj*nk*halox))
  allocate(halotab(next_halo)%from_e(nj*nk*halox))
  allocate(halotab(next_halo)%to_n((ni+2*halox)*nk*haloy))
  allocate(halotab(next_halo)%to_s((ni+2*halox)*nk*haloy))
  allocate(halotab(next_halo)%from_n((ni+2*halox)*nk*haloy))
  allocate(halotab(next_halo)%from_s((ni+2*halox)*nk*haloy))
! create persistent requests
  index = next_halo
  next_halo = next_halo + 1
  end function create_halo_entry

end module rpn_comm_halos
