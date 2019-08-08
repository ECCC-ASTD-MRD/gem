module rpn_comm_grids
  type :: rpn_comm_dist_grid
    integer :: grid_id
    integer :: gni, gnj
    integer :: mini, maxi, lni
    integer :: minj, maxj, lnj
    integer, dimension(:), pointer :: start_i, count_i
    integer, dimension(:), pointer :: start_j, count_j
  end type

  integer, parameter :: MAX_GRIDS = 64
  type(rpn_comm_dist_grid), dimension(MAX_GRIDS), save :: gt

  integer, save :: used_grids = 0

contains
  function init_grid_slot(ix, nx, ny) result(indx)   ! initialize a (maybe new) slot in grid table
    implicit none
    integer, intent(IN) :: ix, nx, ny
    integer :: indx

    if(ix > MAX_GRIDS .or. ix <= 0) then
      indx = -1
      return   ! error, bad slot number
    endif

    gt(ix)%grid_id = -1
    gt(ix)%gni     = 0
    gt(ix)%gni     = 0
    gt(ix)%lni     = 0
    gt(ix)%mini    = 0
    gt(ix)%maxi    = 0
    gt(ix)%lnj     = 0
    gt(ix)%minj    = 0
    gt(ix)%maxj    = 0

    if( .not. associated(gt(ix)%start_i) ) then
      allocate(gt(ix)%start_i(nx), gt(ix)%count_i(nx))
      allocate(gt(ix)%start_j(ny), gt(ix)%count_j(ny))
    endif

    gt(ix)%start_i = 0
    gt(ix)%count_i = 0

    gt(ix)%start_j = 0
    gt(ix)%count_j = 0

    indx = ix

    return
  end function init_grid_slot

  function init_new_grid(nx, ny) result(indx)   ! initialize a (maybe new) slot in grid table
    implicit none
    integer, intent(IN) :: nx, ny
    integer :: indx

    integer :: i, ix

    if(used_grids == 0) then
      do i = 1 , MAX_GRIDS
        nullify(gt(i)%start_i)
        nullify(gt(i)%count_i)
        nullify(gt(i)%start_j)
        nullify(gt(i)%count_j)
        gt(i)%grid_id = -1
        used_grids = 1
      enddo
    endif

    indx = -1
    do i = 1 , MAX_GRIDS
      if(gt(i)%grid_id == -1) then  ! unused grid
        ix = i
        exit
      endif
    enddo
    if(ix == -1) return  ! no free slot found, table is full
    indx = init_grid_slot(ix, nx, ny)  ! fill slot with zeroes

    return
  end function init_new_grid

  function find_grid(id) result(indx)   ! find index in table gt associated to grid_id id
    implicit none
    include 'RPN_COMM_constants.inc'
    integer, intent(IN) :: id
    integer :: indx

    integer :: i

    indx = -1
    i = ieor(id,RPN_COMM_MAGIC)
    if(i <= 0 .or. i > MAX_GRIDS) return  ! invalid index, not a grid_id
    if(gt(i)%grid_id /= id)       return  ! wrong grid id, index is not coherent
    indx = i

  end function find_grid

end module rpn_comm_grids

function rpn_comm_create_2dgrid(gni,gnj,mini,maxi,minj,maxj) result (grid_id)  !InTf!
!
! create a 2D grid rpn_comm grid descriptor
!
! gni, gnj                : dimensions of global grid           global_grid(gni,gnj)
! mini, maxi, minj, maxj  : storage dimensions of local grid    local_grid(mini:maxi,minj:maxj)
!
! the functions returns a grid identifier if successful, -1 if unsucessful
! the number of PEs along x and y is taken from waht was determined by rpn_comm_init...
!
  use rpn_comm
  use rpn_comm_grids
  implicit none
  integer, intent(IN) :: gni, gnj, mini, maxi, minj, maxj              !InTf!
  integer :: grid_id                                                   !InTf!

  integer :: ix, i, j, lni, lnj

  grid_id = RPN_COMM_ERROR
  ix = init_new_grid(pe_nx,pe_ny)
  if(ix <= 0) return                ! id <= 0 is an error

  lni = (gni + pe_nx -1) / pe_nx    ! max lni
  lnj = (gnj + pe_ny -1) / pe_ny    ! max lnj
  if( 1-mini > maxi-lni) return     ! halo x problem (halo x is assumed to be 1-mini)
  if( 1-minj > maxj-lnj) return     ! halo y problem (halo y is assumed to be 1-minj)

  gt(ix)%gni     = gni
  gt(ix)%gnj     = gnj
  gt(ix)%mini    = mini
  gt(ix)%maxi    = maxi
  gt(ix)%minj    = minj
  gt(ix)%maxj    = maxj

  gt(ix)%count_i = lni
  gt(ix)%count_i(pe_nx) = gni - lni * (pe_nx - 1)
  gt(ix)%start_i(1) = 1                           ! (origin 1)
  do i = 2 , pe_nx
    gt(ix)%start_i(i) = gt(ix)%start_i(i-1) + gt(ix)%count_i(i-1)
  enddo
  gt(ix)%lni     = gt(ix)%count_i(pe_mex + 1)    ! adjust local lni

  gt(ix)%count_j = lnj
  gt(ix)%count_j(pe_ny) = gnj - lnj * (pe_ny - 1)
  gt(ix)%start_j(1) = 1                           ! (origin 1)
  do j = 2 , pe_ny
    gt(ix)%start_j(j) = gt(ix)%start_j(j-1) + gt(ix)%count_j(j-1)
  enddo
  gt(ix)%lnj     = gt(ix)%count_j(pe_mey + 1)    ! adjust local lnj

  gt(ix)%grid_id = ieor(ix , RPN_COMM_MAGIC)
  grid_id = gt(ix)%grid_id

end function rpn_comm_create_2dgrid                                    !InTf!

function rpn_comm_get_2dgrid(grid_id,dim_i,dim_j,gni,gnj,mini,maxi,minj,maxj,starti,counti,startj,countj) result (status)  !InTf!
!
! get the full description of 2D rpn_comm grid with id = grid_id
!
! dim_i : dimension of output arrays starti and counti  (starti(1:pe_nx) and counti(1:pe_nx) will be filled)
! dim_j : dimension of output arrays startj and countj  (startj(1:pe_ny) and countj(1:pe_ny) will be filled)
!
! pe_nx is the number of MPI tiles along x, pe_ny is the number of tiles along y
! dim_i < pe_nx or dim_j < pe_ny is an error
! a non valid grid_id is an error
!
! gni, gnj                : dimensions of global grid           global_grid(gni,gnj)
! mini, maxi, minj, maxj  : storage dimensions of local grid    local_grid(mini:maxi,minj:maxj)
! starti(i) contains the start along x in global space of the Ith MPI tile (along x) (origin 1)
! counti(i) contains the number of useful points along x for the Ith MPI tile (along x)
! startj(j) contains the start along y in global space of the Jth MPI tile (along y) (origin 1)
! countj(j) contains the number of useful points along y for the Jth MPI tile (along y)
!
! the value of the function is RPN_COMM_ERROR or RPN_COMM_OK
!
  use rpn_comm
  use rpn_comm_grids
  implicit none
  integer, intent(IN) :: grid_id, dim_i, dim_j                                     !InTf!
  integer, intent(OUT) :: gni,gnj,mini,maxi,minj,maxj                              !InTf!
  integer, intent(OUT), dimension(dim_i) :: starti,counti                          !InTf!
  integer, intent(OUT), dimension(dim_j) :: startj,countj                          !InTf!
  integer :: status                                                                !InTf!

  integer :: indx

  status = RPN_COMM_ERROR   ! precondition outputs for failure
  gni = -1
  gnj = -1
  indx = find_grid(grid_id)
  if(indx <= 0 .or. indx > MAX_GRIDS) return

  if(pe_nx .ne. size(gt(indx)%start_i)) return
  if(pe_ny .ne. size(gt(indx)%start_j)) return
  if(dim_i < pe_nx) return
  if(dim_j < pe_ny) return

  gni    = gt(indx)%gni
  gnj    = gt(indx)%gnj
  mini   = gt(indx)%mini
  maxi   = gt(indx)%maxi
  minj   = gt(indx)%minj
  maxj   = gt(indx)%maxj
  starti(1:pe_nx) = gt(indx)%start_i
  counti(1:pe_nx) = gt(indx)%count_i
  startj(1:pe_ny) = gt(indx)%start_j
  countj(1:pe_ny) = gt(indx)%count_j

  status = RPN_COMM_OK
end function rpn_comm_get_2dgrid                                                   !InTf!
function rpn_comm_2dgrid_test(nparams,params) result(status)
  use rpn_comm
  implicit none
  integer :: status
  include 'RPN_COMM_interfaces.inc'
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params
  integer :: gni, gnj, grid_id, mini, maxi, minj, maxj
  integer :: mini2, maxi2, minj2, maxj2, gni2, gnj2
  integer :: starti(Pe_nx+2), counti(Pe_nx+2)
  integer :: startj(Pe_ny+2), countj(Pe_ny+2)
  integer, parameter :: lni = 11
  integer, parameter :: lnj = 7
  integer :: i
  integer :: gid(10)

  status = RPN_COMM_ERROR
    gni = pe_nx*lni - 7
    gnj = pe_ny*lnj - 5
    
    mini = -1
    minj = -2
    maxi = lni + 3
    maxj = lnj + 4
    do i = 1, 5
      grid_id = rpn_comm_create_2dgrid(gni,gnj,mini,maxi,minj,maxj)
      if(pe_me == 0) print 100,'grid id =',grid_id,ieor(grid_id,RPN_COMM_MAGIC)
      gni = gni + 1
      gnj = gnj + 1
      gid(i) = grid_id
    enddo
    do i = 6, 10
      gid(i) = i - 5
    enddo
    do i = 10, 1, -1
      status = rpn_comm_get_2dgrid(gid(i),Pe_nx+2,Pe_ny+2,gni2,gnj2,mini2,maxi2,minj2,maxj2,starti,counti,startj,countj)
      if(pe_me == 0) then
        if(status == RPN_COMM_OK) then
          print 101,'i,gni2,gnj2,mini2,maxi2,minj2,maxj2',i,gni2,gnj2,mini2,maxi2,minj2,maxj2
          print 101,'                              starti',starti(1:pe_nx)
          print 101,'                              counti',counti(1:pe_nx)
          print 101,'                              startj',startj(1:pe_ny)
          print 101,'                              countj',countj(1:pe_ny)
          print 100,'SUCCESS: id =',gid(i),ieor(gid(i),RPN_COMM_MAGIC)
        else
          print 100,'ERROR: id, status =',gid(i),status
        endif
        print *,''
      endif
    enddo
100 format(A,1x,Z8.8,i8)
101 format(A,1x,20i4)
  status = RPN_COMM_OK
end function rpn_comm_2dgrid_test
