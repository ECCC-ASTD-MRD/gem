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
#define IN_RPN_COMM_io_dist_coll
function RPN_COMM_io_dist_coll_check(gni,gnj,halo_x,halo_y) result(status)
  use rpn_comm
  implicit none
  logical :: status
  integer, intent(IN) :: gni,gnj,halo_x,halo_y

  integer :: lni, lnj, lasti, lastj

  status = .true.
  lni = (gni+pe_nx-1)/pe_nx
  lasti = gni - (pe_nx-1)*lni
  if(lasti < max(1,halo_x)) then
    print 101,"CFG: ERROR: size of last tile along x, halo_x, gni, lni ",lasti, halo_x, gni, lni
    print 101,"CFG: ERROR: not enough points along x for problem, need gni >= ",max(lni,halo_x)*(pe_nx-1)+max(1,halo_x)
    status = .false.
  endif

  lnj = (gnj+pe_ny-1)/pe_ny
  lastj = gnj - (pe_ny-1)*lnj
  if(lastj < max(1,halo_y)) then
    print 101,"CFG: ERROR: size of last tile along y, halo_y, gnj, lnj ",lastj, halo_y, gnj, lnj
    print 101,"CFG: ERROR: not enough points along y for problem, need gnj >= ",max(lnj,halo_y)*(pe_ny-1)+max(1,halo_y)
    status = .false.
  endif

101 format(A,20I5)
end function RPN_COMM_io_dist_coll_check
subroutine RPN_COMM_io_dist_coll_test(nparams,params)
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces_int.inc>
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params
  logical :: periodx, periody
  integer :: setno, me_io, n_io
  integer :: i, j, k, irep
  integer :: gni=10
  integer :: gnj=7
  integer :: dnk=2
  integer :: lnk=6
  integer :: iope_extra=3
  integer :: halo_x = 1
  integer :: halo_y = 1
  integer :: nrep_test=1
  integer, parameter :: MAX_PRINT=80
  integer :: max_dist_k
  integer, dimension(:), allocatable :: liste_k, liste_k2
  logical, dimension(:), allocatable  :: liste_o, liste_o2
  integer, dimension(:,:,:), allocatable :: global,global2
  integer, dimension(:,:,:), allocatable :: glob_l,glob_l2
  integer, dimension(:,:,:), allocatable :: local, loc_l2
  integer :: lni, lnj
  integer :: mini,maxi,minj,maxj,status
  integer, dimension(pe_nx) :: start_x,count_x
  integer, dimension(pe_ny) :: start_y,count_y
  integer :: i0,in,j0,jn,nerrors,nvalid,expected,effective_lnk,nerrors2,nvalid2
  logical, external :: RPN_COMM_io_dist_coll_check
  integer, external :: RPN_COMM_shuf_ezdist, RPN_COMM_shuf_ezdist_h, RPN_COMM_shuf_ezcoll
  integer :: grid_id
!
  periodx = .false.
  periody = .false.
  lni = 0
  lnj = 0
!  goto 1
  if(pe_me == 0 .and. nparams > 2) then
    print 101,"CFG: nparams",nparams
    print 101,"CFG: params",params
  endif
  if(nparams > 2) then
    gni = params(1)
    gnj = params(2)
    dnk = params(3)
  endif
  if(nparams > 3) lnk = params(4)
  max_dist_k = lnk
  if(nparams >= 5) iope_extra = params(5)
  if(nparams >= 7) then
    halo_x = params(6)
    halo_y = params(7)
  endif
  if(nparams >= 8) max_dist_k = params(8)
  if(nparams >= 9) nrep_test = params(9)
  if(pe_me == 0) then
    print 101,"CFG: gni,gnj,dnk,halo_x,halo_y",gni,gnj,dnk,halo_x,halo_y
    print 101,"CFG: lnk,max_dist_k,pe_nx,pe_ny",lnk,max_dist_k,iope_extra,pe_nx,pe_ny
  endif
1 continue
  if(.not. RPN_COMM_io_dist_coll_check(gni,gnj,halo_x,halo_y)) return
!
  allocate( global(gni,gnj,dnk),  glob_l(1-halo_x:gni+halo_x,1-halo_y:gnj+halo_y,dnk))
  allocate(global2(gni,gnj,dnk))
  allocate(liste_k(dnk))
  liste_k = 0
  allocate(liste_k2(dnk))
  liste_k2 = 0
  allocate(liste_o(lnk), liste_o2(lnk))
  liste_o = .false.
  liste_o(max_dist_k) = .true.     ! test of warning for attempt to redistribute a level
  if(max_dist_k<lnk) liste_o(max_dist_k+1:lnk) = .false.
  liste_o2 = liste_o
!
! along X
!
  if(lni == 0) lni = (gni+pe_nx-1)/pe_nx
  mini = 1-halo_x
  maxi = lni+halo_x
  count_x = lni
  start_x(1) = 1
  do i = 2,pe_nx
    start_x(i) = start_x(i-1) + count_x(i-1)
  enddo
  count_x(pe_nx) = gni + 1 - start_x(pe_nx)
  lni = count_x(pe_mex+1)
  i0 = mini
  if(pe_mex == 0 .and. (.not. periodx)) i0 = 1
  in = maxi
  if(pe_mex == pe_nx-1 .and. (.not. periodx)) in = count_x(pe_nx)
  if(pe_me == 0) then
    print *,"start_x =",start_x
    print *,"count_x =",count_x
  endif
!
! along Y
!
  if(lnj == 0) lnj = (gnj+pe_ny-1)/pe_ny
  minj = 1-halo_y
  maxj = lnj+halo_y
  count_y = lnj
  start_y(1) = 1
  do j = 2,pe_ny
    start_y(j) = start_y(j-1) + count_y(j-1)
  enddo
  count_y(pe_ny) = gnj + 1 - start_y(pe_ny)
  lnj = count_y(pe_mey+1)
  j0 = minj
  if(pe_mey == 0 .and. (.not. periody)) j0 = 1
  jn = maxj
  if(pe_mey == pe_ny-1 .and. (.not. periody)) jn = count_y(pe_ny)
  if(pe_me == 0) then
    print *,"start_y =",start_y
    print *,"count_y =",count_y
  endif
!====================================================================================
  grid_id = rpn_comm_create_2dgrid(gni,gnj,mini,maxi,minj,maxj) ! create 2D grid descriptor
!====================================================================================
  if (grid_id == -1) then
    if(pe_me == 0) print 101,"ERROR: cannot create grid id, gni,gnj,mini,maxi,minj,maxj=",gni,gnj,mini,maxi,minj,maxj
    return
  endif
!
  allocate( local(mini:maxi,minj:maxj,lnk))
  allocate(loc_l2(mini:maxi,minj:maxj,lnk))
  local = 999999
  global = 888888
  glob_l = 777777
! create IO PE set
  if(pe_me == 0) then
    print *,'IO PE number of PEs =',min( min(pe_nx,pe_ny)+iope_extra , lnk)
    print *,'pe_nx,pe_ny,iope_extra,lnk',pe_nx,pe_ny,iope_extra,lnk
  endif
!  setno = RPN_COMM_create_io_set( min( min(pe_nx,pe_ny)+iope_extra , lnk) ,0)  ! make sure not to overflow lnk
  setno = RPN_COMM_create_io_set( min(pe_nx,pe_ny)+iope_extra , 0)  ! may overflow lnk
  if(setno <= 0) then
    print *,'ERROR: IO PE set creation error, quitting',setno
    return
  else
    print *,'INFO: IO PE set created :',setno
  endif
  me_io = RPN_COMM_is_io_pe(setno)     ! me in IO_set
  n_io = RPN_COMM_io_pe_size(setno)    ! IO_set population
!
  do irep = 1,nrep_test ! start of test  repetition loop
    liste_o = .false.
    liste_o(max_dist_k) = .true.     ! test of warning for attempt to redistribute a level
    if(max_dist_k<lnk) liste_o(max_dist_k+1:lnk) = .false.
    liste_o2 = liste_o
    if(me_io .ne. -1) then
      print *,"I am a busy IO pe!",me_io+1,' of',n_io
      do k=1,dnk               ! global levels that this PE will distribute
  !      liste_k(k) = lnk - me_io - (k-1)*n_io   !  levels lnk -> (lnk - nio*dnk + 1)
        liste_k(k) = 1 + me_io + (k-1)*n_io     !  levels 1 -> nio*dnk - 1
        if(liste_k(k) > max_dist_k) liste_k(k) = 0
        liste_k(k) = max( liste_k(k) , 0)
        liste_k2 = liste_k
        do j = 1,gnj
        do i = 1,gni
          global(i,j,k) = liste_k(k) + j*100 + i*10000
        enddo
        enddo
        do j = 1-halo_y,gnj+halo_y
        do i = 1-halo_x,gni+halo_x
          glob_l(i,j,k) = liste_k(k) + (j+100)*100 + (i+100)*100000
        enddo
        enddo
      enddo
! print *,'==== global array ====',liste_k
! do j = gnj,1,-1
!   print 102,' ',global(:,j,1)
! enddo
! print *,'==== glob_l array ====',liste_k
! do j = gnj+halo_y,1-halo_y,-1
!   print 102,' ',glob_l(:,j,1)
! enddo
      print *,"level list =",liste_k
      do k= dnk,1,-1
        if(liste_k(k) <= 0) cycle
        print *,"===== source level ==",liste_k(k),"  ====="
        do j=gnj,1,-1
          if(gni*gnj < MAX_PRINT) print 100,j,global(:,j,k)
        enddo
      enddo
    else
      print *,"I am a relaxed  NON-IO pe !"
    endif
! print *,'lni,lnj,mini,maxi,minj,maxj',lni,lnj,mini,maxi,minj,maxj
!return
#define EZ_TEST
#if defined(EZ_TEST)
!====================================================================================
! use previously obtained grid id
    print *,'=== RPN_COMM_shuf_ezdist ==='
    status = RPN_COMM_shuf_ezdist(setno, grid_id, global, dnk, local, lnk, liste_k, liste_o)
    status = RPN_COMM_shuf_ezdist_h(setno, grid_id, glob_l, dnk, loc_l2, lnk, liste_k, liste_o2)
!====================================================================================
#else
    print *,'=== RPN_COMM_shuf_dist ==='
    call RPN_COMM_shuf_dist(setno,  &
                            global,gni,gnj,dnk,  &
                            local,mini,maxi,minj,maxj,lnk,  &
                            liste_k,liste_o,  &
                            start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                            periodx,periody,status)
    call RPN_COMM_shuf_dist(setno,  &
                            glob_l,gni+2*halo_x,gnj+2*halo_y,dnk,  &
                            loc_l2,mini,maxi,minj,maxj,lnk,  &
                            liste_k2,liste_o2,  &
                            start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                            periodx,periody,status)
#endif
    print *,'liste_o, liste_o2 apres=',liste_o, liste_o2
    if(status .ne. RPN_COMM_OK) then
      print 101,"ERROR: RPN_COMM_shuf_dist failure, lnk,dnk,n_io ",lnk,dnk,n_io
      return
    endif
  enddo ! end of test repetition loop
!
  nerrors  = 0
  nerrors2 = 0
  nvalid   = 0
  nvalid2  = 0
  do k = lnk,1,-1
    if(liste_o(k)) then
      do j = j0,jn
      do i = i0,in
        nvalid = nvalid + 1
        expected = k + (start_y(pe_mey+1)+j-1)*100 + (start_x(pe_mex+1)+i-1)*10000
        if(expected .ne. local(i,j,k)) then
          if(nerrors<0)print *,'i,j,k,expected,local(i,j,k)',i,j,k,expected,local(i,j,k)
          nerrors = nerrors + 1
         if(nerrors>3)goto 666
        endif
      enddo
      enddo
      do j = mini,maxi
      do i = minj,maxj
        nvalid2 = nvalid2 + 1
        expected = k + (start_y(pe_mey+1)+j-1+100)*100 + (start_x(pe_mex+1)+i-1+100)*100000
        if(expected .ne. loc_l2(i,j,k)) then
          if(nerrors2<1)print *,'i,j,k,expected,local(i,j,k)',i,j,k,expected,loc_l2(i,j,k)
          nerrors2 = nerrors2 + 1
!          if(nerrors>3)goto 666
        endif
      enddo
      enddo
    else
      print *,'no data at level',k
    endif
  enddo
  print 101,"nerrors,  nvalid,  npts =",nerrors,nvalid,lnk*(in-i0+1)*(jn-j0+1),(in-i0+1),(jn-j0+1),lnk
  print 101,"nerrors2, nvalid2, npts =",nerrors2,nvalid2,lnk*(maxi-mini+1)*(maxj-minj+1),(maxi-mini+1),(maxj-minj+1),lnk

! return
  
  if(nerrors == 0) goto 777
!  goto 777
666 continue
#if  ! defined(DEPRECATED)
  if(gni*gnj < MAX_PRINT) then
    do k = lnk,1,-1
      if(liste_o(k)) then
        print *,"===== level",k," local ====="
        do j = maxj,minj,-1
  !        print 100,j,local(i0:in,j,k)
          print 100,j,local(:,j,k)
        enddo
      else
        print *,'no data at level',k
      endif
    enddo
  endif
#endif
!
100 format(I3,20I7.6)
101 format(A,20I5)
102 format(A,20I9)
!=======================================================================================
777 continue
!=======================================================================================
!
! collect test follows distribute test
!
  do irep = 1,nrep_test ! start of test repeat loop
    global2 = 999999
    liste_k2 = -999999
    effective_lnk = 1
    do k = 1 , lnk
      if(liste_o(k)) effective_lnk = k  ! highest level available
    enddo
    local = 989898
    do k = 1,effective_lnk
    do j = 1,lnj
    do i = 1,lni
      local(i,j,k) = k + (j - 1 + start_y(pe_mey+1))*100 + (i -1 + start_x(pe_mex+1))*10000
    enddo
    enddo
    enddo
    print *,"====== before  shuf_coll, max lnk ======",effective_lnk
#if defined(EZ_TEST)
!====================================================================================
! use previously obtained grid id
    status = RPN_COMM_shuf_ezcoll(setno, grid_id, global2, dnk, local, effective_lnk, liste_k2)
!====================================================================================
#else
    call RPN_COMM_shuf_coll(setno,  &
                            global2,gni,gnj,dnk,  &
                            local,mini,maxi,minj,maxj,effective_lnk,  &
                            liste_k2,  &
                            start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                            status)
#endif
    if(status .ne. RPN_COMM_OK) then
      print 101,"ERROR: RPN_COMM_shuf_coll failure, lnk,dnk,n_io ",effective_lnk,dnk,n_io
      return
    endif
  enddo ! end of test repeat loop
! global2 should be identical to global once k has been adjusted
! expected k + global - mod(global,10)
  if(me_io .ne. -1) then ! I am an IO PE
    print *,"====== after shuf_coll ======"
    print *,"DEBUG: pe_me, liste_k2", pe_me,liste_k2
    nerrors = 0
    nvalid = 0
    do k = 1,dnk
      if(liste_k2(k) <= 0) cycle
      do j = 1,gnj
      do i = 1,gni
        expected = liste_k2(k) + j*100 + i*10000
        nvalid = nvalid + 1
        if(expected .ne. global2(i,j,k)) nerrors = nerrors + 1
        if(nerrors == 1) then
          print *,"ERROR: expected, got =",expected,global2(i,j,k)
        endif
      enddo
      enddo
    enddo
    print *,"nerrors, npoints =",nerrors,nvalid
    do k = 1, dnk
      if(liste_k2(k) > 0) then
        print *,"===== k, collected level",k,liste_k2(k),"  ====="
        do j = gnj , 1 , -1
          if(gni*gnj < MAX_PRINT) print 100,j,global2(:,j,k)
        enddo
      endif
    enddo
  endif
  return
end subroutine RPN_COMM_io_dist_coll_test
!====================================================================================
! global array dimensions are assumed to be gni x gnj
!
!! integer, external :: RPN_COMM_shuf_ezdist  !InTfX!   ! best interface we can provide for the time being
function RPN_COMM_shuf_ezdist(setno, grid_id, global, dnk, local, lnk, liste_i, liste_o) result (status)
! important notes:
!      it is the caller's responsibility to ensure that liste_o is properly initialized to .false. 
!      before calling the distribute function
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces_int.inc>
  integer, intent(IN) :: setno                     ! IO processor set (from RPN_COMM_create_io_set)
  integer, intent(IN) :: grid_id                   ! grid identifier (from rpn_comm_create_2dgrid)
  integer, intent(IN) :: dnk                       ! number of levels to distribute
  integer, intent(IN) :: lnk                       ! number of levels in the "local" array
  integer, intent(IN), dimension(*)     :: global
  integer, intent(OUT), dimension(*)    :: local
  integer, intent(IN), dimension(dnk)   :: liste_i ! needed only on IO Pes, list of levels to distribute
  logical, intent(INOUT), dimension(lnk):: liste_o ! liste_o(k) will be set to .true. if level k received
  integer :: status                                ! RPN_COMM_OK or RPN_COMM_ERROR
!
  integer, dimension(pe_nx)    :: start_x ! PE (i-1,any) points start at start_x(i) in global space (X direction)
  integer, dimension(pe_nx)    :: count_x ! PE (i-1,any) contains count_x(i) points in the X direction
  integer, dimension(pe_ny)    :: start_y ! PE (any,j-1) points start at start_y(j) in global space (Y direction)
  integer, dimension(pe_ny)    :: count_y ! PE (any,j-1) contains count_y(j) points in the Y direction
  integer :: gni,gnj                      ! horizontal dimensions of the "global" array
  integer :: mini,maxi,minj,maxj          ! horizontal dimensions of the "local" array

  status = rpn_comm_get_2dgrid(grid_id,pe_nx,pe_ny,gni,gnj, mini,maxi,minj,maxj,start_x,count_x,start_y,count_y)
!   print *,'DEBUG: rpn_comm_get_2dgrid',mini,maxi,minj,maxj
!
  if(status .ne. RPN_COMM_OK) then
    status = RPN_COMM_ERROR
    return
  endif

  call RPN_COMM_shuf_dist(setno,  &
                          global,gni,gnj,dnk,  &
                          local,mini,maxi,minj,maxj,lnk,  &
                          liste_i,liste_o,  &
                          start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                          .false.,.false.,status)
end function RPN_COMM_shuf_ezdist
!====================================================================================
! variant with assumed global halo
! global array dimensions are assumed to be  gni + 2*halox x gnj + 2*haloy
!
!! integer, external :: RPN_COMM_shuf_ezdist_h  !InTfX!   ! best interface we can provide for the time being
function RPN_COMM_shuf_ezdist_h(setno, grid_id, global, dnk, local, lnk, liste_i, liste_o) result (status)
! important notes:
!      it is the caller's responsibility to ensure that liste_o is properly initialized to .false. 
!      before calling the distribute function
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces_int.inc>
  integer, intent(IN) :: setno                     ! IO processor set (from RPN_COMM_create_io_set)
  integer, intent(IN) :: grid_id                   ! grid identifier (from rpn_comm_create_2dgrid)
  integer, intent(IN) :: dnk                       ! number of levels to distribute
  integer, intent(IN) :: lnk                       ! number of levels in the "local" array
  integer, intent(IN), dimension(*)     :: global
  integer, intent(OUT), dimension(*)    :: local
  integer, intent(IN), dimension(dnk)   :: liste_i ! needed only on IO Pes, list of levels to distribute
  logical, intent(INOUT), dimension(lnk):: liste_o ! liste_o(k) will be set to .true. if level k received
  integer :: status                                ! RPN_COMM_OK or RPN_COMM_ERROR
!
  integer, dimension(pe_nx)    :: start_x ! PE (i-1,any) points start at start_x(i) in global space (X direction)
  integer, dimension(pe_nx)    :: count_x ! PE (i-1,any) contains count_x(i) points in the X direction
  integer, dimension(pe_ny)    :: start_y ! PE (any,j-1) points start at start_y(j) in global space (Y direction)
  integer, dimension(pe_ny)    :: count_y ! PE (any,j-1) contains count_y(j) points in the Y direction
  integer :: gni,gnj                      ! horizontal dimensions of the "global" array
  integer :: mini,maxi,minj,maxj          ! horizontal dimensions of the "local" array
  integer :: ghx, ghy

  status = rpn_comm_get_2dgrid(grid_id,pe_nx,pe_ny,gni,gnj, mini,maxi,minj,maxj,start_x,count_x,start_y,count_y)
!   print *,'DEBUG: rpn_comm_get_2dgrid',mini,maxi,minj,maxj
!
  if(status .ne. RPN_COMM_OK) then
    status = RPN_COMM_ERROR
    return
  endif

  ghx = 1-mini  ! global halo same as local halo
  ghy = 1-minj
  call RPN_COMM_shuf_dist(setno,  &
                          global, gni+2*ghx, gnj+2*ghy, dnk,  &
                          local,mini,maxi,minj,maxj,lnk,  &
                          liste_i,liste_o,  &
                          start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                          .false.,.false.,status)
end function RPN_COMM_shuf_ezdist_h
!====================================================================================
!! integer, external :: RPN_COMM_shuf_ezcoll  !InTfX!   ! best interface we can provide for the time being
function RPN_COMM_shuf_ezcoll(setno, grid_id, global, dnk, local, lnk, liste_o) result (status)
! important notes:
!      it is the caller's responsibility to ensure that liste_o is properly initialized to negative 
!      numbers before calling the collect function
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces_int.inc>
  integer, intent(IN) :: setno                     ! IO processor set (from RPN_COMM_create_io_set)
  integer, intent(IN) :: grid_id                   ! grid identifier (from rpn_comm_create_2dgrid)
  integer, intent(IN) :: dnk                       ! number of levels to collect
  integer, intent(IN) :: lnk                       ! number of levels in the "local" array
  integer, intent(OUT), dimension(*)   :: global
  integer, intent(IN), dimension(*)    :: local
  integer, intent(OUT), dimension(dnk) :: liste_o  ! is set to k when 2D array level k has been received
  integer :: status                                ! RPN_COMM_OK or RPN_COMM_ERROR
!
  integer, dimension(pe_nx)    :: start_x ! PE (i-1,any) points start at start_x(i) in global space (X direction)
  integer, dimension(pe_nx)    :: count_x ! PE (i-1,any) contains count_x(i) points in the X direction
  integer, dimension(pe_ny)    :: start_y ! PE (any,j-1) points start at start_y(j) in global space (Y direction)
  integer, dimension(pe_ny)    :: count_y ! PE (any,j-1) contains count_y(j) points in the Y direction
  integer :: gni,gnj                      ! horizontal dimensions of the "global" array
  integer :: mini,maxi,minj,maxj          ! horizontal dimensions of the "local" array

  status = rpn_comm_get_2dgrid(grid_id,pe_nx,pe_ny,gni,gnj, mini,maxi,minj,maxj,start_x,count_x,start_y,count_y)
!   print *,'DEBUG: rpn_comm_get_2dgrid',mini,maxi,minj,maxj
!
  if(status .ne. RPN_COMM_OK) then
    status = RPN_COMM_ERROR
    return
  endif

  call RPN_COMM_shuf_coll(setno,  &
                          global,gni,gnj,dnk,  &
                          local,mini,maxi,minj,maxj,lnk,  &
                          liste_o,  &
                          start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                          status)
end function RPN_COMM_shuf_ezcoll
!====================================================================================
!
! distribute dnk 2D arrays, destination is dnk of lnk 2D plane of local 3D array
! using IO PE set setno
!
! status will be RPN_COMM_OK or RPN_COMM_ERROR
!
! liste_i contains the list of 2D planes to be inserted into the 3D destination array
! liste_o(k) is set to .true. when 2D plane k has been received
! start_x(i) contains the first global x index of local array in grid PE (i-1,any)
! count_x(i) contains the number of points along x of local array in grid PE (i-1,any)
! nx is the dimension of start_x and count_x (ERROR if not equal to pe_nx)
! start_y(j) contains the first global y index of local array in grid PE (any,j-1)
! count_y(j) contains the number of points along y of local array in grid (any,j-1)
! ny is the dimension of start_y and count_y (ERROR if not equal to pe_ny)
! mini,maxi,minj,maxj are used to determine the halo width in the local array
!                     this is used to save a halo exchange after distribution
!
! start_x, start_y  are in origin(1)
!
! important notes:
!    the following INPUT parameters MUST be the same on ALL PEs
!      setno                    IO set number, as obtained from RPN_COMM_create_io_set
!      mini,maxi,minj,maxj,lnk  (dimensions of the "local" array)
!      gni,gnj,dnk              (dimensions of the "global" array) (may include a halo)
!      nx,ny,start_x,count_x,start_y,count_y  (start_x, start_y : ORIGIN 1)
!      periodx,periody
!
!    the following OUTPUT parameter WILL be the same on ALL PEs
!      liste_o
!
!    even if not used/necessary on a given PE
!      the "global" array MUST exist ( (1,1,1) array is OK on non IO PEs )
!      array liste_i MUST exist  ( a dimension(1) array is OK on non IO PEs )
!
!    attempting to distribute level k where k > lnk is an ERROR
!
!====================================================================================
subroutine RPN_COMM_shuf_dist(setno,  &
                              global,gni,gnj,dnk,  &
                              local,mini,maxi,minj,maxj,lnk,  &
                              liste_i,liste_o,  &
                              start_x,count_x,nx,start_y,count_y,ny,  &
                              periodx,periody,status)
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno                     ! IO processor set (from RPN_COMM_create_io_set)
  integer, intent(IN) :: gni,gnj,dnk               ! dimensions of the "global" array
  integer, intent(IN) :: mini,maxi,minj,maxj,lnk   ! dimensions of the "local" array
  integer, intent(IN) :: nx,ny                     ! number of PEs along X and Y (MUST be the same as pe_nx and pe_ny from module rpn_comm)
  integer, intent(IN), dimension(gni,gnj,dnk) :: global  !  gni and gnj MAY INCLUDE A HALO
  integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
  integer, intent(IN), dimension(nx)    :: start_x ! PE (i-1,any) points start at start_x(i) in global space (X direction)
  integer, intent(IN), dimension(nx)    :: count_x ! PE (i-1,any) contains count_x(i) points in the X direction
  integer, intent(IN), dimension(ny)    :: start_y ! PE (any,j-1) points start at start_y(j) in global space (Y direction)
  integer, intent(IN), dimension(ny)    :: count_y ! PE (any,j-1) contains count_y(j) points in the Y direction
  integer, intent(IN), dimension(dnk)   :: liste_i ! needed only on IO Pes, list of dnk levels to distribute
  logical, intent(OUT), dimension(lnk)  :: liste_o ! liste_o(k) will be set to .true. if level k received
  logical, intent(IN) :: periodx,periody           ! periodicity along X or Y ( periody MUST be .false. for this version)
  integer, intent(OUT) :: status                   ! RPN_COMM_OK or RPN_COMM_ERROR
  integer :: i, k, low, high, listeik, n
  integer :: ghx, ghy                              ! global halos

  status = RPN_COMM_ERROR                     ! preset status to failure
  if(setno < 1 .or. setno > iosets) return    ! setno out of bounds, OUCH !
  if(io_set(setno)%ioset .ne. setno) return   ! IO set no longer valid, OUCH !
  if(pe_nx .ne. nx .or. pe_ny .ne. ny) return ! wrong number of PEs in grid, OUCH !
  if(periody) return                          ! no longer supported in this version
  if(periodx) return                          ! no longer supported in this version
! when calling RPN_COMM_shuf_dist_1a, global halos ghx and ghy are needed
  ghx = (gni - (start_x(nx) + count_x(nx) - 1) ) / 2   !  error if not equal to (1 - mini) or 0
  ghy = (gnj - (start_y(ny) + count_y(ny) - 1) ) / 2   !  error if not equal to (1 - minj) or 0
!   print *,'RPN_COMM_shuf_dist: gni, ghx, gnj, ghy =',gni, ghx, gnj, ghy
!
! distribute one 2D plane at a time, one group of IO PEs at a time
! in a group of IO PEs, no column has more than 1 IO PE, neither has any row
!
  do k= 1 , dnk                         ! loop over 2D planes to distribute
    do i= 1 , io_set(setno)%ngroups     ! loop over groups in this set of IO PEs
      low = 1 + (i-1) * io_set(setno)%groupsize      ! index of first PE in goup
      high = min( io_set(setno)%npe , low+io_set(setno)%groupsize-1 )  ! index of last PE in group
!      print *,"DEBUG: group=",i," of",io_set(setno)%ngroups
!      print *,"DEBUG: , PEs x",io_set(setno)%x(low:high)
!      print *,"DEBUG: , PEs y",io_set(setno)%y(low:high)
      listeik = 0
      do n = low, high  ! is this PE part of this group ?
        if(pe_mex == io_set(setno)%x(n) .and. pe_mey == io_set(setno)%y(n)) listeik = liste_i(k)
      enddo
!      print *,"DEBUG: shuf_dist, k,i=",k,i
!     TEMPORARY REVERSION TO "old" VERSION UNTIL SUSPECTED BUG FOUND IN "new" VERSION
      call RPN_COMM_shuf_dist_old(setno, &
                                global(1,1,k), gni, gnj, listeik,  &
                                local, mini, maxi, minj, maxj, lnk,  &
                                liste_o, io_set(setno)%x(low:high),  io_set(setno)%y(low:high), (high-low+1), &
                                start_x, count_x, start_y, count_y,  &
                                periodx, periody, status)
!       call RPN_COMM_shuf_dist_new(setno, &
!                                 global(1,1,k), gni - 2*ghx, gnj - 2*ghy, listeik, ghx, ghy,  &
!                                 local, mini, maxi, minj, maxj, lnk,  &
!                                 liste_o, io_set(setno)%x(low:high),  io_set(setno)%y(low:high), (high-low+1), &
!                                 start_x, count_x, start_y, count_y,  &
!                                 status)
!      print *,"DEBUG: k,i,status,low,high,liste_i(k)=",k,i,status,low,high,liste_i(k)
      if(status == RPN_COMM_ERROR) return
    enddo
  enddo
  contains
!====================================================================================
!
! distribute one 2D array from one group within an IO set
! assumption: no column has 2 IO PES, neither has any row
!
!====================================================================================
! old version, replaced by RPN_COMM_shuf_dist_new
  subroutine RPN_COMM_shuf_dist_old(setno,  &
                                  global,gni,gnj,gnk,  &
                                  local,mini,maxi,minj,maxj,lnk,  &
                                  liste_o, pe_x, pe_y, npes, &
                                  start_x,count_x,start_y,count_y,  &
                                  periodx,periody,status)
    use rpn_comm
    use RPN_COMM_io_pe_tables
    implicit none
!    real*8, external :: RPN_COMM_wtime
    integer, intent(IN) :: setno
    integer, intent(IN) :: gni,gnj,gnk
    integer, intent(IN), dimension(gni,gnj) :: global
    integer, intent(IN) :: mini,maxi,minj,maxj,lnk
    integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
    logical, intent(OUT), dimension(lnk) :: liste_o
    integer, intent(IN) :: npes
    integer, intent(IN), dimension(npes) :: pe_x, pe_y
    integer, intent(IN), dimension(0:pe_nx-1) :: start_x,count_x
    integer, intent(IN), dimension(0:pe_ny-1) :: start_y,count_y
    logical, intent(IN) :: periodx,periody
    integer, intent(OUT) :: status
    logical :: on_column, error, eff_periodx, duplicate
    integer :: root, ybase, kcol, ierr, halox, haloy, lni, lnj, size2d
    integer, dimension(0:pe_nx-1) :: listofk
    integer :: i, j, k, i0, in, ioff, lni0
    integer, dimension(:,:), pointer :: fullrow
    integer, dimension(:,:,:), pointer :: local_1
    integer, dimension(0:pe_nx-1) :: cxs, dxs, cxr, dxr
    integer, dimension(0:pe_ny-1) :: cy, dy
    real*8, dimension(0:6) :: t
    integer, dimension(0:6) :: it

    on_column = .false.     ! precondition for failure in case a hasty exit is needed
    status = RPN_COMM_ERROR ! precondition for failure in case a hasty exit is needed
    nullify(fullrow)
    nullify(local_1)
!
    root = -1 
    kcol = -1
    haloy = 1 - minj          ! halos are implicitly specified by lower bound of x and y dimensions
    if(pe_ny == 1) haloy = 0  ! no halo along y in this case
    halox = 1 - mini
    size2d = (maxj-minj+1) * (maxi-mini+1)    ! size of a 2D local array
    eff_periodx = periodx .and. (halox > 0)   ! global along x perodioc adjustment needed
    lni = count_x(pe_mex)                     ! useful number of points on local tile
    lnj = count_y(pe_mey)
    if(maxi < lni+halox .or. maxj < lnj+haloy) then
      print 101,"ERROR: upper bound of array too small to accomodate halo"
      print 101,"mini,maxi,lni,minj,maxj,lnj",mini,maxi,lni,minj,maxj,lnj
      return  ! OOPS, upper bound cannot accomodate halo
    endif
    do i = 1 , npes
      if(pe_mex == pe_x(i)) then
        on_column = .true.                       ! there is an IO PE on the column (and only one)
        root = pe_y(i)                           ! y coordinate of IO PE, will be the root for scatterv
        if(root == pe_mey) kcol = gnk             ! level of 2D plane that this PE will distribute (it is the root)
        exit
      endif
    enddo
    t(0) = RPN_COMM_wtime()
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
    t(1:6) = t(0)
    if(on_column) call mpi_bcast(kcol,1,MPI_INTEGER,root,pe_mycol,ierr)    ! send gnk to all PEs on the column
    t(1) = RPN_COMM_wtime()
!
!   get list of which PE has which gnk along row
!   if column had no IO PE, -1 gets transmitted
!   if IO PE on column has nothing to contribute, 0 should get transmitted
!
    listofk = 0
    call mpi_allgather(kcol,1,MPI_INTEGER,listofk,1,MPI_INTEGER,pe_myrow,ierr)
    t(2) = RPN_COMM_wtime()
    if(maxval(listofk) > lnk ) then    ! attempt to store a level > lnk, OOPS
      print 101,"ERROR: 1 or more level to distribute > local nk, max level, local nk ",maxval(listofk),lnk
      return
    endif
!print *,"DEBUG: listofk=",listofk
    if(maxval(listofk) <= 0) then   ! no contribution from any IO PE, job id done for this round
!      print 101,"DEBUG: no work to do on this pass"
      status = RPN_COMM_OK
      return
    endif
!   maybe we should check instead if liste_o(listofk(i)) is already .true. which would point to a possible duplicate
    duplicate = .false.
    do i = 0 , pe_nx-1
       if(listofk(i) > 0) then
         if(liste_o(listofk(i))) then
           print 101,"WARNING: this level has already been distributed :",listofk(i)
         endif
         liste_o(listofk(i)) = .false.
       endif
    enddo
!
!   first and last PE on column get count + haloy rows, others get count + 2*haloy rows
!   periody condition to be dealt with later
!
    cy = count_y + haloy
    cy(1:pe_ny-2) = cy(1:pe_ny-2) + haloy       ! count, adjusted for halo (periody assumed false)
!
!   all except first PE on column have their starting point bumped down by one haloy width
!
    dy = start_y - 1                            ! start of each tile
    dy(1:pe_ny-1) = dy(1:pe_ny-1) - haloy       ! displacement, adjusted for halo (periody assumed false)
!
    if(any(dy < 0) .or. any(cy+dy > gnj)) then
      print 101,"ERROR(RPN_COMM_shuf_dist_1): problem with distribution along y. gnj, min(dy), max(cy+dy)",gnj,minval(dy),maxval(cy+dy)
      return
!    else
!     print 101,"INFO(RPN_COMM_shuf_dist_1): distribution along y. gnj, min(dy), max(cy+dy)",gnj,minval(dy),maxval(cy+dy)
    endif
!
    cy = cy * gni                               ! multiply by row length
    dy = dy * gni                               ! multiply by row length
    ybase = minj
    if(pe_mey == 0) ybase = 1                   ! south PE gets no south halo (periody assumed false)
!
    if(on_column .and. kcol > 0)then    ! this PE is on a column where a member of the IO PE group has something to send
      allocate(fullrow(gni,minj:maxj))
      if(minj < 1)    fullrow(:,minj:0) = 0
      if(maxj > lnj)  fullrow(:,lnj+1:maxj) = 0
!print *,'IN shuffle before scatter'
      t(3) = RPN_COMM_wtime()
      call mpi_scatterv(global,cy,dy,MPI_INTEGER,   &          ! 
                        fullrow(1,ybase),cy(pe_mey),MPI_INTEGER, &
                        root,pe_mycol,ierr)
      t(4) = RPN_COMM_wtime()
!print *,'IN shuffle after scatter'
!print *,"DEBUG: ======= fullrow for level",kcol
!do j=maxj,minj,-1
!  print 100,j,fullrow(:,j)
!enddo
!
!     fullrow now contains what will be redistributed along x, haloy is accounted for
!     we may now process reshaping and halo along x periodicity condition
!
      allocate(local_1(mini:maxi,minj:maxj,pe_nx))    ! reshape for distribution along x
      lni0 = count_x(0)
      do k = 1 , pe_nx
        i0 = mini
        in = maxi
        if(k == 1) i0 = 1                        ! lower bound along x of west PE
        if(k == pe_nx) in = count_x(pe_nx-1)     ! upper boung along x of east PE
        ioff = start_x(k-1) - 1                  ! offset along x in global space for PE no k-1 along x
        do j = minj , maxj
          local_1(  :  ,j,k) = 0
          local_1(i0:in,j,k) = fullrow(i0+ioff:in+ioff,j)
          if(eff_periodx .and. (k == 1)) &
                 local_1(    1-halox:0    ,j,k) = fullrow(gni-halox+1:gni,j)   ! west halo from global east
          if(eff_periodx .and. (k == pe_nx)) &
                 local_1(lni0+1:lni0+halox,j,k) = fullrow(1:halox,j)           ! east halo from global west
        enddo
      enddo
                    ! send sizes and displacements for the final alltoallv
      cxs = size2d  ! a full 2D local slice will be sent to all PEs on row
      dxs(0) = 0
      do i = 1 , pe_nx-1
        dxs(i) = dxs(i-1) + cxs(i-1)
      enddo
!print *,"DEBUG: ======= local_1 for level",kcol
!do k = pe_nx,1,-1
!  print *,"======= PE", k
!  do j = maxj,minj,-1
!    print 100,j,local_1(:,j,k)
!  enddo
!enddo
    else                           ! we are not on a column where there is an IO PE
      allocate(fullrow(1,1))
      allocate(local_1(1,1,1))
      cxs = 0                      ! nothing to send from here, size = displacement = 0
      dxs = 0
    endif
!
!   receive sizes and displacements for the final alltoallv
!
    do i = 0 , pe_nx-1
      if(listofk(i) > 0) then    ! something (2D slice) will be received from column i
        cxr(i) = size2d
        dxr(i) = (listofk(i)-1) * size2d   ! offset into local is listofk(i)-1 2D planes
      else                       ! nothing to receive from column i
        cxr(i) = 0
        dxr(i) = 0
      endif
    enddo
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
!
#if defined(FULL_DEBUG)
if(pe_me==0) print *,"DEBUG: kcol,listofk", kcol,listofk
#endif
!print 100,kcol,cxs,dxs,11111,cxr,dxr,11111,listofk
!do k=lnk,1,-1
!  print *,'=== lv=',k
!  do j=maxj,minj,-1
!    print 100,j,local(:,j,k)
!  enddo
!enddo
    t(5) = RPN_COMM_wtime()
    call mpi_alltoallv(local_1, cxs, dxs, MPI_INTEGER,  &
                       local,   cxr, dxr, MPI_INTEGER,  &
                       pe_myrow, ierr)
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
!
    t(6) = RPN_COMM_wtime()
    it(0) = (t(6) - t(0)) * 1000000
    do i = 1,6
      it(i) = max( 0.0_8 , (t(i) - t(i-1)) * 1000000 ) ! convert to microseconds
    enddo
    
!    print 111,"INFO: distribution timings(microsec)",it(0:6)
!
!do k=lnk,1,-1
!  print *,'=== lv=',k
!  do j=maxj,minj,-1
!    print 100,j,local(:,j,k)
!  enddo
!enddo
!print *,"DEBUG: exiting dist_1"
    if(associated(fullrow)) deallocate(fullrow)
    if(associated(local_1)) deallocate(local_1)
    status = RPN_COMM_OK   ! success at last !
    do i = 0 , pe_nx-1     ! mark 2D array at position listofk(i) as received
       if(listofk(i) > 0) liste_o(listofk(i)) = .true.
    enddo
100 format(I3,20I6.5)
101 format(A,20I5)
111 format(A,20I9)
  end subroutine RPN_COMM_shuf_dist_old
  ! ==================================================================================
! if array global is to have halos, they must be IDENTICAL to the halos of array local
  subroutine RPN_COMM_shuf_dist_new(setno,  &
                                  global,gni,gnj,gnk,ghx,ghy,  &
                                  local,mini,maxi,minj,maxj,lnk,  &
                                  liste_o, pe_x, pe_y, npes, &
                                  start_x,count_x,start_y,count_y,  &
                                  status)
    use rpn_comm
    use RPN_COMM_io_pe_tables
    implicit none
    integer, intent(IN) :: setno                       ! IO set number
    integer, intent(IN) :: gni,gnj,gnk                 ! global dimensions (without halo)
    integer, intent(IN) :: ghx,ghy                     ! global halos (may be zero)
    integer, intent(IN), dimension(1-ghx:gni+ghx,1-ghy:gnj+ghy) :: global   ! (gni,gnj) if ho halos
    integer, intent(IN) :: mini,maxi,minj,maxj,lnk     ! local array dimensions
    integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
    logical, intent(OUT), dimension(lnk) :: liste_o    ! output markers for level k processing
    integer, intent(IN) :: npes                        ! total number of IO PEs
    integer, intent(IN), dimension(npes) :: pe_x, pe_y ! x and y coordinates of IO PEs
    integer, intent(IN), dimension(0:pe_nx-1) :: start_x,count_x  !  start in global space and size of tiles along x
    integer, intent(IN), dimension(0:pe_ny-1) :: start_y,count_y  !  start in global space and size of tiles along y
    integer, intent(OUT) :: status
!
    logical :: on_column, error
    integer :: root, ybase, kcol, ierr, halox, haloy, lni, lnj, size2d
    integer, dimension(0:pe_nx-1) :: listofk            ! list of levels collected on this PE row
    integer :: i, j, k, i0, in, ioff, lni0
    integer, dimension(:,:),   pointer :: fullrow       ! local array to collect a full row (global ni + halo,local nj + halo)
    integer, dimension(:,:,:), pointer :: local_1       ! local array to reshape fullrow for distribution along x
    integer, dimension(0:pe_nx-1) :: cxs, dxs, cxr, dxr
    integer, dimension(0:pe_ny-1) :: cy, dy
    real*8, dimension(0:6)  :: t        ! temporary timing collection array
    integer, dimension(0:6) :: it       ! timing collection array

    on_column = .false.     ! precondition for failure in case a hasty exit is needed
    status = RPN_COMM_ERROR ! precondition for failure in case a hasty exit is needed
    nullify(fullrow)
    nullify(local_1)
!
    root = -1 
    kcol = -1
    haloy = 1 - minj                          ! haloy  is inferred from minj
!     if(pe_ny == 1) haloy = 0                  ! no halo along y in this case
    halox = 1 - mini                          ! halox is inferred from mini
    size2d = (maxj-minj+1) * (maxi-mini+1)    ! size of a 2D local array slice
    lni = count_x(pe_mex)                     ! useful number of points along x on local tile
    lnj = count_y(pe_mey)                     ! useful number of points along y on local tile
    if(maxi < lni+halox .or. maxj < lnj+haloy) then
      print 101,"ERROR: upper bound of array is too small to accomodate halo"
      print 101,"       mini,maxi,lni,minj,maxj,lnj",mini,maxi,lni,minj,maxj,lnj
      return  ! OOPS, upper bound cannot accomodate halo
    endif
    ! global halo may be either 0 or equal to local halo
    if( (ghx .ne. halox .and. ghx .ne. 0) .or. (ghy .ne. haloy .and. ghy .ne.0)) then
      print 101,"ERROR: global halos and local halos are not consistent"
      print 101,"       halox, ghx, haloy, ghy",halox, ghx, haloy, ghy
      return  ! OOPS, inconsistent global and local halos
    endif
    do i = 1 , npes
      if(pe_mex == pe_x(i)) then
        on_column = .true.                       ! there is an IO PE on this column (and only one)
        root = pe_y(i)                           ! y coordinate of IO PE, it will be the root for scatterv
        if(root == pe_mey) kcol = gnk            ! level of the 2D plane that this PE will distribute (it is the root)
        exit
      endif
    enddo
    t(0) = RPN_COMM_wtime()
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
    t(1:6) = t(0)
    if(on_column) call mpi_bcast(kcol,1,MPI_INTEGER,root,pe_mycol,ierr)    ! broadcast gnk to all PEs on this column
    t(1) = RPN_COMM_wtime()
!
!   get list of which PE has which level to distribute along this PE row
!   if there is no IO PE in this column, -1 is sent
!   if IO PE on column has nothing to contribute, 0 should get transmitted
!
    listofk = 0
    call mpi_allgather(kcol,1,MPI_INTEGER,listofk,1,MPI_INTEGER,pe_myrow,ierr)
    t(2) = RPN_COMM_wtime()
    if(maxval(listofk) > lnk ) then    ! attempt to store a level > lnk, OOPS
      print 101,"ERROR: 1 or more level to distribute > local nk, highest level index, local nk ",maxval(listofk),lnk
      return
    endif

    if(maxval(listofk) <= 0) then   ! no contribution from any IO PE, job is done for this round
      status = RPN_COMM_OK
      return
    endif
!   maybe we should check instead if liste_o(listofk(i)) is already .true. which would point to a possible duplicate
    do i = 0 , pe_nx-1
       if(listofk(i) > 0) then
         if(liste_o(listofk(i))) then
           print 101,"WARNING: this level has already been distributed :",listofk(i)
         endif
         liste_o(listofk(i)) = .false.
       endif
    enddo
!
!   south and north PE on column get count + haloy + ghy rows, others get count + 2*haloy rows
!
    cy(0)         = count_y(0)         + haloy + ghy     ! south strip
    cy(1:pe_ny-2) = count_y(1:pe_ny-2) + 2 * haloy       ! count, adjusted for halo
    cy(pe_ny-1)   = count_y(pe_ny-1)   + haloy + ghy     ! north strip
!
!   all except south PE on column have their starting point bumped down by one haloy width
!   south PE on column has its starting point bumped down by ghy (global halo)
!
    dy(0)         = start_y(0)         - 1               ! displacement in gloabal array (south row)
    dy(1:pe_ny-1) = start_y(1:pe_ny-1) - 1 - haloy + ghy ! displacement in gloabal array, adjusted for halo (all but south row)
!
    if(any(dy < 0) .or. any(cy+dy > (gnj+2*ghy))) then   ! negative start, or span along y too large
      print 101,"ERROR(RPN_COMM_shuf_dist_1): problem with distribution along y. gnj, min(dy), max(cy+dy)",gnj,minval(dy),maxval(cy+dy)
      return
!    else
!     print 101,"INFO(RPN_COMM_shuf_dist_1): distribution along y. gnj, min(dy), max(cy+dy)",gnj,minval(dy),maxval(cy+dy)
    endif
!
    cy = cy * (gni + ghx*2)                     ! multiply by row length of global array
    dy = dy * (gni + ghx*2)                     ! multiply by row length of global array
    ybase = minj
    if(pe_mey == 0 .and. ghy == 0) ybase = 1    ! south PE may get no south halo
!
    if(on_column .and. kcol > 0)then    ! this PE is on a column where a member of the IO PE group has something to send
      allocate(fullrow(1-ghx:gni+ghx,minj:maxj))
      if(minj < 1)    fullrow(:,minj:0) = 0     ! set halo along j to 0 in local array
      if(maxj > lnj)  fullrow(:,lnj+1:maxj) = 0
      !  scatter full rows across column, including halos
      t(3) = RPN_COMM_wtime()
      call mpi_scatterv(global,cy,dy,MPI_INTEGER,   &          ! 
                        fullrow(1-ghx,ybase),cy(pe_mey),MPI_INTEGER, &
                        root,pe_mycol,ierr)
      t(4) = RPN_COMM_wtime()
!
!     fullrow now contains what will be redistributed along x, haloy is accounted for
!     we may now proceed to reshaping
!
      allocate(local_1(mini:maxi,minj:maxj,pe_nx))    ! reshape for distribution along x
      lni0 = count_x(0)
      do k = 1 , pe_nx
        i0 = mini
        in = maxi
        if(k == 1)     i0 = 1 - ghx                 ! lower bound along x for west PE (1 if no global halo)
        if(k == pe_nx) in = count_x(pe_nx-1) + ghx  ! upper boung along x for east PE (count_x(pe_nx-1) if no global halo)
        ioff = start_x(k-1) - 1                     ! offset along x in global space for PE no k-1 along x
        do j = minj , maxj
          local_1(  :  ,j,k) = 0                    ! make sure skipped elements are zeroed
          local_1(i0:in,j,k) = fullrow(i0+ioff:in+ioff,j)
        enddo
      enddo
                    ! set sizes and displacements for the final alltoallv
      cxs = size2d  ! a full 2D local slice will be sent to all PEs on row
      dxs(0) = 0
      do i = 1 , pe_nx-1
        dxs(i) = dxs(i-1) + cxs(i-1)    ! source displacement
      enddo
    else                           ! we are not on a column where there is an IO PE
      allocate(fullrow(1,1))       ! dummy allocate
      allocate(local_1(1,1,1))     ! dummy allocate
      cxs = 0                      ! nothing to send from here, size = displacement = 0
      dxs = 0
    endif                          ! (on_column .and. kcol > 0)
!
!   receive sizes and displacements for the final alltoallv along x
!
    do i = 0 , pe_nx-1
      if(listofk(i) > 0) then    ! something (2D slice) will be received from column i
        cxr(i) = size2d          ! size of 2D slice
        dxr(i) = (listofk(i)-1) * size2d   ! offset into local is listofk(i)-1 2D planes
      else                       ! nothing to receive from column i
        cxr(i) = 0
        dxr(i) = 0
      endif
    enddo
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
!
#if defined(FULL_DEBUG)
if(pe_me==0) print *,"DEBUG: kcol,listofk", kcol,listofk
#endif
    t(5) = RPN_COMM_wtime()
    call mpi_alltoallv(local_1, cxs, dxs, MPI_INTEGER,  &
                       local,   cxr, dxr, MPI_INTEGER,  &
                       pe_myrow, ierr)
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
!
    t(6) = RPN_COMM_wtime()
    it(0) = (t(6) - t(0)) * 1000000
    do i = 1,6
      it(i) = max( 0.0_8 , (t(i) - t(i-1)) * 1000000 ) ! convert to microseconds
    enddo
!    
    if(associated(fullrow)) deallocate(fullrow)  ! deallocate temporary arrays
    if(associated(local_1)) deallocate(local_1)
!
    status = RPN_COMM_OK   ! success at last !
    do i = 0 , pe_nx-1     ! mark 2D array at position listofk(i) as received
       if(listofk(i) > 0) liste_o(listofk(i)) = .true.
    enddo

100 format(I3,20I6.5)
101 format(A,20I5)
111 format(A,20I9)
  end subroutine RPN_COMM_shuf_dist_new
end subroutine RPN_COMM_shuf_dist
!====================================================================================
!
! collect nk 2D array sections into IO PEs from the grid PEs using IO PE set setno
!
! nk 2D array sections coming in, reassembled on the IO PEs of set setno
! each IO PE ends up with a small number (nk / size of IO PE set) of full arrays
!
! status will be RPN_COMM_OK or RPN_COMM_ERROR
!
! liste_o(n) is set to k when 2D array level k has been received
! start_x(i) contains the first global x index of local array in grid PE (i-1,any)
! count_x(i) contains the number of points along x of local array in grid PE (i-1,any)
! nx is the dimension of start_x and count_x (ERROR if not equal to pe_nx)
! start_y(j) contains the first global y index of local array in grid PE (any,j-1)
! count_y(j) contains the number of points along y of local array in grid (any,j-1)
! ny is the dimension of start_y and count_y (ERROR if not equal to pe_ny)
! mini,maxi,minj,maxj are used to determine the halo width in the local array
!                     this is used to save a halo exchange after distribution
!
! start_x, start_y  are in origin(1)
!
! important notes:
!    the following parameters MUST de the same on ALL PEs
!      setno                    IO set number, as obtained from RPN_COMM_create_io_set
!      mini,maxi,minj,maxj,nk  (dimensions of the "local" array)
!      gni,gnj,dnk             (dimensions of the "global" array)
!      nx,ny,start_x,count_x,start_y,count_y
!
!    even if not used/necessary on a given PE
!      the "global" array must exist ( (1,1,1) array is OK on non IO PEs )
!      array liste_o must exist  ( a dimension(1) array is OK on non IO PEs )
!
!====================================================================================
subroutine RPN_COMM_shuf_coll(setno,  &
                              global,gni,gnj,dnk,  &
                              local,mini,maxi,minj,maxj,nk,  &
                              liste_o,  &
                              start_x,count_x,nx,start_y,count_y,ny,  &
                              status)
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno,gni,gnj,dnk,mini,maxi,minj,maxj,nk,nx,ny
  integer, intent(OUT), dimension(gni,gnj,dnk) :: global
  integer, intent(IN), dimension(mini:maxi,minj:maxj,nk) :: local
  integer, intent(IN), dimension(1:nx)    :: start_x,count_x
  integer, intent(IN), dimension(1:ny)    :: start_y,count_y
  integer, intent(OUT), dimension(dnk)  :: liste_o
  integer, intent(OUT) :: status
  integer :: iset, npass, setsize, igroup, groupsize
  integer k0, k1, kn, low, high
!
  status = RPN_COMM_ERROR
  if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
  if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
  if(pe_nx .ne. nx .or. pe_ny .ne. ny) return ! wrong number of PEs in grid
!
  if(io_set(setno)%me >= 0) liste_o = 0      ! this PE is a member of the IO set
  setsize = io_set(setno)%npe
  groupsize = io_set(setno)%groupsize
  npass = (nk+setsize-1) / setsize
  if(dnk < npass) then   ! OUCH, cannot collect
    print *,"ERROR: cannot collect. setsize, groupsize,npass,dnk,nk",setsize, groupsize,npass,dnk,nk
    return
  endif
  k0 = 1
!  print *,"DEBUG: npass=",npass
!
! loop over the number of passes necessary to distribute nk over setsize IO PEs
! each IO PE will receive (npass) or (npass - 1) full arrays
!
  do iset = 1 , npass   ! process up to setsize levels per iteration
    k1 = k0             ! base level for the first group of this pass
!
!   loop over the groups in the IO PE set
!   the active routine needs the "no column has 2 IO PES, neither has any row" condition
!   IO PEs low -> high in set will potentially receive a full 2D array
!   from levels k1 -> kn   ( levels k1 -> min(nk , k1+groupsize-1) )
!   if more IO PEs in group than levels to distribute some IO PEs will receive nothing
!
    do igroup = 1 , io_set(setno)%ngroups  ! loop over groups in IO PE set
      low = 1 + (igroup-1) * io_set(setno)%groupsize      ! index of first IO PE in goup
      high = min( io_set(setno)%npe , low+groupsize-1 )   ! index of last PE in group
      kn = min( nk , k1+high-low)
#if defined(FULL_DEBUG)
      print *,"DEBUG: shuf_coll, iset, igroup =",iset, igroup
      print *,"DEBUG: shuf_coll, low, high =",low, high
      print *,"DEBUG: shuf_coll, nk, k1, kn =",nk, k1, kn
#endif
      call RPN_COMM_shuf_coll_1(setno,  &
                                global(1,1,iset),gni,gnj,  &
                                local,mini,maxi,minj,maxj,nk,k1,kn,  &
                                liste_o(iset), io_set(setno)%x(low:high), io_set(setno)%y(low:high), (high-low+1), &
                                start_x,count_x,start_y,count_y,  &
                                status)
      k1 = k1 + groupsize  ! base level for next group
    enddo
    k0 = k0 + setsize      ! base level for the first group of next pass
  enddo
!
! collect one 2D plane at a time
! no column has 2 IO PES, neither has any row
!
  contains
!
!====================================================================================
!
! collect one 2D array at a time from a set of 2D array sections
! if kn-k1+1 smaller than groupsize, some IO PEs will not receive anything
! the Mth PE will receive level k1+M-1
! assumption: no column has 2 IO PES, neither has any row
!
!====================================================================================
  subroutine RPN_COMM_shuf_coll_1(setno,  &
                                  global, gni, gnj, &
                                  local, mini, maxi, minj, maxj, nk, k1, kn, &
                                  levnk, pe_x, pe_y, npes, &
                                  start_x, count_x, start_y, count_y, &
                                  status)
!
    use rpn_comm
    use RPN_COMM_io_pe_tables
    implicit none
    integer, intent(IN) :: setno,gni,gnj,mini,maxi,minj,maxj,nk,k1,kn
    integer, intent(OUT), dimension(gni,gnj) :: global
    integer, intent(IN), dimension(mini:maxi,minj:maxj,nk) :: local
    integer, intent(IN) :: npes
    integer, intent(IN), dimension(npes)  :: pe_x, pe_y
    integer, intent(IN), dimension(0:pe_nx-1)    :: start_x,count_x
    integer, intent(IN), dimension(0:pe_ny-1)    :: start_y,count_y
    integer, intent(OUT)  :: levnk
    integer, intent(OUT) :: status
    integer :: i, j, kexpected
    integer, dimension(0:pe_nx-1) :: cxs, dxs, cxr, dxr
    integer, dimension(0:pe_ny-1) :: cy, dy
    integer, dimension(:,:,:), allocatable :: local_1
    integer, dimension(:,:), allocatable :: local_2
    logical :: io_on_column
    integer :: blocki, blockj, k, nlev, klev, ierr, column_root, dimenj
    real*8, dimension(0:6) :: t
    integer, dimension(0:6) :: it
!
    t(0) = RPN_COMM_wtime()
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
!
    t(1:6) = t(0)
    status = RPN_COMM_ERROR
    nlev = kn - k1 + 1             ! number of levels to distribute
    if(nlev > npes) then           ! nlev cannot be larger than npes
      ! add error message
      print *,"ERROR: nlev cannot be larger than npes. nlev,npes=",nlev,npes
      return
    endif
    kexpected = 0
    do i = 1 , npes  ! if this PE is part of the group, flag level
      if( pe_mex == pe_x(i) .and. pe_mey == pe_y(i) ) then
        kexpected = k1 + i - 1   ! expected level on this PE
        if(kexpected <= kn)   levnk = -kexpected    ! if successful, levnk will be +kexpected
      endif
    enddo
    io_on_column = .false.
    column_root = -1
    klev = 0
    do i = 1 , nlev
      if(pe_x(i) == pe_mex) then
        io_on_column = .true.         ! there a useful IO PE on this column
        klev = k1 + i -1              ! level klev will be collected on this column
        column_root = pe_y(i)         ! by PE having pe_mey == column_root 
      endif
    enddo
#if defined(FULL_DEBUG)
      print *,"DEBUG: shuf_coll1, io_on_column, nk, k1, kn =", io_on_column, nk, k1, kn
      print *,"DEBUG: shuf_coll1, kexpected, klev, column_root =",kexpected, klev, column_root
#endif
!
!   PASS 1 , row alltoallv to get a local (gni,lnj) array with the proper level
!            on PEs in the same column as IO PEs
!            PEs where pe_mex = pe_x(l) receive level (k1+l-1) into (mini:maxi,lnj,pe_nx) array
!            PEs on row will send local(:,1:lnj,k1+l-1) to pe_x(l)
!
    blocki = maxi - mini + 1       ! number of points along x in array local
    blockj = count_y(pe_mey)       ! number of useful points along y in row pe_mey
    dimenj = maxj - minj + 1
!
    if(io_on_column) then               ! this PE will be collecting one level for this row
      allocate( local_1(mini:maxi,blockj,0:pe_nx-1) )
      local_1 = 77777
      allocate( local_2(gni,blockj) )
    else
      allocate( local_1(1,1,1) )
      allocate( local_2(1,1) )
    endif
!
    cxs = 0
    dxs = 0
    do i = 1 , nlev
      j = pe_x(i)                ! PE j will receive level k1 + (i-1)
      cxs(j) =  blocki * blockj
      dxs(j) = (blocki * dimenj) * (k1 + i - 2)  ! offset of level  k1 + (i-1)
      dxs(j) = dxs(j) + blocki * (1-minj)        ! correct for halo along y
    enddo
!
    cxr = 0
    dxr = 0
    if(klev > 0) then          ! this PE collects level klev from all PEs in row
      cxr = blocki * blockj    ! x y 2D block size
      do i = 0 , pe_nx -1
        dxr(i) = (blocki * blockj) * i
      enddo
    endif
!
#if defined(FULL_DEBUG)
    print *,"DEBUG: ========= local ============="
    do k = kn , k1 , -1
      do j = maxj , minj , -1
        print 101,k,j,local(:,j,k)
      enddo
    enddo
    print *,"DEBUG: blocki, blockj",blocki, blockj
    print *,"DEBUG: cxs",cxs
    print *,"DEBUG: dxs",dxs
    print *,"DEBUG: cxr",cxr
    print *,"DEBUG: dxr",dxr
#endif
    t(1) = RPN_COMM_wtime()
    call mpi_alltoallv(local,   cxs, dxs, MPI_INTEGER,  &   ! send from local, size cxs, displacement dxs
                       local_1, cxr, dxr, MPI_INTEGER,  &   ! receive into local_1, 
                       pe_myrow, ierr)
    t(2) = RPN_COMM_wtime()
!
!   REORG    move from (mini:maxi,lnj,pe_nx) to (gni,lnj)
    if(io_on_column) then
#if defined(FULL_DEBUG)
      print *,"DEBUG: ======== local_1 ========="
#endif
      do k = 0 , pe_nx-1
#if defined(FULL_DEBUG)
        do j = blockj , 1 , -1
          print 101,k,j,local_1(:,j,k)
        enddo
#endif
        do j = 1 , blockj
          local_2( start_x(k):start_x(k)+count_x(k)-1 , j ) = local_1( 1:count_x(k) , j , k )
        enddo
      enddo
#if defined(FULL_DEBUG)
      print *,"DEBUG: ======== local_2 ========="
      do j = blockj , 1 , -1
        print 100,j,local_2(:,j)
      enddo
#endif
    endif
    deallocate( local_1 )   ! local_1 no longer useful
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
    t(3) = RPN_COMM_wtime()
!
!
!   PASS 2 , on columns where there is an IO PE, gatherv on IO PE into (gni,gnj) array
!
    if(io_on_column) then   ! column root will collect level  klev from local_2 into global
       cy = 0
       dy = 0
       if(column_root == pe_mey) then  ! this PE is the gather root
         cy = gni * count_y
         dy = gni * (start_y -1)
       endif
#if defined(FULL_DEBUG)
       print *,"DEBUG: ======== before gatherv on",column_root," ========="
       print *,"DEBUG: cy",cy
       print *,"DEBUG: dy",dy
#endif
    t(4) = RPN_COMM_wtime()
       call mpi_gatherv(local_2, gni*blockj, MPI_INTEGER, &
                        global , cy,  dy,  MPI_INTEGER, &
                        column_root, pe_mycol, ierr)
    t(5) = RPN_COMM_wtime()
       if(kexpected .ne. 0) levnk = kexpected
    endif
!
    call mpi_barrier(pe_indomm,ierr)              ! full sync, start / middle / end
!
    t(6) = RPN_COMM_wtime()
    it(0) = (t(6) - t(0)) * 1000000 
    do i = 1,6
      it(i) = max( 0.0_8 , (t(i) - t(i-1)) * 1000000 ) ! convert to microseconds
    enddo
    
!    print 111,"INFO: collect timings(microsec)",it(0:6)
!
    deallocate( local_2 )
    status = RPN_COMM_OK
100 format(I3,20I6.5)
101 format(2I3,20I6.5)
111 format(A,20I9)
  end subroutine RPN_COMM_shuf_coll_1
end subroutine RPN_COMM_shuf_coll
!
