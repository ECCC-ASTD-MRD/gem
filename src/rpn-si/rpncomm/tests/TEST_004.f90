!program test_004
!  call rpn_comm_test_004
!  stop
!end
subroutine rpn_comm_test_004  ! test for two player MPMD update game
  implicit none
  include 'RPN_COMM_constants.inc'
  external sss
  integer npex, npey, ierr, Pelocal, Petotal
!
  integer RPN_COMM_init
  integer :: rank_in_grid, rank_in_universe, size_of_universe, mygrid, my_gridsize
  integer :: i, iter, status, increment, other_pe0
  integer, dimension(3) :: my_grid
  integer, pointer, dimension(:,:) :: the_grids
  integer, pointer, dimension(:) :: g_table, l_table
  integer :: table_entry
  character (len=10) :: my_id
!
!       force a "callback" type initialization
!
  npex=0
  npey=0
  mygrid = RPN_COMM_init(sss,Pelocal,Petotal,npex,npey)
  call RPN_COMM_size(RPN_COMM_UNIVERSE,size_of_universe,ierr)
  call RPN_COMM_rank(RPN_COMM_UNIVERSE,rank_in_universe,ierr)
  call RPN_COMM_size(RPN_COMM_GRID,my_gridsize,ierr)  ! my_gridsize must be equal to npex
  call RPN_COMM_rank(RPN_COMM_GRID,rank_in_grid,ierr)
  
  open(10,file='id.txt',status='OLD',err=111)
  read(10,*)my_id
111 continue
  my_grid(1)=1 ; my_grid(2)=rank_in_grid ; my_grid(3)=rank_in_universe ; increment=1
  if(trim(my_id)=='model') then 
    my_grid(1)=2  ! model is player no 2, the other is player no 1
    increment=-2
  endif
  other_pe0=-1

  allocate(the_grids(3,size_of_universe))
  call RPN_COMM_allgather(my_grid,3,RPN_COMM_INTEGER,the_grids,3,RPN_COMM_INTEGER,RPN_COMM_UNIVERSE,ierr)
  if(rank_in_universe==0) then
    write(6,101)the_grids
101 format(3I5)
  endif
  if(rank_in_grid==0)then ! find other PE 0
    do i=1,size_of_universe
      if(the_grids(2,i)==0 .and. the_grids(1,i)/=my_grid(1)) other_pe0=the_grids(3,i)
    enddo
  endif
!
  do i=0,size_of_universe-1
    call RPN_COMM_barrier(RPN_COMM_UNIVERSE,ierr)
    if(rank_in_universe==i) &
       print *,'RANK in universe=',rank_in_universe,' RANK in grid=',rank_in_grid,' ID=',trim(my_id), &
            my_gridsize,npex,npey
    call RPN_COMM_barrier(RPN_COMM_UNIVERSE,ierr)
  enddo
! and now let's play the game
!
  allocate(g_table(0:size_of_universe-1),l_table(0:my_gridsize-1))
  g_table=100
  do iter=1,10
!
    if(rank_in_grid==0) l_table=g_table(rank_in_universe:rank_in_universe+my_gridsize-1)
    call RPN_COMM_bcast(l_table,my_gridsize,RPN_COMM_INTEGER,0,RPN_COMM_GRID,ierr)
!
    table_entry=l_table(rank_in_grid)+increment
    call RPN_COMM_gather(table_entry,1,RPN_COMM_INTEGER,l_table,1,RPN_COMM_INTEGER,0,RPN_COMM_GRID,ierr)
!
    if(rank_in_grid==0) then ! the model and assim PE 0 only
      if(rank_in_universe ==0) then
        g_table(rank_in_universe:rank_in_universe+my_gridsize-1)=l_table ! update my section of global table
        print *,trim(my_id),rank_in_universe,' sending to', other_pe0,' tag 0'
        call RPN_COMM_send(g_table,size_of_universe,RPN_COMM_INTEGER,other_pe0,0,RPN_COMM_UNIVERSE,ierr)  ! send global update to other PE 0
        print *,trim(my_id),rank_in_universe,' receiving from', other_pe0,' tag 1'
        call RPN_COMM_recv(g_table,size_of_universe,RPN_COMM_INTEGER,other_pe0,1,RPN_COMM_UNIVERSE,status,ierr)  ! get global update from other PE 0
      else
        print *,trim(my_id),rank_in_universe,' receiving from', other_pe0,' tag 0'
        call RPN_COMM_recv(g_table,size_of_universe,RPN_COMM_INTEGER,other_pe0,0,RPN_COMM_UNIVERSE,status,ierr)  ! get global update from other PE 0
        g_table(rank_in_universe:rank_in_universe+my_gridsize-1)=l_table ! update my section of global table
        print *,trim(my_id),rank_in_universe,' sending to', other_pe0,' tag 1'
        call RPN_COMM_send(g_table,size_of_universe,RPN_COMM_INTEGER,other_pe0,1,RPN_COMM_UNIVERSE,ierr)  ! send global update to other PE 0
      endif
    endif
  enddo
  if(rank_in_universe==0) then
    write(6,102)g_table
102 format(I5)
  endif
!
777  call RPN_COMM_finalize(ierr)
  stop
end
!
SUBROUTINE sss(nx,ny)
  implicit none 
  integer nx,ny
!
!        "callback routine" used to get initial topology
!        information (very passive version)
!
  return
end
