subroutine rpn_comm_test_005
implicit none
include 'mpif.h'
integer, parameter :: BUFSIZE=4096*1024
integer, dimension(BUFSIZE) :: rbuf, sbuf
integer :: ierr, increment, i, me, other, irep, irep2
integer :: worldsize, worldrank
integer, allocatable, dimension(:) :: pe_table
real*8 :: t1, t2, best, worst
integer, dimension(MPI_STATUS_SIZE) :: status

call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,worldsize,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,worldrank,ierr)
worldsize=worldsize-mod(worldsize,2) ! make sure worldsize is even
if(worldrank==0) print *,'worldsize=',worldsize

allocate(pe_table(0:worldsize-1))

do irep2=1,3

  best=999999.0
  worst=0.0
  
  do increment=1,worldsize-1,2
    if(increment > 1 .and. mod(worldsize,increment)==0) goto 111 ! i not prime with respect to worldsize
    pe_table = 0 ! redo pe_table
    do i=1,worldsize-1
      pe_table(i)=mod(pe_table(i-1)+increment,worldsize)
    enddo
    do i=1,worldsize-1
      if(pe_table(i)==0) goto 111 ! permutation is incomplete, uninitialized entries in table
    enddo
    if(worldrank==0) print *,'increment=',increment
    me=-1
    do i=0,worldsize-1
      if(pe_table(i) == worldrank) then
        me = i
        other = pe_table( mod(i+worldsize/2,worldsize) )
        exit
      endif
    enddo
    !print *,'me=',worldrank,' other=',other
    t1=MPI_wtime()
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    do irep=1,10
      if(worldrank<worldsize) call mpi_sendrecv(sbuf,BUFSIZE,MPI_INTEGER,other,1,rbuf,BUFSIZE,MPI_INTEGER,other,1,MPI_COMM_WORLD,status,ierr)
    enddo
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    t2=MPI_wtime()
    best=min(best,t2-t1)
    worst=max(worst,t2-t1)
    if(worldrank==0)then
      print *,"time=",t2-t1," for",BUFSIZE*8," bytes, speed=",(BUFSIZE*4)/(t2-t1)*.00001*worldsize," MB/sec"
    endif
  111 continue
  enddo
  if(worldrank==0)then
    print *,"worst time=",worst," for",BUFSIZE*8," bytes, speed=",(BUFSIZE*4)/(worst)*.00001*worldsize," MB/sec"
    print *,"best time=",best," for",BUFSIZE*8," bytes, speed=",(BUFSIZE*4)/(best)*.00001*worldsize," MB/sec"
  endif
  
enddo

call mpi_finalize(ierr)
stop
end
