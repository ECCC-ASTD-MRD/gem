module RPN_COMM_io  ! file replicator across nodes , async file copy, test and helper routines
  use iso_c_binding
  implicit none

  logical, save :: rpn_comm_io_debug=.false.

  public :: rpn_comm_open

  interface
    integer(C_INT) function c_rpn_comm_unlink(name) BIND(C,name='f_RPN_COMM_unlink')  ! interface to libc open function
      use iso_c_binding
      type(C_PTR), value :: name
    end function c_rpn_comm_unlink

    integer(C_INT) function c_rpn_comm_open(name, mode) BIND(C,name='f_RPN_COMM_open')  ! interface to libc open function
      use iso_c_binding
      integer(C_INT), intent(IN), value :: mode
      type(C_PTR), value :: name
    end function c_rpn_comm_open

    integer(C_LONG_LONG) function rpn_comm_read(fd,buffer,nbytes) BIND(C,name='f_RPN_COMM_read')  ! interface to libc read function
    use iso_c_binding
    integer(C_INT), intent(IN), value :: fd
    integer(C_LONG_LONG), intent(IN), value :: nbytes
    type(C_PTR), value :: buffer
    end function rpn_comm_read

    integer(C_LONG_LONG) function rpn_comm_write(fd,buffer,nbytes) BIND(C,name='f_RPN_COMM_write')  ! interface to libc write function
    use iso_c_binding
    integer(C_INT), intent(IN), value :: fd
    integer(C_LONG_LONG), intent(IN), value :: nbytes
    type(C_PTR), value :: buffer
    end function rpn_comm_write

    integer(C_INT) function rpn_comm_close(fd) BIND(C,name='f_RPN_COMM_close')  ! interface to libc close function
    use iso_c_binding
    integer(C_INT), intent(IN), value :: fd
    end function rpn_comm_close

    integer(C_INT) function rpn_comm_wait(id) BIND(C,name='f_RPN_COMM_wait') ! wait for asynchronous copy with tag = id to terminate
    use iso_c_binding
    integer(C_INT), intent(IN), value :: id
    end function rpn_comm_wait

    integer(C_INT) function rpn_comm_copy(fd1,fd2,asynchronous) BIND(C,name='f_RPN_COMM_copy')  ! C descriptor copy, possibly asunchronous
    use iso_c_binding
    integer(C_INT), intent(IN), value :: fd1,fd2,asynchronous
    end function rpn_comm_copy
  end interface

  contains

  integer function rpn_comm_unlink(name) ! interface to interface, handle Fortran char string to C char string conversion
!
! delete a file
!
  implicit none
  character (len=*), intent(IN) ::name
  character (len=1), dimension(len_trim(name)+1), target :: temp
  
  temp = transfer( trim(name)//achar(0) , temp )
  rpn_comm_unlink=c_rpn_comm_unlink(c_loc(temp))
  end function rpn_comm_unlink

  integer function rpn_comm_open(name, mode) ! interface to interface, handle Fortran char string to C char string conversion
!
! open a file for reading or writing
! mode = 0 : read
! mode = 1 : write
!
  implicit none
  character (len=*), intent(IN) ::name
  integer, intent(IN) :: mode
  character (len=1), dimension(len_trim(name)+1), target :: temp
  
  temp = transfer( trim(name)//achar(0) , temp )
  rpn_comm_open=c_rpn_comm_open(c_loc(temp),mode)
  end function rpn_comm_open

  integer function RPN_COMM_file_copy_test()  ! test asynchronous copying in background and file replication across nodes
  implicit none
  include 'mpif.h'
  integer :: fd1, fd2, status, id, ierr, rank
  integer, external :: RPN_COMM_file_copy_start, RPN_COMM_file_bcst

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)                         ! my rank in world
  rpn_comm_io_debug=.true.

  RPN_COMM_file_copy_test = -1
  status = rpn_comm_copy(-1,-1,1)  ! set verbose debug mode

  if(rank /= 0) goto 222  ! first part of test for PE 0 only, asynchronous copying in background 

  print *,"START OF COPY TEST"
  fd1 = rpn_comm_open("existing_input_file",0) ! open for read
  if(fd1 < 0) then
    print *,"ERROR: cannot open existing_input_file for read or open error"
    return
  endif
  fd2 = rpn_comm_open("new_output_file",0) ! open for read
  if(fd2 >0 ) then
    print *,"WARNING: new_output_file exists"
  endif
  fd2 = rpn_comm_open("new_output_file",1) ! open for write
  if(fd2 <0 ) then
    print *,"ERROR: failed to open new_output_file for write"
    return
  endif
!
  status = rpn_comm_copy(fd1,fd2,0)  ! synchronous copy test
  if(status /= 0) then
    print *,"ERROR: synchronous copy failed"
    return
  endif
!
  id = rpn_comm_file_copy_start("existing_input_file","new_output_file_2")  ! asynchronous copy test
  print *,"INFO: asynchronous copy id =",id
  status = rpn_comm_wait(id)
  if(status /= 0) then
    print *,"ERROR: asynchronous copy failed"
    return
  endif
  print *,"INFO: cmp new_output_file new_output_file_2"
  call system("cmp new_output_file new_output_file_2 && echo ' INFO: files compare OK'")
!
  status = rpn_comm_unlink("new_output_file")  ! get rid of created files
  if(status /= 0) then
    print *,"WARNING: failed to delete new_output_file"
  endif
  status = rpn_comm_unlink("new_output_file_2")
  if(status /= 0) then
    print *,"WARNING: failed to delete new_output_file_2"
  endif
  print *,"END OF COPY TEST"

222 call mpi_barrier(MPI_COMM_WORLD,ierr)  ! all PEs wait for PE 0 to start replication test

  if(rank == 0) print *,"START OF REPLICATION TEST target_dir/existing_input_file"
  status = RPN_COMM_file_bcst("target_dir/existing_input_file",MPI_COMM_WORLD)
  if(rank == 0) print *,"INFO: replication status =",status
  if(rank == 0) print *,"END OF REPLICATION TEST"
  call mpi_finalize(ierr)
  return
  end function RPN_COMM_file_copy_test

end module RPN_COMM_io

integer function RPN_COMM_file_copy_start(name1,name2)  ! start asynchronous file copy
  use iso_c_binding
  use rpn_comm_io
  implicit none
  character (len=*), intent(IN) :: name1, name2   ! copy file name1 into file name2
  integer :: fd1, fd2

  RPN_COMM_file_copy_start = 0  ! O.K. status
  fd1 = rpn_comm_open(name1,0)  ! open input file, get C file descriptor fd1
  if(fd1 < 0) goto 999
  fd2 = rpn_comm_open(name2,1)  ! open output file, get C file descriptor fd2
  if(fd2 < 0) goto 999

  RPN_COMM_file_copy_start = rpn_comm_copy(fd1,fd2,1)   ! copy fd1 onto fd2 asynchronously
  
  return

999 continue   ! error exception return
  RPN_COMM_file_copy_start = -1  ! error status
  return
end function RPN_COMM_file_copy_start

integer function RPN_COMM_file_bcst(name,com)
!
! replicate a file across all the hosts in a communication domain
! PE 0 of domain reads the file in chunks and broadcasts the chunks and their length
! one PE per host then writes the chunk into the file ( except PE 0 on hos hosts)
!
! this function is normally called after PE 0 has read the file from shared storage and written it to 
! fast local storage on its host
!
! this is deemed more efficient than having all processes reading the same file on shared storage
! there will be a few processes on the same host reading the same local file subsequently
!
  use iso_c_binding
  use rpn_comm_io
  implicit none
  character (len=*), intent(IN) :: name   ! name of the file to replicate
  integer, intent(IN) :: com              ! communicator

  interface
    integer(C_INT) function c_gethostid()BIND(C,name='gethostid')  ! interface to libc gethostid
      use iso_c_binding
    end function c_gethostid
  end interface

  include 'mpif.h'

  integer(C_LONG_LONG), parameter :: NW8 = 1024*1024
  integer, dimension(0:NW8), target :: buf  ! 4 MegaByte buffer
  integer :: fd, status, rank, rank_on_host, ierr, my_color, same_host, errors, sum_errors
  integer(C_LONG_LONG) :: nwr, nww

  call mpi_comm_rank(com,rank,ierr)                         ! my rank in communicator
  my_color = abs(c_gethostid())                             ! coloring by host identifier
  call MPI_COMM_SPLIT(com,my_color,rank,same_host,ierr)     ! same host communicator
  call mpi_comm_rank(same_host,rank_on_host,ierr)           ! my rank on host

  errors = 0
  fd = 0   ! must initialize because some PEs may not call open
  if(rank == 0) then
    fd = rpn_comm_open(name,0)   ! PE 0, open for read (file MUST exist)
    if(rpn_comm_io_debug) print 111,"rank=",rank," read  fd=",fd," file=",trim(name)," hostid=",my_color
  else if(rank_on_host ==0) then ! open for write, one process per host, but not if same host as PE 0
    fd = rpn_comm_open(name,1) 
    if(rpn_comm_io_debug) print 111,"rank=",rank," write fd=",fd," file=",trim(name)," hostid=",my_color
  else                           ! nothing to do on this PE
    if(rpn_comm_io_debug) print 111,"rank=",rank," noop  fd=",fd," file=",trim(name)," hostid=",my_color
  endif
  if(fd < 0) errors = errors + 1
  call mpi_allreduce(errors,sum_errors,1,MPI_INTEGER,MPI_SUM,com,ierr)  !  get global eror count
  if(sum_errors > 0) goto 999  ! on open failed somewhere

  nwr = 1
  do while(nwr > 0)
    if(rank == 0) buf(0) = rpn_comm_read(fd,c_loc(buf(1)),NW8*4)   ! PE ranked at 0 reads 
    call mpi_bcast(buf,NW8+1,MPI_INTEGER,0,com,ierr)     ! and broadcasts (this is lazy, but intra host broadcast deemed cheap)
    nwr = buf(0)                                                   ! valid everywhere after broadcast
    nww = buf(0)                                                   ! needed for PEs that do not write to no get errors
    if(rank /= 0 .and. rank_on_host == 0 .and. nwr > 0) then
      if(rpn_comm_io_debug) print 111,"DEBUG: rank no",rank," writing"
      nww = rpn_comm_write(fd,c_loc(buf(1)),nwr)  ! one PE per node writes
    endif
    if(nww /= nwr) errors = errors + 1                             ! not everything written, OOPS
  enddo
  call mpi_allreduce(errors,sum_errors,1,MPI_INTEGER,MPI_SUM,com,ierr)  ! OOPS anywhere ?
  if(sum_errors > 0) goto 999  ! on open failed somewhere

  status = rpn_comm_close(fd)
  RPN_COMM_file_bcst = status
  return

111 format(A,I4,A,I4,A,A,A,Z9)
999 continue   ! error exception
  if(rank == 0) print 111,"ERROR: RPN_COMM_file_bcst errors =", sum_errors," rank =",rank
  RPN_COMM_file_bcst = sum_errors
  return

end function RPN_COMM_file_bcst
