subroutine rpn_comm_test_009
! call RPN_COMM_haloflip_test(halox, haloy, ni, nj, nk)
  real*8 :: t1,t2,mpi_wtime
  integer :: msec
  external :: mpi_wtime
  call mpi_init(ierr)
  call RPN_COMM_haloflip_test(1,     2,      5,  6,  1)
  t1 = mpi_wtime()
  call RPN_COMM_haloflip_test(3,     3,     75, 27, 80)
  t2 = mpi_wtime()
  msec=(t2-t1)*1000000
  print *,'time=',msec,' microseconds'
  t1 = mpi_wtime()
  call RPN_COMM_haloflip_test(3,     3,     75, 27, 80)
  t2 = mpi_wtime()
  msec=(t2-t1)*1000000
  print *,'time=',msec,' microseconds'
  call mpi_finalize(ierr)
  stop
end
