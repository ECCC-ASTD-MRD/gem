#if defined(TEST)
program my_test
call mpi_init(ierr)
call rpn_comm_test_012
call mpi_finalize(ierr)
stop
end
#endif
subroutine rpn_comm_test_012
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'

  type(rpncomm_datatype) :: my_type
  type(rpncomm_group) :: my_group
  type(rpncomm_window) :: my_win
  type(rpncomm_communicator) :: my_com
  type(rpncomm_file) :: my_file
  type(rpncomm_request) :: my_req
  type(rpncomm_info) :: my_info
  type(rpncomm_operator) :: my_op
  integer(C_INT) :: result

  my_type = NULL_rpncomm_datatype
  my_type%t1 = -99999
  my_type%t2 = -99999
  my_type%p = RPN_COMM_f2c(my_type,MPI_INTEGER)
  result = RPN_COMM_c2f(my_type)
  print 100,'type result =',result,MPI_INTEGER,my_type%t2

  my_group = NULL_rpncomm_group
  my_group%t1 = -99999
  my_group%t2 = -99999
  my_group%p = RPN_COMM_f2c(my_group,MPI_GROUP_NULL)
  result = RPN_COMM_c2f(my_group)
  print 100,'group result =',result,MPI_GROUP_NULL,my_group%t2

  my_win = NULL_rpncomm_window
  my_win%t1 = -99999
  my_win%t2 = -99999
  my_win%p = RPN_COMM_f2c(my_win,MPI_WIN_NULL)
  result = RPN_COMM_c2f(my_win)
  print 100,'window result =',result,MPI_WIN_NULL,my_win%t2

  my_com = NULL_rpncomm_communicator
  my_com%t1 = -99999
  my_com%t2 = -99999
  my_com%p = RPN_COMM_f2c(my_com,MPI_COMM_WORLD)
  result = RPN_COMM_c2f(my_com)
  print 100,'com result =',result,MPI_COMM_WORLD,my_com%t2

!  my_file = NULL_rpncomm_file
!  result = RPN_COMM_c2f(my_file)
!  print 100,'file result =',result

  my_req = NULL_rpncomm_request
  my_req%t1 = -99999
  my_req%t2 = -99999
  my_req%p = RPN_COMM_f2c(my_req,MPI_REQUEST_NULL)
  result = RPN_COMM_c2f(my_req)
  print 100,'request result =',result,MPI_REQUEST_NULL,my_req%t2

  my_info = NULL_rpncomm_info
  my_info%t1 = -99999
  my_info%t2 = -99999
  my_info%p = RPN_COMM_f2c(my_info,MPI_INFO_NULL)
  result = RPN_COMM_c2f(my_info)
  print 100,'info result =',result,MPI_INFO_NULL,my_info%t2

  my_op = NULL_rpncomm_operator
  my_op%t1 = -99999
  my_op%t2 = -99999
  my_op%p = RPN_COMM_f2c(my_op,MPI_SUM)
  result = RPN_COMM_c2f(my_op)
  print 100,'operator result =',result,MPI_SUM,my_op%t2
100 format(A20,3I5)
  return
end subroutine rpn_comm_test_012