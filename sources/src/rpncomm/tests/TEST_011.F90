subroutine rpn_comm_test_011
  use rpn_comm
  implicit none
  integer, external :: rpn_comm_2dgrid_test
  integer :: status, params

  pe_nx = 6
  pe_ny = 7
  pe_me = 0
  status = rpn_comm_2dgrid_test(1,params)
end subroutine rpn_comm_test_011
