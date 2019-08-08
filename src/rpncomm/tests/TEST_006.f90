subroutine rpn_comm_test_006
use rpn_comm_io, ONLY : RPN_COMM_file_copy_test
implicit none
integer status
status = RPN_COMM_file_copy_test()
return
end
