module rpn_comm_types
use ISO_C_BINDING
include 'RPN_COMM_types.inc'

interface rpncomm_associated
  module procedure rpncomm_context_associated
  module procedure rpncomm_communicator_associated
  module procedure rpncomm_datatype_associated
  module procedure rpncomm_operator_associated
  module procedure rpncomm_group_associated
  module procedure rpncomm_ptr_associated
end interface

private :: marker

contains
function rpncomm_context_associated(object) result(value)
  implicit none
  type(rpncomm_context), intent(IN) :: object
  logical :: value
  value = c_associated(object%p)
end function rpncomm_context_associated

function rpncomm_communicator_associated(object) result(value)
  implicit none
  type(rpncomm_communicator), intent(IN) :: object
  logical :: value
  value = c_associated(object%p)
end function rpncomm_communicator_associated

function rpncomm_datatype_associated(object) result(value)
  implicit none
  type(rpncomm_datatype), intent(IN) :: object
  logical :: value
  value = c_associated(object%p)
end function rpncomm_datatype_associated

function rpncomm_operator_associated(object) result(value)
  implicit none
  type(rpncomm_operator), intent(IN) :: object
  logical :: value
  value = c_associated(object%p)
end function rpncomm_operator_associated

function rpncomm_group_associated(object) result(value)
  implicit none
  type(rpncomm_group), intent(IN) :: object
  logical :: value
  value = c_associated(object%p)
end function rpncomm_group_associated

function rpncomm_ptr_associated(object) result(value)
  implicit none
  type(rpncomm_ptr), intent(IN) :: object
  logical :: value
  value = c_associated(object%p)
end function rpncomm_ptr_associated

subroutine marker
call MARKER
return
end subroutine marker

end module rpn_comm_types
