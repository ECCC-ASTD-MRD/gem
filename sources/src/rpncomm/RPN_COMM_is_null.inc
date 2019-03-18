interface RPN_COMM_is_null
  function RPN_COMM_null_context(what) result(status)
    import
    type(rpncomm_context), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_context

  function RPN_COMM_null_window(what) result(status)
    import
    type(rpncomm_window), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_window

  function RPN_COMM_null_communicator(what) result(status)
    import
    type(rpncomm_communicator), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_communicator

  function RPN_COMM_null_info(what) result(status)
    import
    type(rpncomm_info), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_info

  function RPN_COMM_null_file(what) result(status)
    import
    type(rpncomm_file), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_file

  function RPN_COMM_null_request(what) result(status)
    import
    type(rpncomm_request), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_request

  function RPN_COMM_null_datatype(what) result(status)
    import
    type(rpncomm_datatype), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_datatype

  function RPN_COMM_null_operator(what) result(status)
    import
    type(rpncomm_operator), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_operator

  function RPN_COMM_null_group(what) result(status)
    import
    type(rpncomm_group), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_group

  function RPN_COMM_null_ptr(what) result(status)
    import
    type(rpncomm_ptr), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_ptr

  function RPN_COMM_null_pattern(what) result(status)
    import
    type(rpncomm_pattern), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_pattern

  function RPN_COMM_null_array(what) result(status)
    import
    type(rpncomm_array), intent(IN) :: what
    logical :: status
  end function RPN_COMM_null_array
end interface

