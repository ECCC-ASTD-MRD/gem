function model_timeout_alarm(seconds) result(seconds_since)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: seconds
  integer :: seconds_since

  interface
  function c_alarm(seconds) result(seconds_since) BIND(C,name='alarm')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), value :: seconds
    integer(C_INT) :: seconds_since
  end function c_alarm
  end interface

  seconds_since = c_alarm(seconds)
  return
end function model_timeout_alarm
