subroutine time2sec_main()
   use, intrinsic :: iso_fortran_env, only: REAL64
  use timestr_mod, only: timestr2sec
  use str_mod, only: str_toreal
  implicit none
  !  Argument-free interface for time2sec time string converter
#include <rmnlib_basics.hf>

  !  Internal parameterrs
  integer, parameter :: STDOUT=6,STDERR=0,STRLEN=1024

  !  CCARD preparation
  integer, parameter :: NARGS=2
  character(len=STRLEN), dimension(NARGS) :: &
       arglist, &
       def, &
       val

  !  Internal variables
  integer :: npos,istat
  real :: sec,dt
  real(REAL64) :: dt_8
  character(len=STRLEN) :: timestr,dtstr

  !  Obtain arguments
  arglist(1) = 'time.'
  arglist(2) = 'dt.'
  def(1) = 'UNDEFINED'
  def(2) = '1.'
  val(1) = 'UNDEFINED'
  val(2) = '1.'
  call ccard(arglist,def,val,NARGS,npos)
  timestr = val(1)
  dtstr   = val(2)
  
  !  Enforce mandatory arguments
  if (timestr == def(1)) then
     write(STDERR, '(a)') "Error: a time string (-time) must be provided "//trim(timestr)
     return
  endif

  istat = str_toreal(dt,dtstr)
  if (.not.RMN_IS_OK(istat)) then
     write(STDERR, '(a)') "Error: invalid number provided for delta_T (-dt): "//trim(dtstr)
     return
  endif
  dt_8 = dt

  !  Call library function to retrieve number of seconds from time string
  istat = timestr2sec(sec,timestr,dt_8)
  if (.not.RMN_IS_OK(istat)) then
     write(STDERR, '(a)') "Error: unable to convert date string "//trim(timestr)
     return
  else
     write(STDOUT, '(i0)') sec
  endif
  
  return
end subroutine time2sec_main
