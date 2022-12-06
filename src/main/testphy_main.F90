subroutine testphy_main()
  ! Test stub top-level for the RPN physics package
  use testphy, only: PT_TEST_TYPES, pt_run
  implicit none

#include <rmnlib_basics.hf>

  ! Internal variables
  integer :: istat
  character(len=2048) :: name,test_type,usage

  ! Retrieve and check options
  call get_command_argument(0,value=name)
  usage = 'Usage: '//trim(name)//' [TYPE]'
  call get_command_argument(1,value=test_type,status=istat)
  if (istat > 0) test_type = 'full'
  if (.not.any(PT_TEST_TYPES == test_type)) then
     write(RMN_STDERR,*) trim(usage)
     stop
  endif

  ! Call test stub driver
  istat = pt_run(test_type)
  if (istat /= RMN_OK) then
     write(RMN_STDERR,*) '(testphy): Error returned by pt_run'
     stop
  endif

end subroutine testphy_main
