subroutine ccpl_timeflip(F_moins,F_plus,F_n,F_nk)

  implicit none
!!!#include <arch_specific.hf>

  ! Timeflip for column coupling:
  !   F_moins = F_plus

  ! Argument declarations
  integer :: F_n                                !Horizontal dimension
  integer :: F_nk                               !Number of levels
  real, dimension(F_n,F_nk) :: F_moins          !Coupled component state at previous step to be updated
  real, dimension(F_n,F_nk) :: F_plus           !Coupled component state at the end of current component step

  ! Flip time-plus into time-minus for next step
  F_moins = F_plus

  ! End of subprogram
  return
end subroutine ccpl_timeflip
