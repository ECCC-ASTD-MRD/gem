
module phy_restart_mod
  use clib_itf_mod, only: clib_toupper

#ifdef HAVE_NEMO
  use cpl_itf, only: cpl_restart
#endif

  implicit none
  private
  public :: phy_restart

  character(len=4), parameter, public :: PHY_RESTART_READ  = 'R'
  character(len=4), parameter, public :: PHY_RESTART_WRITE = 'W'

contains

  !/@*
  function phy_restart(F_WorR_S, F_spin_L) result(F_istat)
    implicit none
 !#TODO: F_spin_L never used
    character(len=*),intent(in) :: F_WorR_S
    logical,         intent(in) :: F_spin_L
    integer :: F_istat  !Return status (RMN_OK or RMN_ERR)

    !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2015
    !@revision
    !*@/

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

    integer :: istat
    character(len=32) :: WorR_S

    ! ------------------------------------------------------------------
    WorR_S = adjustl(F_WorR_S)
    istat = clib_toupper(WorR_S)

! physics has currently nothing to do for restart

    if (WorR_S == PHY_RESTART_READ) then
    endif

    if (WorR_S == PHY_RESTART_WRITE) then
    endif

    if (F_spin_L) then
       print *, 'F_spin_L ignored'
    endif

#ifdef HAVE_NEMO
! coupling may have something to do for restart
    call cpl_restart(WorR_S)
#endif

    F_istat = RMN_OK

    ! ------------------------------------------------------------------
    return
  end function phy_restart

end module phy_restart_mod
