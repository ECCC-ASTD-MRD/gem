
module phy_terminate_mod

  private
  public :: phy_terminate

contains

  !/@*
  function phy_terminate() result(F_istat)

#ifdef HAVE_NEMO
      use cpl_itf, only: cpl_terminate
#endif

      use series_mod, only: series_terminate
      use phy_status, only: PHY_NONE, PHY_CTRL_INI_OK, phy_error_L, phy_init_ctrl
      implicit none

    integer :: F_istat  !Return status (RMN_OK or RMN_ERR)

    !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014
    !@revision

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

    integer :: istat
    ! ---------------------------------------------------------------------
    F_istat = RMN_ERR
    if (phy_init_ctrl == PHY_NONE) then
       F_istat = PHY_NONE
       return
    else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
       call msg(MSG_ERROR,'(phy_terminate) Physics not properly initialized.')
       return
    endif

    ! Set error return value
    F_istat = RMN_ERR

    ! Shut down time series
    istat = series_terminate()

#ifdef HAVE_NEMO
    ! Shut down coupler
    call cpl_terminate (.true.)
#endif

    if (phy_error_L) return

    ! Successful completion
    F_istat = RMN_OK

    ! ---------------------------------------------------------------------
    return
  end function phy_terminate

end module phy_terminate_mod
