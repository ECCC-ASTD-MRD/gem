
! Declaration and initialization of indices for phybus


module phy_status
   public
   save

   integer, parameter :: PHY_ERROR = -1
   integer, parameter :: PHY_NONE  = 0
   integer, parameter :: PHY_OK    = 1
   integer, parameter :: PHY_CTRL_NML_OK   = 1001
   integer, parameter :: PHY_CTRL_INI_OK   = 1003

   logical :: phy_error_L = .false.

   integer :: phy_init_ctrl = PHY_ERROR

contains

   !# This s/r is needed for the compiler to produce a .o of this file (make requirement)
   subroutine phy_status_dummy()
      print *,'dummy'
   end subroutine phy_status_dummy

end module phy_status
