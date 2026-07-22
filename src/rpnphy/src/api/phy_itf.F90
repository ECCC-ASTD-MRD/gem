
module phy_itf

  ! Subprogram elements of the API
  use phy_typedef
  use phy_nml_mod
  use phy_init_mod
  use phy_input
  use phy_step_mod
  use phy_get_mod
  use phy_getmeta_mod
  use phy_put_mod
  use phy_putmeta_mod
  use phy_terminate_mod
  use phy_restart_mod
  use phy_snapshot_mod

  implicit none

  !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014
  !@revision

  !@objective 
  ! Fill the external interface (API) for the physics package

  ! Allow access to the full API
  public
  private :: PATHLENGTH, NPATH_DEFAULT, BPATH_DEFAULT

  ! Externally-accessible constants and variables
  integer, parameter :: PHY_COMPATIBILITY_LVL = 21

end module phy_itf
