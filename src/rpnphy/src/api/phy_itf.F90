!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer, 
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms 
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer 
!version 3 or (at your option) any later version that should be found at: 
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html 
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software; 
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec), 
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

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
  integer, parameter :: PHY_COMPATIBILITY_LVL = 20

end module phy_itf
