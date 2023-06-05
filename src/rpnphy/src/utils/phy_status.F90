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
