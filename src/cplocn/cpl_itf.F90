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

module cpl_itf

  ! Subprogram elements of the API
  use cpl_nml_mod
  use cpl_init_mod
  use cpl_step_mod
  use cpl_update_mod
  use cpl_restart_mod
  use cpl_terminate_mod
  use cpl_snapshot_mod

  implicit none

  !@authors Desgagne, Roy -- Spring 2015
  !@revision

  !@objective 
  ! Fill the external interface (API) for the coupler package

  ! Allow access to the full API
  public

  ! Externally-accessible constants and variables
  integer, parameter :: CPL_COMPATIBILITY_LVL = 1

end module cpl_itf
