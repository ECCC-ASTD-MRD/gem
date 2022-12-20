!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

!/@*
subroutine physeterror(F_from_S, F_msg_S)
   use phy_status, only: phy_error_L
   implicit none
   !@object Set physics error flag and print message
   !@arguments
   character(len=*), intent(in) :: &
        F_from_S, &  !# caller name
        F_msg_S      !# error message
   !@author 
   !*@/
#include <rmn/msg.h>
   !----------------------------------------------------------------
   call msg_toall(MSG_ERROR, '('//trim(F_from_S)//') '//trim(F_msg_S))
   phy_error_L = .true.
   !----------------------------------------------------------------
   return
end subroutine physeterror
