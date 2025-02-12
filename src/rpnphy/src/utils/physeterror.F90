
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
