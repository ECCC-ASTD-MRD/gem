   
module phy_status
   public
   save

   integer, parameter :: PHY_ERROR = -1
   integer, parameter :: PHY_NONE  = 0
   integer, parameter :: PHY_OK    = 1
   integer, parameter :: PHY_CTRL_NML_OK   = 1001
   integer, parameter :: PHY_CTRL_INI_OK   = 1003

   integer :: phy_init_ctrl = PHY_ERROR

   logical :: phy_error_L = .false.

contains

   !/@*
   subroutine physeterror(F_from_S, F_msg_S)
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
   
   !/@*
   function phyhaserror() result(F_stat)
      implicit none
      !@object return true is phy has errors
      !@return
      logical :: F_stat
      !@author 
      !*@/
      !----------------------------------------------------------------
      F_stat = phy_error_L
      !----------------------------------------------------------------
      return
   end function phyhaserror

end module phy_status
