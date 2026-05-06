module cpl_restart_mod

   private
   public :: cpl_restart

contains

      subroutine cpl_restart (F_WorR_S)
      use cpl_mod
      implicit none
#include <arch_specific.hf>
      character(len=*),intent(in) :: F_WorR_S

      if ( trim(F_WorR_S) == 'R' ) cpl_rstn_S='THIS_IS_A_RESTART'

      return
      end subroutine cpl_restart

end module cpl_restart_mod
