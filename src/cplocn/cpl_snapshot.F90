module cpl_snapshot_mod

   private
   public :: cpl_snapshot

contains

      subroutine cpl_snapshot (F_mode)
      implicit none
#include <arch_specific.hf>
      character(len=*), intent(in) :: F_mode

      include "cpl.cdk"

      if ( F_mode == 'W' ) cpl_dgflt_H = .true.
      if ( F_mode == 'R' ) cpl_dgflt_H = .false.

      return
      end subroutine cpl_snapshot

end module cpl_snapshot_mod
