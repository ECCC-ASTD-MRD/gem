      subroutine out_cfile()
      use out_mod
      use out3
      implicit none
#include <arch_specific.hf>
!
!
      integer, external :: fstfrm
      integer err
      real dummy
!
!----------------------------------------------------------------------
!
      call out_fstecr ( dummy,dummy,dummy,dummy,dummy,dummy,&
                         dummy,dummy,dummy,dummy,dummy,dummy,&
                         dummy,dummy,dummy, .true. )

      if ( (Out3_iome >= 0) .and. (Out_unf > 0) ) then
         err = fstfrm(Out_unf)
         call fclos(Out_unf)
         Out_unf = 0
      end if
!----------------------------------------------------------------------
      return
      end
