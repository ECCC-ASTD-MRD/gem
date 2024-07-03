      subroutine out_cfile()
      use out_mod
      use out3
      implicit none
#include <arch_specific.hf>
!
!
      real dummy
!
!----------------------------------------------------------------------
!
      call out_fstecr ( dummy,dummy,dummy,dummy,dummy,dummy,&
                         dummy,dummy,dummy,dummy,dummy,dummy,&
                         dummy,dummy,dummy, .true. )

      if ( (Out3_iome >= 0) .and. Out_file%is_open() ) then
         success = Out_file%close()
      end if
!----------------------------------------------------------------------
      return
      end
