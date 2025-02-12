
!/@*
subroutine serdbu()
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@object: Initializes series/profiles buffers to zero
   !*@/
   include "series.cdk"
   !---------------------------------------------------------------
   surface(:,2) = '        '
   profils(:,2) = '        '
   sers    = 0.
   serp    = 0.
   !---------------------------------------------------------------
   return
end subroutine serdbu
