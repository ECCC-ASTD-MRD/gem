!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------

!/@*
subroutine drv_print_banner(F_begin_L)
   use, intrinsic :: iso_fortran_env, only: INT64
   implicit none
   !@objective Print exdb banner
   !@arguments
   logical,intent(in) :: F_begin_L
   !@author Stephane Chamberland, 2012-02
!*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   character(len=512),save :: version_number_S = '',version_title_S = ''

   integer ::istat,msgUnit
   character(len=512) :: model_name_S,dstp_S,arch_S,compil_S,user_S
   logical :: is_official_L
   !---------------------------------------------------------------------
   msgUnit = msg_getUnit(MSG_INFO)
   if (msgUnit > 0) then
      if (version_title_S == '') then
         call atm_model_getversion2(model_name_S, version_number_S, dstp_S, &
              arch_S, compil_S, user_S, is_official_L)
         if (is_official_L) then
            version_title_S = trim(model_name_S)//' --- Release of: '//trim(dstp_S)
         else
            version_title_S = trim(model_name_S)//' --- '//trim(user_S)//' Build: '//trim(dstp_S)
         endif

      endif
      if (F_begin_L) then
         istat = exdb(trim(version_title_S),trim(version_number_S),'NON')
      else
         istat = exfin(trim(Version_title_S),trim(Version_number_S),'OK')
      endif
   endif
   !---------------------------------------------------------------------
   return
end subroutine drv_print_banner


