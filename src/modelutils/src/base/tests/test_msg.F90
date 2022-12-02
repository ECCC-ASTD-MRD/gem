!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 
#include "rmn/msg.h"

!/@*
subroutine test_msg()
   implicit none
#include <utils4tests.hf>
   !@author  Stephane Chamberland, 2009-11
!*@/
   integer :: iErr,mype,i,iUnit
   integer, external :: msg_getUnit
   !---------------------------------------------------------------------
   !TODO: find a way to check for not printed msg that should be printed
   mype = utils4test_toy_mpi_init('test_msg')

   call msg(MSG_DEBUG,TEST_FAILED_MSG)
   call msg(MSG_INFO,TEST_FAILED_MSG)
   call msg(MSG_WARNING,TEST_PASSED_MSG)
   call msg(MSG_ERROR,TEST_PASSED_MSG)
   call msg(MSG_CRITICAL,TEST_PASSED_MSG)

   call msg_set_minMessageLevel(MSG_INFO)
   call msg(MSG_DEBUG,TEST_FAILED_MSG)
   call msg(MSG_INFO,TEST_PASSED_MSG)
   call msg(MSG_WARNING,TEST_PASSED_MSG)
   call msg(MSG_ERROR,TEST_PASSED_MSG)
   call msg(MSG_CRITICAL,TEST_PASSED_MSG)

   call msg_set_minMessageLevel(MSG_CRITICAL)
   call msg(MSG_DEBUG,TEST_FAILED_MSG)
   call msg(MSG_INFO,TEST_FAILED_MSG)
   call msg(MSG_WARNING,TEST_FAILED_MSG)
   call msg(MSG_ERROR,TEST_FAILED_MSG)
   call msg(MSG_CRITICAL,TEST_PASSED_MSG)

   call msg_set_p0only(mype)

   call msg_set_minMessageLevel(MSG_INFO)
   if (mype == 0) then
      call msg(MSG_DEBUG,TEST_FAILED_MSG)
      call msg(MSG_INFO,TEST_PASSED_MSG)
      call msg(MSG_WARNING,TEST_PASSED_MSG)
      call msg(MSG_ERROR,TEST_PASSED_MSG)
      call msg(MSG_CRITICAL,TEST_PASSED_MSG)
   else
      do i=1,MSG_NBLEVELS
         call msg(i,TEST_FAILED_MSG)
      enddo
   endif

   call msg_set_minMessageLevel(MSG_CRITICAL)
   if (mype == 0) then
      call msg(MSG_DEBUG,TEST_FAILED_MSG)
      call msg(MSG_INFO,TEST_FAILED_MSG)
      call msg(MSG_WARNING,TEST_FAILED_MSG)
      call msg(MSG_ERROR,TEST_FAILED_MSG)
      call msg(MSG_CRITICAL,TEST_PASSED_MSG)
   else
      do i=1,MSG_NBLEVELS
         call msg(i,TEST_FAILED_MSG)
      enddo
   endif

   call msg_set_redirect2fileUnit(MSG_STDERR,MSG_STDOUT)
   call msg_set_minMessageLevel(MSG_DEBUG)
   if (mype == 0) then
      call msg(MSG_DEBUG,TEST_PASSED_MSG)
      call msg(MSG_INFO,TEST_PASSED_MSG)
      call msg(MSG_WARNING,TEST_PASSED_MSG)
      call msg(MSG_ERROR,TEST_PASSED_MSG)
      call msg(MSG_CRITICAL,TEST_PASSED_MSG)
   else
      do i=1,MSG_NBLEVELS
         call msg(i,TEST_FAILED_MSG)
      enddo
   endif

   call msg_set_minMessageLevel(MSG_WARNING)
   iUnit = msg_getUnit(MSG_INFO)
   call assertFailed(iUnit,'should not print msg below min msgLevel')

   call msg_set_minMessageLevel(MSG_WARNING)
   iUnit = msg_getUnit(MSG_WARNING)
   if (mype == 0) then
      call assertOK(iUnit,'should print msg for proc==0')
   else
      call assertFailed(iUnit,'should not print msg for proc!=0')
   endif

   call utils4test_toy_mpi_end()
   !---------------------------------------------------------------------
   return
end subroutine test_msg


