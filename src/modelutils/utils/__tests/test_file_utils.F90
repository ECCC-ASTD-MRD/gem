!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 
!/@*
subroutine test_file_utils()
   implicit none
#include <utils4tests.hf>
   !@author  Stephane Chamberland, 2009-11
!*@/
   integer :: fileUnit,iErr,mype
   integer, external :: file_open_existing
   !---------------------------------------------------------------------
   mype = utils4test_toy_mpi_init('test_file_utils')

   fileUnit = file_open_existing('__DOES_NOT_EXIST__','FTN')
   call assertFailed(fileUnit,'should not open a non existing file')

   fileUnit = file_open_existing('../tests','FTN')
   call assertFailed(fileUnit,'should not open a dir')

   fileUnit = file_open_existing('test_file_utils.ftn90','FTN')
   call assertOK(fileUnit,'should open an existing file')

   !TODO: check for not existing but not readable file

   call utils4test_toy_mpi_end()
   !---------------------------------------------------------------------
   return
end subroutine test_file_utils


