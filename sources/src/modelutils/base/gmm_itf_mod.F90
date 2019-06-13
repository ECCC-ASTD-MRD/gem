!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------
!/@*
module gmm_itf_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   implicit none
   public
   !@objective Module form of gmm include file
   !@author  Stephane Chamberland, 2010-06
!*@/

#include <mu_gmm.hf>

#undef GMM_IS_OK
#undef GMM_IS_ERROR

   contains

      logical function GMM_IS_OK(errcode)
         implicit none
         integer, intent(in) :: errcode
         GMM_IS_OK = (errcode >= 0)
      end function

      logical function GMM_IS_ERROR(errcode)
         implicit none
         integer, intent(in) :: errcode
         GMM_IS_ERROR = (errcode < 0)
      end function

end module gmm_itf_mod
