!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!-------------------------------------------------------------------------- 
!/@*
module clib_itf_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   implicit none
   public
   !@objective Module form of clib_interfaces include file
   !@author  Stephane Chamberland, 2010-06
!*@/

#include <clib_interface_mu.hf>

   integer,parameter,public :: CLIB_ok = CLIB_OK
   integer,parameter,public :: CLIB_error = CLIB_ERROR

   contains

#undef CLIB_IS_OK
#undef CLIB_IS_ERROR

      logical function CLIB_IS_OK(errcode)
         implicit none
         integer, intent(in) :: errcode
         CLIB_IS_OK = (errcode >= 0)
      end function

      logical function CLIB_IS_ERROR(errcode)
         implicit none
         integer, intent(in) :: errcode
         CLIB_IS_ERROR = (errcode < 0)
      end function

end module clib_itf_mod
