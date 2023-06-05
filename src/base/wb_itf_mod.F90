!--------------------------------------------------------------------------
! This is free software, you can use/redistribute/modify it under the terms of
! the EC-RPN License v2 or any later version found (if not provided) at:
! - http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
! - EC-RPN License, 2121 TransCanada, suite 500, Dorval (Qc), CANADA, H9P 1J3
! - service.rpn@ec.gc.ca
! It is distributed WITHOUT ANY WARRANTY of FITNESS FOR ANY PARTICULAR PURPOSE.
!--------------------------------------------------------------------------
!/@*
module wb_itf_mod
   implicit none
   public
   !@objective Module form of whiteboard include file
   !@author  Stephane Chamberland, 2010-06
!*@/

#include <rmn/WhiteBoard.hf>

   integer,parameter :: wb_OK = WB_OK
   integer,parameter :: wb_ERROR = WB_ERROR

   integer,parameter :: wb_MAXNAMELENGTH = WB_MAXNAMELENGTH
   integer,parameter :: wb_MAXSTRINGLENGTH = WB_MAXSTRINGLENGTH

   integer,parameter :: wb_DICT_FILE = WB_STRICT_DICTIONARY
   integer,parameter :: wb_CONF_FILE = WB_FORBID_DEFINE
   integer,parameter :: wb_STRICT_DICTIONARY = WB_STRICT_DICTIONARY
   integer,parameter :: wb_FORBID_DEFINE = WB_FORBID_DEFINE
   integer,parameter :: wb_REWRITE_AT_RESTART = WB_REWRITE_AT_RESTART
   integer,parameter :: wb_REWRITE_MANY = WB_REWRITE_MANY
   integer,parameter :: wb_REWRITE_NONE = WB_REWRITE_NONE
   integer,parameter :: wb_IS_LOCAL = WB_IS_LOCAL
   integer,parameter :: wb_INITIALIZED = WB_INITIALIZED

   integer,parameter :: wb_MSG_DEBUG = WB_MSG_DEBUG
   integer,parameter :: wb_MSG_INFO = WB_MSG_INFO
   integer,parameter :: wb_MSG_WARN = WB_MSG_WARN
   integer,parameter :: wb_MSG_ERROR = WB_MSG_ERROR
   integer,parameter :: wb_MSG_SEVERE = WB_MSG_SEVERE
   integer,parameter :: wb_MSG_FATAL = WB_MSG_FATAL

   integer,parameter :: wb_FORTRAN_REAL = WB_FORTRAN_REAL
   integer,parameter :: wb_FORTRAN_INT = WB_FORTRAN_INT
   integer,parameter :: wb_FORTRAN_CHAR = WB_FORTRAN_CHAR
   integer,parameter :: wb_FORTRAN_BOOL = WB_FORTRAN_BOOL

   !Note: Other macro consts avail in <WhiteBoard.h>

#undef WB_IS_OK
#undef WB_IS_ERROR

   contains

      logical function WB_IS_OK(errcode)
         implicit none
         integer, intent(in) :: errcode
         WB_IS_OK = (errcode >= 0)
      end function

      logical function WB_IS_ERROR(errcode)
         implicit none
         integer, intent(in) :: errcode
         WB_IS_ERROR = (errcode < 0)
      end function

end module wb_itf_mod
