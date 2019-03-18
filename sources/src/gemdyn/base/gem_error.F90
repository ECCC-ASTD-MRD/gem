      subroutine gem_error (F_errorCode, F_FromSubName, F_Message)
use iso_c_binding
      use lun
      implicit none
#include <arch_specific.hf>
      integer :: F_errorCode
      character(len=*) :: F_FromSubName
      character(len=*) :: F_Message
      !@author  Michel Desgagne

   include "rpn_comm.inc"
      integer :: errcode, err
!
!     ---------------------------------------------------------------
!
      call rpn_comm_allreduce (F_errorCode, errcode,1,RPN_COMM_INTEGER,&
                               "MPI_MIN",RPN_COMM_MULTIGRID,err)

      if (errcode < 0) then
         if (Lun_out > 0) write(Lun_out,2000) F_FromSubName, F_Message
         call rpn_comm_FINALIZE(err)
         stop
      end if

2000  format (/60('#')/2x,a': ',a/2x,'ABORT'/60('#')/)
   !---------------------------------------------------------------------
   return
end subroutine gem_error
