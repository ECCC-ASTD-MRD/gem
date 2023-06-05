      subroutine gem_error (F_errorCode, F_FromSubName, F_Message)
      use iso_c_binding
      use app
      implicit none

      integer :: F_errorCode
      character(len=*) :: F_FromSubName
      character(len=*) :: F_Message

   include "rpn_comm.inc"
      integer :: errcode, err
!
!     ---------------------------------------------------------------
!
      call rpn_comm_allreduce (F_errorCode, errcode,1,RPN_COMM_INTEGER,&
                               "MPI_MIN",RPN_COMM_MULTIGRID,err)

      if (errcode < 0) then
         call app_log(APP_FATAL,F_FromSubName//': '//F_Message)
         call rpn_comm_FINALIZE(err)
         app_status=app_end(errcode)
         stop app_status
      end if

   !---------------------------------------------------------------------
   return
end subroutine gem_error
