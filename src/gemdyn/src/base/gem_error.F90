      subroutine gem_error (F_errorCode, F_FromSubName, F_Message)
      use iso_c_binding
      use app
      use inp_mod
      implicit none

      integer :: F_errorCode
      character(len=*) :: F_FromSubName
      character(len=*) :: F_Message

   include 'mpif.h'
   include "rpn_comm.inc"
      integer :: errcode, err
!
!     ---------------------------------------------------------------
!
      call rpn_comm_allreduce (F_errorCode, errcode,1,RPN_COMM_INTEGER,&
                               "MPI_MIN",RPN_COMM_MULTIGRID,err)

      if (errcode < 0) then
         call app_log(APP_FATAL,F_FromSubName//': '//F_Message)
         app_status=app_end(errcode)
         if (associated(Inp_recv)) call MPI_Win_free  (Inp_window,err)
         call rpn_comm_FINALIZE(err)
         stop! app_status
      end if

   !---------------------------------------------------------------------
   return
end subroutine gem_error
