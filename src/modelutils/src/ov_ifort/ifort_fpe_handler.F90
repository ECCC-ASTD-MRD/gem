
#include <rmn/msg.h>

subroutine fpe_setup()

#ifdef __INTEL_COMPILER

   use ifport
   implicit none
   interface
      subroutine fpe_handler(signo, siginfo)
         integer(4), intent(in) :: signo, siginfo
      end subroutine fpe_handler
   end interface
   integer :: ir
   ir = ieee_handler('set', 'invalid',  fpe_handler)
   ir = ieee_handler('set', 'division', fpe_handler)
   ir = ieee_handler('set', 'overflow', fpe_handler)
!!$   ir = ieee_handler('set', 'undeflow', fpe_handler)
!!$   ir = ieee_handler('set', 'inexact',  fpe_handler)

#else

   call msg(MSG_WARNING,'(fpe_setup) Not yet available')

#endif

   return
end subroutine fpe_setup


subroutine fpe_handler(sig, code)
#ifdef __INTEL_COMPILER
   use ifcore
   use ifport
#endif
   implicit none
   integer :: sig, code
   character(len=128) :: msg_s

#ifdef __INTEL_COMPILER

   if (code == FPE$INVALID .or. code == FPE$ZERODIVIDE .or. code == FPE$OVERFLOW) then
      call msg_buffer_flush()
   endif
   write(msg_s, '(a,1x,i0,1x,a,1x,i0,1x,a)') 'Traceback: Application SIGFPE error! sig=', sig,', code=', code, ':'
   if (code == FPE$INVALID) msg_s = trim(msg_s)//'INVALID'
   if (code == FPE$ZERODIVIDE) msg_s = trim(msg_s)//'ZERODIVIDE'
   if (code == FPE$OVERFLOW) msg_s = trim(msg_s)//'OVERFLOW'
   if (code == FPE$UNDERFLOW) msg_s = trim(msg_s)//'UNDERFLOW'
   if (code == FPE$INEXACT) msg_s = trim(msg_s)//'INEXACT'
   call tracebackqq(msg_s, user_exit_code=-1) !#, eptr=eptrs)

#else

   call msg(MSG_WARNING,'(fpe_handler) Not yet available')
   write(msg_s, '(a,1x,i0,1x,a,1x,i0)') 'called with: sig=', sig,', code=', code
   call msg(MSG_WARNING,'(fpe_handler) '//msg_s)
   call flush(6)

#endif

   call abort()
end subroutine fpe_handler
