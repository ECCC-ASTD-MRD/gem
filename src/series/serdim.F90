
!/@*
function serdim(nstn,n,nk) result(F_dim)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@object Determine vector dimensions for time-series
   !@params
   ! nstn     number of stations
   ! n        number of surface or profile variables requested
   ! nk       vertical dimension
   integer,intent(in) :: nstn,n,nk
   !@return
   ! F_dim    dimension to be allocated by the dynamics
   integer :: F_dim
   !@author m. desgagne (mar 99)
   !*@/
   include "series.cdk"
   !---------------------------------------------------------------
   F_dim = 0
   if (nstn <= 0) return
   mxstt = nstn
   if (nk == 1) then
      mxsrf = n
      F_dim = mxstt * mxsrf
   else
      mxprf = n
      mxnvo = nk
      F_dim = mxstt * mxprf * mxnvo
   endif
   !---------------------------------------------------------------
   return
end function serdim
