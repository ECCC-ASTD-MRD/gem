

function phyput_input_param() result(istat)
   use, intrinsic :: iso_fortran_env, only: INT64
   use wb_itf_mod, only: WB_OK, WB_MSG_INFO, WB_REWRITE_MANY, WB_IS_LOCAL, wb_put, wb_verbosity
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !@objective initialization of the surface input parameters
   !@Author L. Spacek (Fall 2013)
   !@Revisions
   ! 001      L. Spacek (Fall 2013)   - initial revision
   !@return
   integer :: istat

#include <rmn/msg.h>

   include "phyinput.inc"
   include "physteps.cdk"

   integer :: options,iverb
   !-------------------------------------
   options = WB_REWRITE_MANY+WB_IS_LOCAL
   iverb = wb_verbosity(WB_MSG_INFO)
   istat = WB_OK
   istat = min(wb_put('phyinput/phyinread_max'   ,phyinread_max   ,options),istat)
   istat = min(wb_put('phyinput/phyinread_n'     ,phyinread_n     ,options),istat)
   istat = min(wb_put('phyinput/phyinread_jdateo',phyinread_jdateo ,options),istat)
   istat = min(wb_put('phyinput/phyinread_dt'    ,phyinread_dt    ,options),istat)
   istat = min(wb_put('phyinput/phyinread_step'  ,phyinread_step  ,options),istat)
   istat = min(wb_put('phyinput/phyinread_list_nk',phyinread_list_nk,options),istat)
   istat = min(wb_put('phyinput/phyinread_list_S',phyinread_list_S,options),istat)
   istat = min(wb_put('physteps/step_kount'      ,step_kount      ,options),istat)
   istat = min(wb_put('physteps/step_driver'     ,step_driver     ,options),istat)
   iverb = wb_verbosity(iverb)
   if (istat /= WB_OK) then
      istat = -1
      call msg(MSG_ERROR,'(phyput_input_param)')
   endif
   !-------------------------------------
   return
end function  phyput_input_param

