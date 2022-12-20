!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------


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

