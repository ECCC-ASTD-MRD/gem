!-------------------------------------- LICENCE BEGIN -------------------------
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

function sfc_get_input_param()  result(istat)
   use, intrinsic :: iso_fortran_env, only: INT64
   use wb_itf_mod, only: WB_OK, WB_MSG_INFO, wb_verbosity, wb_get
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@Object initialization of the surface input parameters
   !@returns
   integer :: istat
   !@Author L. Spacek (Fall 2013)
   !@Revisions
   !*
#include <rmn/msg.h>

   include "sfcinput.cdk"

   integer :: nv,iverb
   !-------------------------------------
   iverb = wb_verbosity(WB_MSG_INFO)
   istat = WB_OK
   istat = min(wb_get('phyinput/phyinread_max'   ,phyinread_dim   ),istat)
   istat = min(wb_get('phyinput/phyinread_n'     ,phyinread_n     ),istat)
   istat = min(wb_get('phyinput/phyinread_jdateo',phyinread_jdateo ),istat)
   istat = min(wb_get('phyinput/phyinread_dt'    ,phyinread_dt    ),istat)
   istat = min(wb_get('phyinput/phyinread_step'  ,phyinread_step  ),istat)
   istat = min(wb_get('phyinput/phyinread_list_nk',phyinread_list_nk,nv),istat)
   istat = min(wb_get('phyinput/phyinread_list_S',phyinread_list_S,nv),istat)
!!$   istat = min(wb_get('physteps/step_kount'      ,step_kount      ),istat)
!!$   istat = min(wb_get('physteps/step_driver'     ,step_driver     ),istat)
   iverb = wb_verbosity(iverb)
   if (istat /= WB_OK) then
      call msg(MSG_ERROR,'(sfc_get_input_param)')
   endif
   !-------------------------------------
   return
end function sfc_get_input_param


