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

!@*
function sfc_init1() result(F_istat)
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>
   !@Object initialization of the surface parameters at the beginning
   !        of each execution of the model
   !@Arguments
   integer :: F_istat

   !@Author L. Spacek (Fall 2012)
   !@Revisions
   ! 001      L. Spacek (Fall 2012)   - initial revision
   !@Description
   !          sfcdebu does the following :
   !          1) it initializes a few constants necessary
   !             for the execution of the surface package.
   !          2) call iniptsurf to construct the surface bus
!*@/
#include <rmnlib_basics.hf>
   integer, external :: iniptsurf5
   !---------------------------------------------------------------------
   !#TODO: remove this useless step
   F_istat = iniptsurf5()
   !----------------------------------------------------------------------
end function sfc_init1
