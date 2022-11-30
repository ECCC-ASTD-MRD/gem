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
!-------------------------------------- LICENCE END ---------------------------

!/@*
subroutine litozon(F_file_S, F_myproc)
   use mu_jdate_mod, only: jdate_month
   use phy_options
   implicit none
!!!#include <arch_specific.hf>

   character(len=*), intent(in) :: F_file_S
   integer, intent(in) :: F_myproc

   !@Author B. Bilodeau (April 1994) - From lirozon
   !@Revision
   !
   ! 001      M. Desgagne (Oct 98) - call back to rdradf_d (from dynamics)
   ! 002      M. Desgagne (Mar 08) - optional ozone file content
   !@Object   read in the ozone climatology table
   !@Arguments
   !          - Input -
   ! F_file_S  full path of the radiation table file
   !@Notes
   !     1 -  The ozone climatology file provided by J.P.Blanchet
   !          is based on Kita and Sumi 1986.
   !     2 -  Interpolation to the day of the month is done and all
   !          monthly climatological values are supposed to be valid
   !          on the 15th. No interpolation is needed (done) when
   !          the input jour is 15
   !*@/

#include <rmn/msg.h>
#include "ozopnt.cdk"

   external :: rd_ozone
   !-----------------------------------------------------------------
   call msg(MSG_INFO,'litozon')

   call phyrdfile(F_file_S, rd_ozone, 'OZONE', F_myproc)

   fozon = 1
   clat = fozon + nlacl*npcl
   pref = clat + nlacl
   !-----------------------------------------------------------------
   return
end subroutine litozon

