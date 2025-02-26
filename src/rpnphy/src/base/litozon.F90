
!/@*
subroutine litozon(F_file_S, F_myproc)
   use mu_jdate_mod, only: jdate_month
   use phy_options
   use phyrdfile, only: phyrdfile1, READOZO
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

   call phyrdfile1(F_file_S, READOZO, 'OZONE', F_myproc)

   fozon = 1
   clat = fozon + nlacl*npcl
   pref = clat + nlacl
   !-----------------------------------------------------------------
   return
end subroutine litozon

