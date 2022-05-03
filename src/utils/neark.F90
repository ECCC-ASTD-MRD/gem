!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!**S/P  NEARK
!
function neark(sig,ps,delp,ni,nk,kindex) result(stat)
   implicit none
!!!#include <arch_specific.hf>

   ! Input arguments
   integer, intent(in) :: ni,nk                         !Slab dimensions
   real, dimension(ni,nk), intent(in) :: sig            !Sigma level values for profile
   real, dimension(ni), intent(in) :: ps                !Surface pressure (Pa)
   real, intent(in) :: delp                             !Pressure depth for level search (Pa)

   ! Output arguments
   integer, dimension(ni), intent(out) :: kindex        !Level index closes to "delp" above the surface
   integer :: stat                                      !Return status (0 success; -1 error)

   ! Internal variables
   integer :: i,k
   real, dimension(ni) :: sigtarget

   ! Set error return values
   kindex = nk
   stat = -1

   ! Compute nearest level index
   sigtarget = 1.-delp/ps
   do i=1,ni
      k = 2      
      do while (kindex(i) == nk .and. k <= nk)
         if (sig(i,k) > sigtarget(i) .and. sig(i,k-1) <= sigtarget(i)) kindex(i) = k
         k = k+1
      enddo
   enddo

   ! End of subprogram
   stat = 0
   return
end function neark
