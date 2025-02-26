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
