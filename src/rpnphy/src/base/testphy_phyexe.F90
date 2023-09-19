module testphy_phyexe
   private
   public :: testphy_phyexe1

contains

   subroutine testphy_phyexe1(pvars, kount, ni, nk, trnch)
      ! Test double for RPN Physics parameterizations: increment time-plus 
      ! temperatures throughout the column by 1K and put 100 into the
      ! boundary layer height.
      use phymem, only: phyvar, PHY_BUSIDXV, PHY_DBUSIDX, PHY_PBUSIDX, PHY_VBUSIDX
      use phybusidx
      implicit none

      ! Architecture-dependent content
!!!#include <arch_specific.hf>

      ! Input arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - input -
      ! trnch    slice number
      ! kount    timestep number
      ! ni       horizontal running length
      ! nk       vertical dimension

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, kount, ni, nk

      ! External declarations
#include "phymkptr.hf"

      ! Local variables
      integer :: i,k
      real, pointer, contiguous :: ztplus(:,:), zh(:)

      ! Associate pointers
      MKPTR2D(ztplus, tplus, pvars)
      MKPTR1D(zh, h, pvars)

      ! Increment temperatures
      do k=1,nk
         do i=1,ni
            ztplus(i,k) = ztplus(i,k) + 1.
         enddo
      enddo

      ! Specify boundary layer height
      zh(:) = 100.

      if (ni<0) print *,'called testphy_phyexe with:',trnch,kount,ni,nk
      ! End of test double
      return
   end subroutine testphy_phyexe1

end module testphy_phyexe
