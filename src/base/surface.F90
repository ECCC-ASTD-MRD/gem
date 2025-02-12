
module surface
   implicit none
   private
   public :: surface1

contains

   subroutine surface1(pvars, delt, kount, ni, nk, trnch)
      ! Call API for surface schemes
      use sfc_main, only: sfc_main2
      use phy_options, only: fluvert, timings_L
      use phy_status, only: phy_error_L
      use phybusidx
      use phymem, only: phyvar
      use tdpack_const, only: CAPPA
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      ! Input arguments
      type(phyvar), pointer, contiguous :: pvars(:)  !all phy vars (meta + slab data)
      integer, intent(in) :: trnch                     !Slice number
      integer, intent(in) :: kount                     !Time step number
      real, intent(in) :: delt                         !Time step (sec)
      integer, intent(in) :: ni                        !Number of columns
      integer, intent(in) :: nk                        !Number of levels

      ! Internal variables
      integer :: istat, i, nkm1
      real :: sc
      real, dimension(:), pointer, contiguous :: zthetaa, zthetaap
      real, dimension(:,:), pointer, contiguous :: zsigt, ztmoins, &
           ztplus
#include "phymkptr.hf"

      ! Do not run surface schemes unless the PBL is active
      if (fluvert == 'NIL') return

      ! Associate pointers to bus entries
      nkm1 = nk-1
      MKPTR2Dm1(zsigt, sigt, pvars)
      MKPTR2Dm1(ztmoins, tmoins, pvars)
      MKPTR2Dm1(ztplus, tplus, pvars)
      MKPTR1D(zthetaa, thetaa, pvars)
      MKPTR1D(zthetaap, thetaap, pvars)     

      ! Compute first-level potential temperatures
      do i=1,ni
         sc = zsigt(i,nkm1)**(-CAPPA)
         zthetaa(i) = sc*ztmoins(i,nkm1)
         zthetaap(i) = sc*ztplus(i,nkm1)
      enddo

      ! Call main surface driver
      if (timings_L) call timing_start_omp(425, 'surface', 46)
      istat = sfc_main2(pvars, trnch, kount, delt, ni, nk)
      if (timings_L) call timing_stop_omp(425)
      if (phy_error_L .or. .not.RMN_IS_OK(istat)) then
         call physeterror('surface', 'Problem in sfc_main')
         return
      endif

   end subroutine surface1

end module surface
