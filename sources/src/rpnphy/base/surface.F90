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

module surface
   implicit none
   private
   public :: surface1

contains

   subroutine surface1(trnch, kount, delt, ni, nk)
      ! Call API for surface schemes
      use sfc_main, only: sfc_main2
      use phy_options, only: fluvert, timings_L
      use phy_status, only: phy_error_L
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      ! Input arguments
      integer, intent(in) :: trnch                      !Slice number
      integer, intent(in) :: kount                      !Time step number
      real, intent(in) :: delt                          !Time step (sec)
      integer, intent(in) :: ni                         !Number of columns
      integer, intent(in) :: nk                         !Number of levels

      ! Internal variables
      integer :: istat
 
      ! Do not run surface schemes unless the PBL is active
      if (fluvert == 'NIL') return
 
      ! Call main surface driver
      if (timings_L) call timing_start_omp(425, 'surface', 46)
      istat = sfc_main2(trnch, kount, delt, ni, nk)
      if (timings_L) call timing_stop_omp(425)
      if (phy_error_L .or. .not.RMN_IS_OK(istat)) then
         call physeterror('surface', 'Problem in sfc_main')
         return
      endif

   end subroutine surface1

end module surface
