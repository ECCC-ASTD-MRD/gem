!-------------------------------------- LICENCE BEGIN --------------------------
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

module metox
   implicit none
   private
   public :: metox3

contains

   !/@*
   subroutine metox3(dbus, vbus, fbus, ni, nk)
      use debug_mod, only: init2nan
      use phy_options
      use phybus
      use tendency, only: apply_tendencies
      implicit none
!!!#include <arch_specific.hf>
      !@object Methane oxidation
      ! Produces the specific humidity tendency due to methane
      ! oxidation (based on ECMWF scheme)
      !@arguments
      ! ni       horizonal index
      ! nk       vertical  index
      ! vbus     volatile bus
      ! fbus     permanent bus
      ! dbus     dynamics bus

      integer, intent(in) :: ni, nk
      real, dimension(:), pointer, contiguous :: dbus, fbus, vbus

      !@author M. Charron (RPN): November 2005
      !@revisions
      ! * L. Spacek (oct 2011) - Complete remake
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"

      real, parameter :: QQ = 4.25e-6 
      real, parameter :: ALPHA1 = 0.5431969

      real, pointer, dimension(:), contiguous :: psp, ztdmask
      real, pointer, dimension(:,:), contiguous :: sigma, oxme, qqp, zhuplus
      real, dimension(ni,nk) :: press, kmetox
      !----------------------------------------------------------------

      if (.not.lmetox) return
      call msg_toall(MSG_DEBUG, 'metox [BEGIN]')
      if (timings_L) call timing_start_omp(415, 'metox', 46)

      MKPTR1D(psp, pmoins, fbus)
      MKPTR1D(ztdmask, tdmask, fbus)
      MKPTR2D(zhuplus, huplus, dbus)
      MKPTR2D(sigma, sigw, dbus)
      MKPTR2D(qqp, huplus, dbus)
      MKPTR2D(oxme, qmetox, vbus)

      call init2nan(press, kmetox)

      press(:,:) = spread(psp(:), dim=2, ncopies=nk)
      press(:,:) = sigma(:,1:nk)*press(:,:)

      where(press >= 10000.)
         kmetox = 0.
      elsewhere (press > 50..and. press < 10000.)
         kmetox = 1.157407e-7/(1.+ALPHA1*(alog(press/50.))**4/(alog(10000./press)))
      elsewhere
         kmetox = 1.157407e-7
      endwhere

      oxme = kmetox*(QQ-qqp)

      call apply_tendencies(zhuplus, oxme, ztdmask, ni, nk, nk-1)

      if (timings_L) call timing_stop_omp(415)
      call msg_toall(MSG_DEBUG, 'metox [END]')
      !----------------------------------------------------------------
      return
   end subroutine metox3

end module metox
