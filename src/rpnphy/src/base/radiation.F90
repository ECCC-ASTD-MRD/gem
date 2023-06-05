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

module radiation
   implicit none
   private
   public :: radiation3

contains

   !/@*
   subroutine radiation3(dbus, fbus, vbus, &
        ni, nk, kount, trnch)
      use iso_c_binding
      use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
      use debug_mod, only: init2nan
      use ccc1_cccmarad, only: cccmarad1
      use ccc2_cccmarad, only: ccc2_cccmarad2
      use phy_options
      use phy_status, only: phy_error_L
      use phybus
      use radslop, only: radslop3
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@Object Interface to radiation
      !@Arguments
      !          - Input -
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! kount    timestep number
      ! trnch    slice number
      !          - Input/Output -
      ! dbus     dynamics input field
      ! fbus     historic variables for the physics
      ! vbus     physics tendencies and other output fields from the physics

      integer, intent(in) :: ni, nk, kount, trnch
      real, pointer, contiguous :: dbus(:), fbus(:), vbus(:)

      !@Author L.Spacek, November 2011
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"

      integer :: yy, mo, dd, hh, mn, ss, nkm1
      real :: hz0, hz, julien

      real, dimension(ni,nk) :: cldfrac, liqwcin, icewcin, liqwp, icewp, trav2d
      real, pointer, dimension(:), contiguous :: zpmoins
      real, pointer, dimension(:,:), contiguous :: ztmoins, zhumoins, zsigw
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'radiation [BEGIN]')
      if (timings_L) call timing_start_omp(410, 'radiation', 46)

      MKPTR1D(zpmoins, pmoins, fbus)

      MKPTR2D(ztmoins, tmoins, dbus)
      MKPTR2D(zhumoins, humoins, dbus)
      MKPTR2D(zsigw, sigw, dbus)

      call init2nan(cldfrac, liqwcin, icewcin, liqwp, icewp, trav2d)

      nkm1 = nk - 1

      IF_RAD: if (radia /= 'NIL') then

         select case (radia)
         case('CCCMARAD')

            call cccmarad1(dbus, fbus, vbus, &
                 ztmoins, zhumoins, &
                 zpmoins, zsigw, delt, kount, &
                 trnch, ni, nkm1, nk, &
                 liqwcin, icewcin, liqwp, icewp, cldfrac)
            if (phy_error_L) return

         case('CCCMARAD2')

            call ccc2_cccmarad2(dbus, fbus, vbus, &
                 ztmoins, zhumoins, &
                 zpmoins, zsigw, delt, kount, &
                 trnch, ni, nkm1, nk, &
                 liqwcin, icewcin, liqwp, icewp, cldfrac)
            if (phy_error_L) return

         end select

      endif IF_RAD

      call mu_js2ymdhms(jdateo, yy, mo, dd, hh, mn, ss)
      hz0 = hh + float(mn)/60. + float(ss)/3600.
      hz = amod(hz0 + (float(kount)*delt)/3600., 24.)
      julien = real(jdate_day_of_year(jdateo + kount*int(delt) + MU_JDATE_HALFDAY))

      call radslop3(fbus, vbus, ni, hz, julien, trnch)

      if (timings_L) call timing_stop_omp(410)
      call msg_toall(MSG_DEBUG, 'radiation [END]')
      !----------------------------------------------------------------
      return
   end subroutine radiation3

end module radiation
