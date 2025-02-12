
module radiation
   implicit none
   private
   public :: radiation3

contains

   !/@*
   subroutine radiation3(pvars, kount, ni, nk, trnch)
      use iso_c_binding
      use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
      use debug_mod, only: init2nan
      use ccc1_cccmarad, only: cccmarad1
      use ccc2_cccmarad, only: ccc2_cccmarad2
      use phy_options
      use phy_status, only: phy_error_L
      use phybusidx
      use phymem, only: phyvar
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
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)

      integer, intent(in) :: ni, nk, kount, trnch
      type(phyvar), pointer, contiguous :: pvars(:)

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

      MKPTR1D(zpmoins, pmoins, pvars)

      MKPTR2D(ztmoins, tmoins, pvars)
      MKPTR2D(zhumoins, humoins, pvars)
      MKPTR2D(zsigw, sigw, pvars)

      call init2nan(cldfrac, liqwcin, icewcin, liqwp, icewp, trav2d)

      nkm1 = nk - 1

      IF_RAD: if (radia /= 'NIL') then

         select case (radia)
         case('CCCMARAD')

            call cccmarad1(pvars, &
                 ztmoins, zhumoins, &
                 zpmoins, zsigw, delt, kount, &
                 trnch, ni, nkm1, nk, &
                 liqwcin, icewcin, liqwp, icewp, cldfrac)
            if (phy_error_L) return

         case('CCCMARAD2')

            call ccc2_cccmarad2(pvars, &
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

      call radslop3(pvars, hz, julien, ni, trnch)

      if (timings_L) call timing_stop_omp(410)
      call msg_toall(MSG_DEBUG, 'radiation [END]')
      !----------------------------------------------------------------
      return
   end subroutine radiation3

end module radiation
