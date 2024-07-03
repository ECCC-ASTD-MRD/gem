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

module extdiag
   implicit none
   private
   public :: extdiag3

contains
   !/@*
   subroutine extdiag3(pvars, ni, nk, trnch)
      use tdpack_const, only: CAPPA
      use series_mod, only: series_xst, series_isstep, series_isvar
      use phy_options
      use phybusidx
      use phymem, only: phyvar, nphyvars
      implicit none
!!!#include <arch_specific.hf>
      !@Object calculate averages and accumulators of tendencies and diagnostics
      !@Arguments
      !          - input/output -
      ! pvars    list of all phy vars (meta + slab data)
      !          - input -
      ! trnch    slice number
      ! ni       horizontal running length
      ! nk       vertical dimension

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: trnch, ni, nk

      !@Author B. Bilodeau Feb 2003 - from serdyn5 and phyexe1
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"

      logical, parameter :: OVERWRITE_L = .true.

      character(len=4) :: oname
      integer :: i, k, kam, ivar, nkm1

      real, dimension(ni)     :: work1d
      real, dimension(ni, nk) :: p, work2d1, work2d2

      real, pointer, dimension(:) , contiguous :: zpmoins, busptr2d, ztdew, zflusolaf, &
           zuvsmax, zuvsavg, zhrsmax, zhrsmin, zhusavg, zttmins1, zttmaxs1, zfl
      real, pointer, dimension(:, :), contiguous :: zhuplus, zsigw, ztplus, &
           zuplus, zvplus, zwplus, ztcond, zze, busptr3d, &
           zqcplus, zftot, zfbl, zqtbl, zfdc
      !----------------------------------------------------------------

      if (.not.series_isstep()) return

      call msg_toall(MSG_DEBUG, 'extdiag [BEGIN]')
      if (timings_L) call timing_start_omp(495, 'extdiag', 46)

      nkm1 = nk-1

      MKPTR2D(zfbl, fbl, pvars)
      MKPTR2D(zfdc, fdc, pvars)
      MKPTR1D(zfl, fl, pvars)
      MKPTR1D(zflusolaf, flusolaf, pvars)
      MKPTR2D(zftot, ftot, pvars)
      MKPTR1D(zhrsmax, hrsmax, pvars)
      MKPTR1D(zhrsmin, hrsmin, pvars)
      MKPTR2D(zhuplus, huplus, pvars)
      MKPTR1D(zhusavg, husavg, pvars)
      MKPTR1D(zpmoins, pmoins, pvars)
      MKPTR2D(zqcplus, qcplus, pvars)
      MKPTR2D(zqtbl, qtbl, pvars)
      MKPTR2D(zsigw, sigw, pvars)
      MKPTR2D(ztcond, tcond, pvars)
      MKPTR1D(ztdew, tdew, pvars)
      MKPTR2D(ztplus, tplus, pvars)
      MKPTR1DK(zttmaxs1, ttmax, nk, pvars)
      MKPTR1DK(zttmins1, ttmin, nk, pvars)
      MKPTR2D(zuplus, uplus, pvars)
      MKPTR1D(zuvsavg, uvsavg, pvars)
      MKPTR1D(zuvsmax, uvsmax, pvars)
      MKPTR2D(zvplus, vplus, pvars)
      MKPTR2D(zwplus, wplus, pvars)
      MKPTR2D(zze, ze, pvars)

      if (associated(zfl)) call series_xst(zfl, 'fl', trnch)

      ! Extract time series and zonal diagnostics on nk levels

      if (associated(ztdew)) call series_xst(ztdew,     'TW', trnch) !TDK
      if (associated(zflusolaf)) call series_xst(zflusolaf, 'AF', trnch) !N4

      if (moyhr > 0) then
         if (associated(zuvsmax)) call series_xst(zuvsmax, 'UX', trnch) !UVMX
         if (associated(zuvsavg)) call series_xst(zuvsavg, 'WA', trnch) !UVAV
         if (associated(zhrsmax)) call series_xst(zhrsmax, 'HX', trnch) !HRMX
         if (associated(zhrsmin)) call series_xst(zhrsmin, 'HN', trnch) !HRMN
         if (associated(zhusavg)) call series_xst(zhusavg, 'QV', trnch) !QSAV
         if (associated(zttmins1)) call series_xst(zttmins1, 'T5', trnch) !T5.. +S1!
         if (associated(zttmaxs1)) call series_xst(zttmaxs1, 'T9', trnch) !T9.. +S1!
      endif

      ! Clouds
      if (associated(zftot)) call series_xst(zftot, 'NU', trnch) !FN
      if (associated(zfbl)) call series_xst(zfbl,  'NJ', trnch) !NC
      if (associated(zqtbl)) call series_xst(zqtbl, 'QB', trnch) !QTBL

      if (.not.any(convec == (/'NIL   ', 'SEC   '/))) then
         ! Nuages convectifs
         if (associated(zfdc)) call series_xst(zfdc, 'NC', trnch) !CK
      endif

      ! Effect of precipitation evaporation
      if (series_isvar('EP')) then
         do k=1, nk-1
            do i=1, ni
               work2d1(i, k)=min(0., ztcond(i, k))
            enddo
         enddo
         call series_xst(work2d1, 'EP', trnch)
      endif

      if (associated(zuplus)) call series_xst(zuplus, 'UUWE', trnch) !UP
      if (associated(zvplus)) call series_xst(zvplus, 'VVSN', trnch) !VP
      if (associated(zwplus)) call series_xst(zwplus, 'WW', trnch)   !WP
      if (associated(ztplus)) call series_xst(ztplus, 'TT', trnch) !2T

      ! Potential temperature
      if (series_isvar('TH')) then
         do k=1, nk
!vdir nodep
            do i=1, ni
               p(i, k) = zsigw(i, k)*zpmoins(i)
               work2d1(i, k) = (1.e-5*p(i, k))**(-CAPPA)
               work2d1(i, k) = ztplus(i, k)*work2d1(i, k)
            end do
         end do
         call series_xst(work2d1, 'TH', trnch)
         call series_xst(P, 'PX', trnch)
      endif

      !# Dew point temperature
      if (series_isvar('TD')) then
         ! Find dew point depression first (es). saturation calculated with
         ! respect to water only (since td may be compared to observed tephigram).
         call mhuaes3(work2d1, zhuplus, ztplus, p, .false., ni, nk, ni)
         do k=1, nk
!vdir nodep
            do i=1, ni
               work2d1(i, k) = min( &
                    ztplus(i, k), &
                    ztplus(i, k) - work2d1(i, k) &
                    )
            end do
         end do
         call series_xst(work2d1, 'TD', trnch)
      endif

      !# Moisture: Eliminate negative values
      if (series_isvar('HR')) then
!vdir nodep
         do k=1, nk
            do i=1, ni
               work2d1(i, k) = max(0.0, zhuplus(i, k))
            end do
         end do
         call series_xst(work2d1, 'HU', trnch)
      endif

      !# Relative humidity
      if (series_isvar('HR')) then
         do k=1, nk
            do i=1, ni
               work2d2(i, k) = zsigw(i, k)*zpmoins(i)
            end do
         end do
         call mfohr4 (work2d1, zhuplus, ztplus, &
              work2d2, ni, nk, ni, satuco)
         call series_xst(work2d1, 'HR', trnch)
      endif

      if (stcond /= 'NIL' ) then
         !       Total cloud water
         if (associated(zqcplus)) call series_xst(zqcplus, 'QC', trnch) !CQ
      endif

      if (associated(zpmoins)) call series_xst(zpmoins, 'P0', trnch) !P8

      ! Modulus of the surface wind
      if (series_isvar('VE')) then
         do i = 1, ni
            work1d(i) = sqrt(zuplus(i, nk)*zuplus(i, nk) + &
                 zvplus(i, nk)*zvplus(i, nk))
         end do
         call series_xst(work1d, 'VE', trnch)
      endif

      ! Geopotential height (in dam)
      if (series_isvar('GZ')) then
         do k=1, nk
            do i=1, ni
               work2d1(i, k) = 0.1 * zze(i, k)
            end do
         enddo
         call series_xst(work2d1, 'GZ', trnch)
      endif

      ! Prepare only lowest-level component visibilities for series
      if (stcond(1:6) == 'MP_MY2') then
         do i=1,ni
            work1d(i) = pvars(vis)%data(i+(nkm1-1)*ni)
         enddo
         call series_xst(work1d,    'VIS' , trnch)
         do i=1,ni
            work1d(i) = pvars(vis1)%data(i+(nkm1-1)*ni)
         enddo
         call series_xst(work1d,    'VIS1', trnch)
         do i=1,ni
            work1d(i) = pvars(vis2)%data(i+(nkm1-1)*ni)
         enddo
         call series_xst(work1d,    'VIS2', trnch)
         do i=1,ni
            work1d(i) = pvars(vis3)%data(i+(nkm1-1)*ni)
         enddo
         call series_xst(work1d,    'VIS3', trnch)
      endif

      !# Extract Series for all bus var
      do ivar=1,nphyvars
         oname = pvars(ivar)%meta%oname
         if (series_isvar(oname)) then
            kam    = pvars(ivar)%meta%nk
            if (.not.associated(pvars(ivar)%data)) cycle
            if (kam == 1) then
               busptr2d(1:ni) => pvars(ivar)%data(:)
               call series_xst(busptr2d, oname, trnch, &
                    F_overwrite_L=.not.OVERWRITE_L)
            else
               busptr3d(1:ni,1:kam) => pvars(ivar)%data(:)
               call series_xst(busptr3d, oname, trnch, &
                    F_overwrite_L=.not.OVERWRITE_L)
            endif
        endif
      enddo
      
      if (timings_L) call timing_stop_omp(495)
      call msg_toall(MSG_DEBUG, 'extdiag [END]')
      !----------------------------------------------------------------
      return
   end subroutine extdiag3

end module extdiag
