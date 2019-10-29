!-------------------------------------- LICENCE BEGIN ---------------------------
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
!-------------------------------------- LICENCE END -----------------------------

module precipitation
   implicit none
   private
   public :: precipitation4

contains

   !/@*
   subroutine precipitation4(tplus0, huplus0, dbus, dsiz, fbus, fsiz, vbus, vsiz,&
        dt, ni, nk, kount, trnch)
      use debug_mod, only: init2nan
      use tdpack_const, only: RAUW
      use cnv_main, only: cnv_main3
      use condensation, only: condensation3
      use phy_status, only: phy_error_L
      use phy_options
      use phybus
      use tendency, only: apply_tendencies
      implicit none
!!!#include <arch_specific.hf>
      !@Object Interface to convection/condensation
      !@Arguments
      !
      !          - Input -
      ! dsiz     dimension of dbus
      ! fsiz     dimension of fbus
      ! vsiz     dimension of vbus
      ! dt       timestep (sec.)
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! kount    timestep number
      ! trnch    slice number
      !
      !          - Input/Output -
      ! dbus     dynamics input field
      ! fbus     historic variables for the physics
      ! vbus     physics tendencies and other output fields from the physics

      integer, intent(in) :: fsiz,vsiz,dsiz,ni,nk,kount,trnch
      real,    intent(in) :: dt
      real, dimension(ni,nk), intent(inout) :: tplus0,huplus0
      real, target, intent(inout) :: dbus(dsiz), fbus(fsiz), vbus(vsiz)

      !@Author L.Spacek, November 2011
      !*@/
#include <msg.h>
#include "phymkptr.hf"
      include "nocld.cdk"

      ! Local variable declarations
      integer :: nkm1
      integer, dimension(ni,nk-1) :: ilab
      real :: irhow
      real, dimension(ni) :: beta
      real, dimension(ni,nk) :: press
      real, dimension(ni,nk-1) :: t0,q0,qc0

      real, dimension(:), pointer :: zrckfc,ztlc,ztlcs,ztls,ztsc,ztscs,ztss,zpmoins,ztlcm
      real, dimension(:,:), pointer :: ztcond,zhucond,ztshal,zhushal,zqcphytd,zqrphytd, &
           zcte,zcqe,zste,zsqe,zcqce,zsqce,zprcten,zsqre,zhumoins,zsigt,zmte,zmqe, &
           zmqce,zprctnm
      logical :: lkfbe
      !----------------------------------------------------------------

      ! Early return moist processes are inactive
      if (convec == 'NIL' .and. stcond == 'NIL' .and. conv_shal == 'NIL') return
      call msg_toall(MSG_DEBUG, 'precipitation [BEGIN]')

      ! Assign local bus pointers
      MKPTR2D(zcqce, cqce, vbus)
      MKPTR2D(zcqe, cqe, vbus)
      MKPTR2D(zcte, cte, vbus)
      MKPTR2D(zhucond, hucond, vbus)
      MKPTR2D(zmqe, mqe, vbus)
      MKPTR2D(zhumoins, humoins, dbus)
      MKPTR2D(zhushal, hushal, vbus)
      MKPTR1D(zpmoins, pmoins, fbus)
      MKPTR2D(zprcten, prcten, fbus)
      MKPTR2D(zprctnm, prctnm, vbus)
      MKPTR2D(zmqce, mqce, vbus)
      MKPTR2D(zqcphytd, qcphytd, vbus)
      MKPTR2D(zqrphytd, qrphytd, vbus)
      MKPTR1D(zrckfc, rckfc, fbus)
      MKPTR2D(zsigt, sigt, dbus)
      MKPTR2D(zsqce, sqce, vbus)
      MKPTR2D(zsqe, sqe, vbus)
      MKPTR2D(zsqre, sqre, vbus)
      MKPTR2D(zste, ste, vbus)
      MKPTR2D(ztcond, tcond, vbus)
      MKPTR1D(ztlc, tlc, fbus)
      MKPTR1D(ztlcm, tlcm, vbus)
      MKPTR1D(ztlcs, tlcs, vbus)
      MKPTR1D(ztls, tls, fbus)
      MKPTR2D(zmte, mte, vbus)
      MKPTR1D(ztsc, tsc, fbus)
      MKPTR1D(ztscs, tscs, vbus)
      MKPTR2D(ztshal, tshal, vbus)
      MKPTR1D(ztss, tss, fbus)

      call init2nan(beta)
      call init2nan(t0, q0, qc0, press)

      ! Local initialization
      nkm1 = nk-1
      lkfbe  = any(convec == (/ &
           'BECHTOLD', &
           'KFC     ', &
           'KFC2    ', &
           'KFC3    ' &
           /))

      ! Run deep and shallow convective schemes
      if (timings_L) call timing_start_omp(435, 'cnv_main', 46)
      call cnv_main3(dbus, dsiz, fbus, fsiz, vbus, vsiz,   &
           t0, q0, qc0, ilab, beta, dt, ni, nk, &
           kount)
      if (timings_L) call timing_stop_omp(435)
      if (phy_error_L) return

      ! Run the gridscale condensation / microphysics scheme
      if (timings_L) call timing_start_omp(440, 'condensation', 46)
      call condensation3(dbus, dsiz, fbus, fsiz, vbus, vsiz, &
           tplus0, t0, huplus0, q0, qc0, ilab, beta, &
           dt, ni, nk, kount, trnch)
      if (timings_L) call timing_stop_omp(440)
      if (phy_error_L) return

      ! Eliminate cloud tendencies at upper levels and in dry regions on user request
      if (climat .or. stratos) then
         press(:,:) = zsigt(:,:) * spread(zpmoins(:), dim=2, ncopies=nk)
         where(press < TOPC .or. zhumoins <= MINQ)
            zcte = 0.
            zste = 0.
            zmte = 0.
            zcqe = 0.
            zsqe = 0.
            zmqe = 0.
            zcqce = 0.
            zsqce = 0.
            zmqce = 0.
            zsqre = 0.
         endwhere
         if (associated(zprcten)) then
            where(press < TOPC .or. zhumoins <= MINQ)
               zprcten = 0.
            endwhere
         endif
      endif

      ! Total tendencies are the sum of convective and stratiform
      ztcond(:,1:nkm1) = zcte(:,1:nkm1) + zste(:,1:nkm1) + zmte(:,1:nkm1)
      zhucond(:,1:nkm1) = zcqe(:,1:nkm1) + zsqe(:,1:nkm1) + zmqe(:,1:nkm1)
      if (lkfbe.and.stcond(1:3)=='MP_') then !Use only convective liquid fraction for MP schemes
         zqcphytd(:,1:nkm1) = zprcten(:,1:nkm1) + zsqce(:,1:nkm1) + zprctnm(:,1:nkm1)
      else
         zqcphytd(:,1:nkm1) = zcqce(:,1:nkm1) + zsqce(:,1:nkm1) + zmqce(:,1:nkm1)
      endif
      zqrphytd(:,1:nkm1) = zsqre(:,1:nkm1)

      ! Apply shallow convective tendencies when computed before deep
      if (conv_shal == 'BECHTOLD') then
         call apply_tendencies(dbus,vbus,fbus,tplus,tshal,ni,nk,nkm1)
         call apply_tendencies(dbus,vbus,fbus,huplus,hushal,ni,nk,nkm1)
         call apply_tendencies(dbus,vbus,fbus,qcplus,prctns,ni,nk,nkm1)
         call apply_tendencies(dbus,vbus,fbus,qcplus,pritns,ni,nk,nkm1)
      endif

      ! Apply deep CPS and condensation tendencies for non-microphysics schemes
      if (.not.any(stcond == (/'MP_MY2','MP_P3 '/))) then
         call apply_tendencies(dbus,vbus,fbus,tplus,tcond,ni,nk,nkm1)
         call apply_tendencies(dbus,vbus,fbus,huplus,hucond,ni,nk,nkm1)
         call apply_tendencies(dbus,vbus,fbus,qcplus,qcphytd,ni,nk,nkm1)
      endif

      ! Add shallow convection tendencies to convection/condensation tendencies
      ztcond(:,1:nkm1) = ztcond(:,1:nkm1) + ztshal(:,1:nkm1)
      zhucond(:,1:nkm1) = zhucond(:,1:nkm1) + zhushal(:,1:nkm1)

      ! Convert from flux to liquid-equivalent precipitation rates by dividing by
      ! the density of water
      irhow = 1./RAUW
      zrckfc = zrckfc * irhow
      ztlc   = ztlc   * irhow
      ztlcm  = ztlcm  * irhow
      ztls   = ztls   * irhow
      ztsc   = ztsc   * irhow
      ztss   = ztss   * irhow

      ! Reset precipitation rates to 0 for kount==0 physics initialization, unless
      ! using a microphysics scheme
      if (kount == 0 .and. stcond(1:2) /= 'MP') then
         zrckfc = 0.
         ztlc   = 0.
         ztlcm  = 0.
         ztlcs  = 0.
         ztls   = 0.
         ztsc   = 0.
         ztscs  = 0.
         ztss   = 0.
      endif

      call msg_toall(MSG_DEBUG, 'precipitation [END]')
      !----------------------------------------------------------------
      return
   end subroutine precipitation4

end module precipitation
