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
   subroutine precipitation4(dbus, fbus, vbus, dt, ni, nk, kount, trnch)
      use debug_mod, only: init2nan
      use tdpack_const, only: RAUW
      use shal_ktrsnt, only: sc_ktrsnt
      use cnv_main, only: cnv_main4
      use condensation, only: condensation4
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

      integer, intent(in) :: ni,nk,kount,trnch
      real,    intent(in) :: dt
      real, dimension(:), pointer, contiguous :: dbus, fbus, vbus

      !@Author L.Spacek, November 2011
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"
      include "nocld.cdk"

      ! Local variable declarations
      integer :: nkm1
      real :: irhow
      real, dimension(ni,nk) :: press, t0, q0, qc0

      real, dimension(:), pointer, contiguous :: zrckfc,ztlc,ztlcs,ztls,ztsc,ztscs,ztss,zpmoins,ztlcm
      real, dimension(:,:), pointer, contiguous :: ztcond,zhucond,ztshal,zhushal,zqcphytd,zqrphytd, &
           zcte,zcqe,zste,zsqe,zcqce,zsqce,zprcten,zsqre,zhumoins,zsigt,zmte,zmqe, &
           zmqce,zprctnm,ztplus,zhuplus,zqcplus,zqcmoins
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
      MKPTR2D(zhuplus, huplus, dbus)
      MKPTR2D(zhushal, hushal, vbus)
      MKPTR1D(zpmoins, pmoins, fbus)
      MKPTR2D(zprcten, prcten, fbus)
      MKPTR2D(zprctnm, prctnm, vbus)
      MKPTR2D(zmqce, mqce, vbus)
      MKPTR2D(zqcmoins, qcmoins, dbus)
      MKPTR2D(zqcphytd, qcphytd, vbus)
      MKPTR2D(zqcplus, qcplus, dbus)
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
      MKPTR2D(ztplus, tplus, dbus)
      MKPTR1D(ztsc, tsc, fbus)
      MKPTR1D(ztscs, tscs, vbus)
      MKPTR2D(ztshal, tshal, vbus)
      MKPTR1D(ztss, tss, fbus)

      ! Local initialization
      nkm1 = nk-1
      call init2nan(t0, q0, qc0, press)
      
      ! Run "special" pre-convection shallow scheme
      call sc_ktrsnt(dbus, fbus, vbus, dt, ni, nkm1)
      if (phy_error_L) return

      ! Store current state
      t0(:,:) = ztplus(:,:)
      q0(:,:) = zhuplus(:,:)
      where (zhuplus(:,:) < 0.) zhuplus(:,:) = 0.
      if (stcond /= 'NIL' ) then
         qc0(:,:) = zqcplus(:,:)
         where(zqcplus(:,:) < 0.) zqcplus(:,:) = 0.
         where(zqcmoins(:,:) < 0.) zqcmoins(:,:) = 0.
      else
         qc0(:,:) = 0.
      endif
      
      ! Run deep and shallow convective schemes
      if (timings_L) call timing_start_omp(435, 'cnv_main', 46)
      call cnv_main4(dbus, fbus, vbus, dt, ni, nk, kount)
      if (timings_L) call timing_stop_omp(435)
      if (phy_error_L) return

      ! Run the gridscale condensation / microphysics scheme
      if (timings_L) call timing_start_omp(440, 'condensation', 46)
      call condensation4(dbus, fbus, vbus, dt, ni, nk, kount, trnch)
      if (timings_L) call timing_stop_omp(440)
      if (phy_error_L) return

      ! Eliminate cloud,deep and mid  tendencies at upper levels and in dry regions on user request
      if (climat .or. stratos) then
         press(:,:) = zsigt(:,:) * spread(zpmoins(:), dim=2, ncopies=nk)
         where(press < TOPC)
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
            where(press < TOPC)
               zprcten = 0.
            endwhere
         endif
      endif

      ! Sum of convective and stratiform 
      ztcond(:,1:nkm1) = zcte(:,1:nkm1) + zste(:,1:nkm1) + zmte(:,1:nkm1)
      zhucond(:,1:nkm1) = zcqe(:,1:nkm1) + zsqe(:,1:nkm1) + zmqe(:,1:nkm1)
      if (associated(zqcphytd)) zqcphytd(:,1:nkm1) = 0.
      if (stcond(1:3)=='MP_') then
         ! Use only convective liquid fraction for MP schemes
         if (associated(zprcten)) zqcphytd(:,1:nkm1) = zqcphytd(:,1:nkm1) + zprcten(:,1:nkm1)
         if (associated(zsqce)) zqcphytd(:,1:nkm1) = zqcphytd(:,1:nkm1) + zsqce(:,1:nkm1)
         if (associated(zprctnm)) zqcphytd(:,1:nkm1) = zqcphytd(:,1:nkm1) + zprctnm(:,1:nkm1)
      elseif (associated(zqcphytd)) then
         ! Use full convective condensate tendency for condensation schemes
         if (associated(zcqce)) zqcphytd(:,1:nkm1) = zqcphytd(:,1:nkm1) + zcqce(:,1:nkm1)
         if (associated(zsqce)) zqcphytd(:,1:nkm1) = zqcphytd(:,1:nkm1) + zsqce(:,1:nkm1)
         if (associated(zmqce)) zqcphytd(:,1:nkm1) = zqcphytd(:,1:nkm1) + zmqce(:,1:nkm1)
      endif
      if (associated(zqrphytd)) zqrphytd(:,1:nkm1) = zsqre(:,1:nkm1)
      ! note: zqcphytd does not contain shal contributions; ok because will be re-calculated in tendency

      ! Note: at end of condensation we do: if (.not.any(stcond == (/'MP_MY2','MP_P3 '/))) then 
      ! tplus=t0,huplus=hu0,qcp=q0 i.e. go back to thermo state AFTER KTRSNT
      ! must now re-apply masked(if climat .or.stratos) conv and cond tendencies
      ! must NOT re-apply shallow conv tendencies if shal == ktrsnt
      if (.not.any(stcond == (/'MP_MY2', 'MP_P3 '/))) then
         ztplus(:,1:nkm1) = t0(:,1:nkm1)
         zhuplus(:,1:nkm1) = q0(:,1:nkm1)
         call apply_tendencies(dbus,vbus,fbus,tplus,tcond,ni,nk,nkm1)
         call apply_tendencies(dbus,vbus,fbus,huplus,hucond,ni,nk,nkm1)
         if (associated(zqcplus)) then
            zqcplus(:,1:nkm1) = qc0(:,1:nkm1)
            call apply_tendencies(dbus,vbus,fbus,qcplus,qcphytd,ni,nk,nkm1)
         endif
        ! Apply BECHTOLD shallow convective tendencies 
        if (conv_shal == 'BECHTOLD') then
           call apply_tendencies(dbus,vbus,fbus,tplus,tshal,ni,nk,nkm1)
           call apply_tendencies(dbus,vbus,fbus,huplus,hushal,ni,nk,nkm1)
           if (associated(zqcplus)) then
              call apply_tendencies(dbus,vbus,fbus,qcplus,prctns,ni,nk,nkm1)
              call apply_tendencies(dbus,vbus,fbus,qcplus,pritns,ni,nk,nkm1)
           endif
        endif
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
