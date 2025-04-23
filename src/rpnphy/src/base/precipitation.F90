
module precipitation
   implicit none
   private
   public :: precipitation4

contains

   !/@*
   subroutine precipitation4(pvars, dt, kount, ni, nk)
      use debug_mod, only: init2nan
      use tdpack_const, only: RAUW
      use shal_ktrsnt, only: sc_ktrsnt
      use cnv_main, only: cnv_main4
      use condensation, only: condensation4
      use microphy_statcond, only: sc_adjust
      use phy_status, only: phy_error_L
      use phy_options
      use phybusidx
      use phymem, only: phyvar
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
      !
      !          - Input/Output -
      ! pvars    list of all phy vars (meta + slab data)

      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni,nk,kount
      real,    intent(in) :: dt

      !@Author L.Spacek, November 2011
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"
      include "nocld.cdk"

      ! Local variable declarations
      integer :: nkm1
      real :: irhow
      real, dimension(ni,nk) :: press, t0, q0, qc0

      real, dimension(:), pointer, contiguous :: zrckfc,ztlc,ztlcs,ztls,ztsc,ztscs,ztss,zpmoins,ztlcm, ztdmaskxdt
      real, dimension(:,:), pointer, contiguous :: ztcond,zhucond,ztshal,zhushal,zqcphytd,zqrphytd, &
           zcte,zcqe,zste,zsqe,zcqce,zsqce,zprcten,zsqre,zhumoins,zsigt,zmte,zmqe, &
           zmqce,zprctnm,ztplus,zhuplus,zqcplus,zqcmoins,zprctns,zpritns,zqccondc2, &
           zqcondc2,ztcondc2
      !----------------------------------------------------------------

      ! Early return moist processes are inactive
      if (convec == 'NIL' .and. stcond == 'NIL' .and. conv_shal == 'NIL') return
      call msg_toall(MSG_DEBUG, 'precipitation [BEGIN]')

      MKPTR2D(zcqce, cqce, pvars)
      MKPTR2D(zcqe, cqe, pvars)
      MKPTR2D(zcte, cte, pvars)
      MKPTR2D(zhucond, hucond, pvars)
      MKPTR2D(zmqe, mqe, pvars)
      MKPTR2D(zhumoins, humoins, pvars)
      MKPTR2D(zhuplus, huplus, pvars)
      MKPTR2D(zhushal, hushal, pvars)
      MKPTR1D(zpmoins, pmoins, pvars)
      MKPTR2D(zprcten, prcten, pvars)
      MKPTR2D(zprctnm, prctnm, pvars)
      MKPTR2D(zmqce, mqce, pvars)
      MKPTR2D(zqccondc2, qccondc2, pvars)
      MKPTR2D(zqcmoins, qcmoins, pvars)
      MKPTR2D(zqcondc2, qcondc2, pvars)
      MKPTR2D(zqcphytd, qcphytd, pvars)
      MKPTR2D(zqcplus, qcplus, pvars)
      MKPTR2D(zqrphytd, qrphytd, pvars)
      MKPTR1D(zrckfc, rckfc, pvars)
      MKPTR2D(zsigt, sigt, pvars)
      MKPTR2D(zsqce, sqce, pvars)
      MKPTR2D(zsqe, sqe, pvars)
      MKPTR2D(zsqre, sqre, pvars)
      MKPTR2D(zste, ste, pvars)
      MKPTR2D(ztcond, tcond, pvars)
      MKPTR2D(ztcondc2, tcondc2, pvars)
      MKPTR1D(ztlc, tlc, pvars)
      MKPTR1D(ztlcm, tlcm, pvars)
      MKPTR1D(ztlcs, tlcs, pvars)
      MKPTR1D(ztls, tls, pvars)
      MKPTR2D(zmte, mte, pvars)
      MKPTR2D(ztplus, tplus, pvars)
      MKPTR1D(ztsc, tsc, pvars)
      MKPTR1D(ztscs, tscs, pvars)
      MKPTR2D(ztshal, tshal, pvars)
      MKPTR1D(ztss, tss, pvars)
      
      MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)
      MKPTR2D(zprctns, prctns, pvars)
      MKPTR2D(zpritns, pritns, pvars)

      ! Local initialization
      nkm1 = nk-1
      call init2nan(t0, q0, qc0, press)
      
      ! Run "special" pre-convection shallow scheme
      call sc_ktrsnt(pvars, dt, ni, nkm1)
      if (phy_error_L) return

      ! Store current state
      if (.not.any(stcond == (/&
           'MP_MY2 ', &
           'MP_P3  ', &
           'MP_P3V3'/))) then
         t0(:,:) = ztplus(:,:)
         q0(:,:) = zhuplus(:,:)
         if (stcond /= 'NIL' .and. associated(zqcplus)) then
            qc0(:,:) = zqcplus(:,:)
         else
            qc0(:,:) = 0.
         endif
      endif

      ! Clipping
      where (zhuplus(:,:) < 0.) zhuplus(:,:) = 0.
      if (stcond /= 'NIL') then
         where(zqcplus(:,:) < 0.) zqcplus(:,:) = 0.
         where(zqcmoins(:,:) < 0.) zqcmoins(:,:) = 0.
      endif
      
      ! Run deep and shallow convective schemes
      if (timings_L) call timing_start_omp(435, 'cnv_main', 46)
      call cnv_main4(pvars, dt, kount, ni, nk)
      if (timings_L) call timing_stop_omp(435)
      if (phy_error_L) return

      ! Run the gridscale condensation / microphysics scheme
      if (timings_L) call timing_start_omp(440, 'condensation', 46)
      call condensation4(pvars, dt, kount, ni, nk)
      if (timings_L) call timing_stop_omp(440)
      if (phy_error_L) return

      ! Eliminate cloud,deep and mid  tendencies at upper levels and in dry regions on user request
      if (climat .or. stratos) then
         press(:,:) = zsigt(:,:) * spread(zpmoins(:), dim=2, ncopies=nk)
         where(press < TOPC)
            zcte = 0.
            zste = 0.
            zcqe = 0.
            zsqe = 0.
            zcqce = 0.
            zsqce = 0.
            zsqre = 0.
         endwhere
         if (conv_mid /= 'NIL') then
            where(press < TOPC)
               zmte = 0.
               zmqe = 0.
               zmqce = 0.
            endwhere
         endif
         if (associated(zprcten)) then
            where(press < TOPC)
               zprcten = 0.
            endwhere
         endif
      endif

      ! Sum of convective and stratiform 
      ztcond(:,1:nkm1) = zcte(:,1:nkm1) + zste(:,1:nkm1)
      if (associated(zmte)) &
         ztcond(:,1:nkm1) = ztcond(:,1:nkm1) + zmte(:,1:nkm1)
      zhucond(:,1:nkm1) = zcqe(:,1:nkm1) + zsqe(:,1:nkm1)
      if (associated(zmqe)) &
           zhucond(:,1:nkm1) = zhucond(:,1:nkm1) + zmqe(:,1:nkm1)
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
      if (.not.any(stcond == (/&
           'MP_MY2 ', &
           'MP_P3  ', &
           'MP_P3V3'/))) then
         ztplus(:,1:nkm1) = t0(:,1:nkm1)
         zhuplus(:,1:nkm1) = q0(:,1:nkm1)
         call apply_tendencies(ztplus, zhuplus, &
              &                ztcond, zhucond, ztdmaskxdt,ni,nk,nkm1)
         if (associated(zqcplus)) then
            zqcplus(:,1:nkm1) = qc0(:,1:nkm1)
            call apply_tendencies(zqcplus,zqcphytd,ztdmaskxdt,ni,nk,nkm1)
         endif
        ! Apply BECHTOLD shallow convective tendencies 
        if (conv_shal == 'BECHTOLD') then
           call apply_tendencies(ztplus, zhuplus,&
                &                ztshal, zhushal, ztdmaskxdt,ni,nk,nkm1)
           if (associated(zqcplus)) then
              call apply_tendencies(zqcplus, zprctns, ztdmaskxdt,ni,nk,nkm1)
              call apply_tendencies(zqcplus, zpritns, ztdmaskxdt,ni,nk,nkm1)
           endif
        endif
      endif

      ! Post-scheme condensation adjustment
      if (stcond == 'S2') &
           call sc_adjust(ztcondc2, zqcondc2, zqccondc2, pvars, delt, ni, nkm1)
      
      ! Add shallow convection tendencies to convection/condensation tendencies 
      if (associated(ztshal))  ztcond(:,1:nkm1)  = ztcond(:,1:nkm1)  + ztshal(:,1:nkm1)
      if (associated(zhushal)) zhucond(:,1:nkm1) = zhucond(:,1:nkm1) + zhushal(:,1:nkm1)
      
           
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
