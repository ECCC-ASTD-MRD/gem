
module cnv_main
   implicit none
   private
   public :: cnv_main4

contains

   !/@*
   subroutine cnv_main4(pvars, dt, kount, ni, nk)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: CHLC, CPD, GRAV, RAUW
      use phy_status, only: phy_error_L, PHY_OK
      use phybusidx
      use phymem, only: phyvar
      use phy_options
      use cnv_options
      use kfrpn, only: kfrpn4
      use phybudget, only: pb_compute, pb_conserve, pb_residual
      use shallconv, only: shallconv5
      use kfmid, only: kfmid1
      use tendency, only: apply_tendencies
      use conv_mp_tendencies, only: conv_mp_tendencies1
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@objective Interface to convection/condensation
      !@Arguments
      !          - Input -
      ! dt       timestep (parameter) (sec.)
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! kount    timestep number
      !
      !          - Input/Output -
      ! pvars    list of all phy vars (meta + slab data)
 
      type(phyvar), pointer, contiguous :: pvars(:)
      integer, intent(in) :: ni, nk, kount
      real, intent(in) :: dt
      !@Author L.Spacek, November 2011
      !@Revisions
      ! 001   -PV/JM-nov2014- fix communication between deep convection and MY_DM
      ! 002   -PV-dec2014 -
      ! 003   -JY-PV-Jul2015- Separate deep and shallow convection scheme in unified bechtold scheme
      ! 004   -JY-PV-Nov2016- Bkf-shallow only checked if no deep convection is active
      !*@/
#include <rmn/msg.h>
#include "phymkptr.hf"

      include "surface.cdk"

      real    :: cdt1, rcdt1
      integer :: i,k,niter, ier, nkm1

      real(REAL64), dimension(ni) :: l_en0, l_pw0
      !*Set Dimensions for convection call
      !
      !INTEGER KENS = 0  ! number of additional ensemble members for convection
      ! (apart from base scheme) [ 1 - 3 ]
      ! NOTA: for base scheme (1 deep + 1 shallow draft)
      !       just put KENS = 0. This will be more
      !       efficient numerically but results will be
      !       less smooth
      !       RECOMMENDATION: Use KENS = 3 for smooth
      !           climate and NWP prediction runs
      integer, parameter :: kchm = 1 ! maximum number of chemical tracers

      integer:: kidia = 1        ! start index for  horizontal computations !#TODO - replace by "1"
      integer:: kbdia = 1        ! start index for vertical computations    !#TODO - replace by "1"
      integer:: ktdia = 1        ! end index for vertical computations      !#TODO - replace by "1"
      ! defined as klev + 1 - ktdia

      !*Define convective in and output arrays
      !
      integer, dimension(ni)           :: kcount
      real,    dimension(ni)           :: cnv_active, cape
      real,    dimension(ni,nk-1)      :: dummy1, dummy2, geop, cond_prflux, qtl, qts
      real,    dimension(ni,nk-1)      :: ppres,pudr, pddr, phsflx, dmsedt
      real,    dimension(ni,nk-1,kchm) :: pch1, pch1ten


      ! Switches for convection call
      !
      !integer:: kice = 0           ! take ice phase into account or not (0)
      !logical:: ldown    = .true.  ! switch for convective downdrafts
      !logical:: luvconv  = .false. ! account for transport of hor. wind
      ! attention: should be set to false as not yet sufficiently validated
      !logical:: lch1conv = .false. ! account for transport of tracers
      ! in the convection code
      !logical:: lsettadj = .false. ! set adjustment times
      ! (otherwise it is computed automatically)
      !real   :: xtadjd=3600., xtadjs=10800. ! only used if lsettadj = .true.


      ! Pointer to buses

      real, pointer, dimension(:), contiguous :: psm, psp, ztdmaskxdt, zfcpflg, zkkfc, zrckfc, &
           ztlc, ztsc, zkshal, zwstar, zconedc, zconesc, zconqdc, zconqsc, ztlcs, ztscs, &
           ztstar, zcapekfc, ztauckfc, zcinkfc, zmg, zml, zdlat, zdxdy, zkmid, &
           zabekfc, zpeffkfc, zrice_int, zrliq_int, zwumaxkfc, zzbasekfc, zztopkfc, &
           zcoadvu,zcoadvv,zcoage,zcowlcl,zwklcl,zcozlcl,ztlcm, zconemc, zconqmc, &
           zmcd,zmpeff,zmainc
      real, pointer, dimension(:,:), contiguous :: ncp, nip, qcm, qcp, qip, qqm, qqp, qrp, &
           sigma, ttm, ttp, uu, vv, wz, zfdc, zgztherm, zhufcp, zhushal, &
           zprcten, zpriten, zqckfc, ztfcp, ztshal, ztusc, ztvsc,  &
           zufcp, zvfcp, zufcp1, zvfcp1, zufcp2, zvfcp2, zufcp3, zvfcp3,  zsufcp, zsvfcp, &
           zprctns,zpritns,qti1p, nti1p, zti1p, zfsc, zqlsc, zqssc, &
           zumfs, ztpostshal, zhupostshal, zcqce, zcqe, zcte, zen, zkt, &
           zareaup, zdmfkfc, zkfcrf, zkfcsf, zqldi, zqrkfc, zqsdi, zumfkfc, &
           zqdifv, zwklclplus, zkfmrf, zkfmsf, zqlmi, zqsmi, zfmc, zprctnm, zpritnm, &
           zmqce, zmte, zmqe, zumid, zvmid, zrnflx, zsnoflx, zmrk2

      logical :: kfcmom
      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'cnv_main [BEGIN]')

      nkm1 = nk - 1
      kfcmom = (cmt_type_i /= CMT_NONE)

      MKPTR1D(psm, pmoins, pvars)
      MKPTR1D(psp, pmoins, pvars)
      MKPTR1D(zabekfc, abekfc, pvars)
      MKPTR1D(zcapekfc, capekfc, pvars)
      MKPTR1D(ztauckfc, tauckfc, pvars)
      MKPTR1D(zcinkfc, cinkfc, pvars)
      MKPTR1D(zcoadvu, coadvu, pvars)
      MKPTR1D(zcoadvv, coadvv, pvars)
      MKPTR1D(zcoage, coage, pvars)
      MKPTR1D(zconedc, conedc, pvars)
      MKPTR1D(zconemc, conemc, pvars)
      MKPTR1D(zconesc, conesc, pvars)
      MKPTR1D(zconqdc, conqdc, pvars)
      MKPTR1D(zconqmc, conqmc, pvars)
      MKPTR1D(zconqsc, conqsc, pvars)
      MKPTR1D(zcowlcl, cowlcl, pvars)
      MKPTR1D(zcozlcl, cozlcl, pvars)
      MKPTR1D(zdlat, dlat, pvars)
      MKPTR1D(zdxdy, dxdy, pvars)
      MKPTR1D(zfcpflg, fcpflg, pvars)
      MKPTR1D(zkkfc, kkfc, pvars)
      MKPTR1D(zkmid, kmid, pvars)
      MKPTR1D(zkshal, kshal, pvars)
      MKPTR1D(zmainc, mainc, pvars)
      MKPTR1D(zmcd, mcd, pvars)
      MKPTR1D(zmg, mg, pvars)
      MKPTR1D(zml, ml, pvars)
      MKPTR1D(zmpeff, mpeff, pvars)
      MKPTR1D(zpeffkfc, peffkfc, pvars)
      MKPTR1D(zrckfc, rckfc, pvars)
      MKPTR1D(ztlcm, tlcm, pvars)
      MKPTR1D(zrice_int, rice_int, pvars)
      MKPTR1D(zrliq_int, rliq_int, pvars)
      MKPTR1D(ztdmaskxdt, tdmaskxdt, pvars)
      MKPTR1D(ztlc, tlc, pvars)
      MKPTR1D(ztlcs, tlcs, pvars)
      MKPTR1D(ztsc, tsc, pvars)
      MKPTR1D(ztscs, tscs, pvars)
      MKPTR1D(ztstar, tstar, pvars)
      MKPTR1D(zwklcl, wklcl, pvars)
      MKPTR1D(zwstar, wstar, pvars)
      MKPTR1D(zwumaxkfc, wumaxkfc, pvars)
      MKPTR1D(zzbasekfc, zbasekfc, pvars)
      MKPTR1D(zztopkfc, ztopkfc, pvars)

      MKPTR2D(sigma, sigw, pvars)

      MKPTR2Dm1(zti1p, zti1plus, pvars)
      MKPTR2Dm1(nti1p, nti1plus, pvars)
      MKPTR2Dm1(qti1p, qti1plus, pvars)
      MKPTR2Dm1(ncp, ncplus, pvars)
      MKPTR2Dm1(nip, niplus, pvars)
      MKPTR2Dm1(qcm, qcmoins, pvars)
      MKPTR2Dm1(qcp, qcplus, pvars)
      MKPTR2Dm1(qip, qiplus, pvars)
      MKPTR2Dm1(qqm, humoins, pvars)
      MKPTR2Dm1(qqp, huplus, pvars)
      MKPTR2Dm1(ttm, tmoins, pvars)
      MKPTR2Dm1(qrp, qrplus, pvars)
      MKPTR2Dm1(ttp, tplus, pvars)
      MKPTR2Dm1(uu, uplus, pvars)
      MKPTR2Dm1(vv, vplus, pvars)
      MKPTR2Dm1(wz, wplus, pvars)
      MKPTR2Dm1(zcqce, cqce, pvars)
      MKPTR2Dm1(zcqe, cqe, pvars)
      MKPTR2Dm1(zcte, cte, pvars)
      MKPTR2Dm1(zfdc, fdc, pvars)
      MKPTR2Dm1(zfmc, fmc, pvars)
      MKPTR2Dm1(zfsc, fsc, pvars)
      MKPTR2Dm1(zgztherm, gztherm, pvars)
      MKPTR2Dm1(zhufcp, hufcp, pvars)
      MKPTR2Dm1(zmqe, mqe, pvars)
      MKPTR2Dm1(zhupostshal, hupostshal, pvars)
      MKPTR2Dm1(zhushal, hushal, pvars)
      MKPTR2Dm1(zkt, kt, pvars)
      MKPTR2Dm1(zprcten, prcten, pvars)
      MKPTR2Dm1(zprctnm, prctnm, pvars)
      MKPTR2Dm1(zprctns, prctns, pvars)
      MKPTR2Dm1(zpriten, priten, pvars)
      MKPTR2Dm1(zpritnm, pritnm, pvars)
      MKPTR2Dm1(zpritns, pritns, pvars)
      MKPTR2Dm1(zqckfc, qckfc, pvars)
      MKPTR2Dm1(zmqce, mqce, pvars)
      MKPTR2Dm1(zqdifv, qdifv, pvars)
      MKPTR2Dm1(zqlsc, qlsc, pvars)
      MKPTR2Dm1(zqssc, qssc, pvars)
      
      MKPTR2Dm1(zsufcp, sufcp, pvars)
      MKPTR2Dm1(zsvfcp, svfcp, pvars)
      
      MKPTR2Dm1(ztfcp, tfcp, pvars)
      MKPTR2Dm1(zmte, mte, pvars)
      MKPTR2Dm1(zrnflx, rnflx, pvars)
      MKPTR2Dm1(zsnoflx, snoflx, pvars)
      MKPTR2Dm1(ztpostshal, tpostshal, pvars)
      MKPTR2Dm1(ztshal, tshal, pvars)
      MKPTR2Dm1(ztusc, tusc, pvars)
      MKPTR2Dm1(ztvsc, tvsc, pvars)
      
      MKPTR2Dm1(zufcp, ufcp, pvars)
      MKPTR2Dm1(zufcp1, ufcp1, pvars)
      MKPTR2Dm1(zufcp2, ufcp2, pvars)
      MKPTR2Dm1(zufcp3, ufcp3, pvars)
      
      MKPTR2Dm1(zumfs, umfs, pvars)
      MKPTR2Dm1(zumid, umid, pvars)
      
      MKPTR2Dm1(zvfcp, vfcp, pvars)
      MKPTR2Dm1(zvfcp1, vfcp1, pvars)
      MKPTR2Dm1(zvfcp2, vfcp2, pvars)
      MKPTR2Dm1(zvfcp3, vfcp3, pvars)
      
      MKPTR2Dm1(zwklclplus, wklclplus, pvars)
      MKPTR2Dm1(zvmid, vmid, pvars)
      MKPTR2Dm1(zareaup, areaup, pvars)
      MKPTR2Dm1(zdmfkfc, dmfkfc, pvars)
      MKPTR2Dm1(zkfcrf, kfcrf, pvars)
      MKPTR2Dm1(zkfcsf, kfcsf, pvars)
      MKPTR2Dm1(zkfmrf, kfmrf, pvars)
      MKPTR2Dm1(zkfmsf, kfmsf, pvars)
      MKPTR2Dm1(zqldi, qldi, pvars)
      MKPTR2Dm1(zqlmi, qlmi, pvars)
      MKPTR2Dm1(zqrkfc, qrkfc, pvars)
      MKPTR2Dm1(zqsdi, qsdi, pvars)
      MKPTR2Dm1(zqsmi, qsmi, pvars)
      MKPTR2Dm1(zumfkfc, umfkfc, pvars)

      MKPTR2D(zmrk2, mrk2, pvars)

      if (advectke) then
         MKPTR2Dm1(zen, enplus, pvars)
      else
         MKPTR2Dm1(zen, en, pvars)
      endif
      
      call init2nan(l_en0, l_pw0)
      call init2nan(dummy1, dummy2, geop) 
      call init2nan(ppres, pudr, pddr, phsflx, dmsedt)
      call init2nan(pch1, pch1ten)
      call init2nan(qtl, qts)

      ! Prepare for deep CPS calls
      cdt1  = dt
      rcdt1 = 1./cdt1
      geop  = zgztherm*GRAV
      ier  = 0
      zcte = 0.
      zcqe = 0.
      zcqce= 0.
      pch1 = 0.
      pch1ten = 0.
      zkkfc = 0.

      ! Run selected deep convective scheme
      DEEP_CONVECTION: if (convec == 'SEC') then            ! ajustement convectif sec
         call secajus(zcte, ttp, sigma, psp, niter, 0.1, cdt1, ni, nkm1)
         if (phy_error_L) return
         call apply_tendencies(ttp, zcte, ztdmaskxdt, ni, nk, nkm1)

      else if (convec == 'KFC') then

         call kfcp8(ni,nkm1,zfcpflg,zkkfc,psp,ttp,qqp, &
              uu,vv,wz, &
              ztfcp,zhufcp,zufcp,zvfcp, &
              zprcten,zpriten, zqrkfc, &
              sigma,zdxdy,zrckfc,geop, &
              zabekfc, zareaup, zfdc, zdmfkfc, &
              zpeffkfc, zumfkfc, zzbasekfc, &
              zztopkfc, zwumaxkfc, &
              zqldi, zqsdi, &
              zrliq_int, zrice_int, &
              zkfcrf, zkfcsf, &
              kount,zdlat,zmg,zml,zwstar,critmask, delt)
         if (phy_error_L) return
         ztlc = zrckfc
         zqckfc = zprcten + zpriten

      else if (convec == 'KFC2') then

         call kfrpn4(ni,nkm1,zfcpflg,zkkfc,psp,ttp,qqp, &
              uu,vv,wz, &
              ztfcp,zhufcp,zufcp,zvfcp, &
              zufcp1,zvfcp1,zufcp2,zvfcp2, zufcp3,zvfcp3, zsufcp,zsvfcp, &
              zprcten,zpriten, zqrkfc, &
              sigma,zdxdy,zrckfc,geop, &
              zabekfc,zcapekfc,ztauckfc, zcinkfc, &
              zwklclplus, zwklcl, zareaup, zfdc, zdmfkfc, &
              zpeffkfc, zumfkfc, zzbasekfc, &
              zztopkfc, zwumaxkfc, &
              zqldi, zqsdi, &
              zrliq_int, zrice_int, &
              zkfcrf, zkfcsf, &
              kount,zdlat,zmg,zml,zwstar,ztstar,zen,zkt, &
              zcoadvu,zcoadvv,zcoage,zcowlcl,zcozlcl,zmrk2,critmask,delt)
         if (phy_error_L) return
         ztlc = zrckfc
         zqckfc = zprcten + zpriten

      else if (convec == 'BECHTOLD') then

         kcount(:) = int(zfcpflg(:))
         phsflx(:,:) = 0.
         do k = 1,nkm1
            do i = 1,ni
               ppres(i,k) = sigma(i,k)*psp(i)
            enddo
         enddo
         dummy1 = 0.
         dummy2 = 0.
         call bkf_deep4(ni, nkm1,  kidia, ni, kbdia, ktdia,      &
              dt, bkf_ldown, bkf_kice,                           &
              bkf_kens,                                          &
              ppres, zgztherm,zdxdy, phsflx,                     &
              ttp, qqp, dummy1, dummy2, uu, vv, wz,              &
              kcount, ztfcp, zhufcp, zprcten, zpriten,           &
              ztlc,ztsc,                                         &
              zumfkfc, zdmfkfc, zkfcrf, zkfcsf, zcapekfc, zztopkfc, zzbasekfc, &
              zqldi, zqsdi, zfdc,                          &
              zareaup, zwumaxkfc,                                &
              kfcmom, zufcp, zvfcp,                              &
              bkf_lch1conv, bkf_kch, pch1, pch1ten,              &
              pudr, pddr, zkkfc, psp, zdlat, zmg, zml)
         if (phy_error_L) return
         zqckfc = zprcten + zpriten
         if (kount == 0) kcount(:) = 0
         zfcpflg(:) = real(kcount(:))

      endif DEEP_CONVECTION
      if (phy_error_L) return

      ! Call post-deep shallow convective schemes
      IF_CONV_SHAL: if (conv_shal == 'BECHTOLD') then
         if (kount == 0) then
            ztpostshal = ttm
            zhupostshal = qqm
         endif
         do k = 1,nkm1
            do i = 1,ni
               ppres(i,k) = sigma(i,k)*psp(i)
               dmsedt(i,k) = (CPD*(ttp(i,k)-ztpostshal(i,k)) + CHLC*(qqp(i,k)-zhupostshal(i,k))) / dt
            enddo
         enddo
         dummy1 = 0.
         dummy2 = 0.

         ! Pre-shallow scheme state for budget
         if (pb_compute(zconesc, zconqsc, l_en0, l_pw0, &
              pvars, nkm1) /= PHY_OK) then
            call physeterror('cnv_main', &
                 'Problem computing preliminary budget for '//trim(conv_shal))
            return
         endif

         ! Set mask to avoid activating with deep
         cnv_active = 0.
         if (associated(zkkfc)) cnv_active = cnv_active + zkkfc

         ! Shallow convective scheme (mass-flux)
         call bkf_shallow6(ni, nkm1, dt, &
              ppres, zgztherm,                                     &
              ttp, qqp, dummy1, dummy2, uu, vv, wz, dmsedt,        &
              zfsc, zqlsc, zqssc,                      &
              ztshal, zhushal, zprctns,zpritns,ztusc, ztvsc,zumfs, &
              pch1, pch1ten,                &
              zkshal, psp, zwstar, zdxdy, zmrk2, cnv_active)
         if (phy_error_L) return

         ! Adjust tendencies to impose conservation
         if (pb_conserve(shal_conserve, ztshal, zhushal, pvars, &
              F_dqc=zprctns, F_dqi=zpritns) /= PHY_OK) then
            call physeterror('cnv_main', &
                 'Cannot correct conservation for '//trim(conv_shal))
            return
         endif
         
         ! Apply shallow convective tendencies
         if (bkf_lshalm) then
            call apply_tendencies(ttp,    qqp,     uu,    vv, &
                 &                ztshal, zhushal, ztusc, ztvsc, &
                 &                ztdmaskxdt, ni, nk, nkm1)
         else
            call apply_tendencies(ttp,    qqp, &
                 &                ztshal, zhushal, ztdmaskxdt, ni, nk, nkm1)
            
         endif

         ! Apply shallow convective tendencies for consdensed variables
         call conv_mp_tendencies1(zprctns, zpritns, ttp, qcp, ncp, qip, nip, &
              qti1p, nti1p, zti1p, ztdmaskxdt, ni, nk, nkm1)

         ! Post-scheme budget analysis
         if (pb_residual(zconesc, zconqsc, l_en0, l_pw0, pvars, &
              delt, nkm1) /= PHY_OK) then
            call physeterror('cnv_main', &
                 'Problem computing final budget for '//trim(conv_shal))
            return
         endif

         ! Store temperature and moisture profiles for MSE tendency calculation at next step
         ztpostshal = ttp
         zhupostshal = qqp

      endif IF_CONV_SHAL
      if (phy_error_L) return

      ! Pre-deep CPS state for budget
      if (pb_compute(zconedc, zconqdc, l_en0, l_pw0, &
           pvars, nkm1) /= PHY_OK) then
         call physeterror('cnv_main', &
              'Problem computing preliminary budget for '//trim(convec))
         return
      endif

      ! Translate deep convective tendency names if required
      if (associated(ztfcp)) zcte = zcte + ztfcp
      if (associated(zhufcp)) zcqe = zcqe + zhufcp
      if (associated(zqckfc)) zcqce = zcqce + zqckfc

      ! Adjust deep tendencies to impose conservation
      if (pb_conserve(deep_conserve, zcte, zcqe, pvars, &
           F_dqc=zcqce, F_dqi=zpriten, F_rain=ztlc, F_snow=ztsc) /= PHY_OK) then
         call physeterror('cnv_main', &
              'Cannot correct conservation for '//trim(convec))
         return
      endif

      ! Apply deep convective tendencies for non-condensed state variables
      if (kfcmom) then
         call apply_tendencies(ttp,  qqp,  uu,    vv, &
              &                zcte, zcqe, zufcp, zvfcp, &
              &                ztdmaskxdt, ni, nk, nkm1)
      else
         call apply_tendencies(ttp,  qqp, &
              &                zcte, zcqe, ztdmaskxdt, ni, nk, nkm1)         
      endif

      ! Apply deep convective tendencies for consdensed variables
      call conv_mp_tendencies1(zprcten, zpriten, ttp, qcp, ncp, qip, nip, &
            qti1p, nti1p, zti1p, ztdmaskxdt, ni, nk, nkm1)

      ! Post-deep budget analysis
      if (pb_residual(zconedc, zconqdc, l_en0, l_pw0, pvars, &
           delt, nkm1, F_rain=ztlc, F_snow=ztsc) /= PHY_OK) then
         call physeterror('cnv_main', &
              'Problem computing final budget for '//trim(convec))
         return
      endif

      ! Pre-midlevel CPS state for budget
      if (pb_compute(zconemc, zconqmc, l_en0, l_pw0, &
           pvars, nkm1) /= PHY_OK) then
         call physeterror('cnv_main', &
              'Problem computing preliminary budget for '//trim(conv_mid))
         return
      endif
      
      ! Represent midlevel elevated / low-CAPE convection
      MIDLEVEL_CONVECTION: if (conv_mid == 'KF') then

         cape = 0.
         if (associated(zcapekfc)) cape = zcapekfc
         cnv_active = 0.
         if (associated(zkkfc)) cnv_active = zkkfc
         cond_prflux = zrnflx + zsnoflx
         call kfmid1(ttp, qqp, uu, vv, wz, zgztherm, sigma, psp, zdxdy, zdlat, &
              cape, cnv_active, cond_prflux, &
              zmte, zmqe, zumid, zvmid, zprctnm, zpritnm, &
              zkmid, ztlcm, zfmc, zqlmi, zqsmi, zkfmrf, zkfmsf, zmcd, zmpeff, zmainc, &
              zmrk2, delt, ni, nkm1)
         if (phy_error_L) return
         zmqce = zprctnm + zpritnm

      endif MIDLEVEL_CONVECTION
      if (phy_error_L) return

      MIDLEVEL_CONS_TEND: if (conv_mid /= 'NIL') then
         
         ! Adjust midlevel tendencies to impose conservation
         if (pb_conserve(mid_conserve, zmte, zmqe, pvars, &
              F_dqc=zprctnm, F_dqi=zpritnm, F_rain=ztlcm) /= PHY_OK) then
            call physeterror('cnv_main', &
                 'Cannot correct conservation for '//trim(conv_mid))
            return
         endif

         ! Apply midlevel convective tendencies for non-condensed state variables
         if (associated(zumid) .and. associated(zvmid)) then
            call apply_tendencies(ttp,  qqp,  uu,    vv, &
                 &                zmte, zmqe, zumid, zvmid, &
                 &                ztdmaskxdt, ni, nk, nkm1)
         else
            call apply_tendencies(ttp,  qqp, &
                 &                zmte, zmqe, ztdmaskxdt, ni, nk, nkm1)
         endif

         ! Apply mid-level convective tendencies for condensed variables
         call conv_mp_tendencies1(zprctnm, zpritnm, ttp, qcp, ncp, qip, nip, &
              qti1p, nti1p, zti1p, ztdmaskxdt, ni, nk, nkm1)

         ! Post-midlevel budget analysis
         if (pb_residual(zconemc, zconqmc, l_en0, l_pw0, pvars, &
              delt, nkm1, F_rain=ztlcm) /= PHY_OK) then
            call physeterror('cnv_main', &
                 'Problem computing final budget for '//trim(conv_mid))
            return
         endif
         
      endif MIDLEVEL_CONS_TEND

      call msg_toall(MSG_DEBUG, 'cnv_main [END]')
      !----------------------------------------------------------------
      return
   end subroutine cnv_main4

end module cnv_main
