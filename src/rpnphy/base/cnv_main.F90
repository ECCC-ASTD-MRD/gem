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
!-------------------------------------- LICENCE END ----------------------------

module cnv_main
   implicit none
   private
   public :: cnv_main3

contains

   !/@*
   subroutine cnv_main3(d, dsiz, f, fsiz, v, vsiz, t0, q0, qc0, ilab, bett, &
        dt, ni, nk, kount)
      use, intrinsic :: iso_fortran_env, only: REAL64
      use debug_mod, only: init2nan
      use tdpack_const, only: CHLC, CPD, GRAV, RAUW
      use phy_status, only: phy_error_L
      use phybus
      use phy_options
      use cnv_options
      use kfrpn, only: kfrpn3
      use energy_budget, only: eb_en,eb_pw,eb_residual_en,eb_residual_pw,eb_conserve_en,eb_conserve_pw,EB_OK
      use shallconv, only: shallconv5
      use kfmid, only: kfmid1
      use tendency, only: apply_tendencies
      use conv_mp_tendencies, only: conv_mp_tendencies1
      use ens_perturb, only: ens_nc2d
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
      !@objective Interface to convection/condensation
      !@Arguments
      !          - Input -
      ! dsiz     dimension of dbus
      ! fsiz     dimension of fbus
      ! vsiz     dimension of vbus
      ! dt       timestep (parameter) (sec.)
      ! ni       horizontal running length
      ! nk       vertical dimension
      ! kount    timestep number
      !
      !          - Input/Output -
      ! d        dynamics input field
      ! f        historic variables for the physics
      ! v        physics tendencies and other output fields from the physics
      !
      !          - Output -
      ! ilab     flag array: an indication of convective activity from Kuo schemes
      ! bett     estimated averaged cloud fraction growth rate for kuostd

      integer, intent(in) :: fsiz,vsiz,dsiz,ni,nk,kount
      integer,dimension(ni,nk-1), intent(out) :: ilab
      real,   dimension(ni), intent(out)      :: bett
      real,   dimension(ni,nk-1), intent(inout) :: t0,q0,qc0
      real,   target, intent(inout) :: f(fsiz), v(vsiz), d(dsiz)
      real dt
      !@Author L.Spacek, November 2011
      !@Revisions
      ! 001   -PV/JM-nov2014- fix communication between deep convection and MY_DM
      ! 002   -PV-dec2014 -
      ! 003   -JY-PV-Jul2015- Separate deep and shallow convection scheme in unified bechtold scheme
      ! 004   -JY-PV-Nov2016- Bkf-shallow only checked if no deep convection is active
      !*@/
#include <msg.h>
#include "phymkptr.hf"

      include "surface.cdk"
      include "comphy.cdk"

      real    :: cdt1, rcdt1
      integer :: i,k,niter, ier, nkm1, ier2

      real(REAL64), dimension(ni) :: l_en0, l_en, l_pw0, l_pw, l_enr, l_pwr
      real, dimension(ni,nk-1) :: zfm, dotp
      !
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
      real,    dimension(ni,nk-1)      :: dummy1, dummy2, geop, cond_prflux
      real,    dimension(ni,nk-1)      :: ppres,pudr, pddr, phsflx, purv, dmsedt
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

      real, pointer, dimension(:) :: psm, psp, ztdmask, zfcpflg, zkkfc, zrckfc, &
           ztlc, ztsc, zkshal, zwstar, zconedc, zconesc, zconqdc, zconqsc, ztlcs, ztscs, &
           ztstar, zcapekfc, ztauckfc, zcinkfc, zmg, zml, zdlat, zdxdy, zkmid, &
           zabekfc, zpeffkfc, zrice_int, zrliq_int, zwumaxkfc, zzbasekfc, zztopkfc, &
           zcoadvu,zcoadvv,zcoage,zcowlcl,zwklcl,zcozlcl,ztlcm, zconemc, zconqmc, &
           zmcd,zmpeff,zmainc
      real, pointer, dimension(:,:) :: ncp, nip, qcm, qcp, qip, qrp, qqm, qqp, &
           sigma, ttm, ttp, uu, vv, wz, zfdc, zgztherm, zhufcp, zhushal, &
           zprcten, zpriten, zqckfc, ztfcp, ztshal, ztusc, ztvsc, zufcp, zvfcp, &
           zprctns,zpritns,qti1p, nti1p, zfsc, zqlsc, zqssc, zfbl, &
           zumfs, ztpostshal, zhupostshal, zcqce, zcqe, zcte, zen, zkt, &
           zareaup, zdmfkfc, zkfcrf, zkfcsf, zqldi, zqrkfc, zqsdi, zumfkfc, &
           zqcz, zqdifv, zwklclplus, zkfmrf, zkfmsf, zqlmi, zqsmi, zfmc, zprctnm, zpritnm, &
           zmqce, zmte, zmqe, zumid, zvmid, zrnflx, zsnoflx, zmrk2

      !----------------------------------------------------------------
      call msg_toall(MSG_DEBUG, 'cnv_main [BEGIN]')

      nkm1 = nk - 1

      MKPTR1D(psm, pmoins, f)
      MKPTR1D(psp, pmoins, f)
      MKPTR1D(zabekfc, abekfc, f)
      MKPTR1D(zcapekfc, capekfc, f)
      MKPTR1D(ztauckfc, tauckfc, f)
      MKPTR1D(zcinkfc, cinkfc, f)
      MKPTR1D(zcoadvu, coadvu, f)
      MKPTR1D(zcoadvv, coadvv, f)
      MKPTR1D(zcoage, coage, f)
      MKPTR1D(zconedc, conedc, v)
      MKPTR1D(zconemc, conemc, v)
      MKPTR1D(zconesc, conesc, v)
      MKPTR1D(zconqdc, conqdc, v)
      MKPTR1D(zconqmc, conqmc, v)
      MKPTR1D(zconqsc, conqsc, v)
      MKPTR1D(zcowlcl, cowlcl, f)
      MKPTR1D(zcozlcl, cozlcl, f)
      MKPTR1D(zdlat, dlat, f)
      MKPTR1D(zdxdy, dxdy, f)
      MKPTR1D(zfcpflg, fcpflg, f)
      MKPTR1D(zkkfc, kkfc, v)
      MKPTR1D(zkmid, kmid, v)
      MKPTR1D(zkshal, kshal, v)
      MKPTR1D(zmainc, mainc, v)
      MKPTR1D(zmcd, mcd, v)
      MKPTR1D(zmg, mg, f)
      MKPTR1D(zml, ml, f)
      MKPTR1D(zmpeff, mpeff, v)
      MKPTR1D(zpeffkfc, peffkfc, f)
      MKPTR1D(zrckfc, rckfc, f)
      MKPTR1D(ztlcm, tlcm, v)
      MKPTR1D(zrice_int, rice_int, f)
      MKPTR1D(zrliq_int, rliq_int, f)
      MKPTR1D(ztdmask, tdmask, f)
      MKPTR1D(ztlc, tlc, f)
      MKPTR1D(ztlcs, tlcs, v)
      MKPTR1D(ztsc, tsc, f)
      MKPTR1D(ztscs, tscs, v)
      MKPTR1D(ztstar, tstar, v)
      MKPTR1D(zwklcl, wklcl, f)
      MKPTR1D(zwstar, wstar, v)
      MKPTR1D(zwumaxkfc, wumaxkfc, f)
      MKPTR1D(zzbasekfc, zbasekfc, f)
      MKPTR1D(zztopkfc, ztopkfc, f)

      MKPTR2D(sigma, sigw, d)

      MKPTR2Dm1(nti1p, nti1plus, d)
      MKPTR2Dm1(qti1p, qti1plus, d)
      MKPTR2Dm1(ncp, ncplus, d)
      MKPTR2Dm1(nip, niplus, d)
      MKPTR2Dm1(qcm, qcmoins, d)
      MKPTR2Dm1(qcp, qcplus, d)
      MKPTR2Dm1(qip, qiplus, d)
      MKPTR2Dm1(qqm, humoins, d)
      MKPTR2Dm1(qqp, huplus, d)
      MKPTR2Dm1(qrp, qrplus, d)
      MKPTR2Dm1(ttm, tmoins, d)
      MKPTR2Dm1(ttp, tplus, d)
      MKPTR2Dm1(uu, uplus, d)
      MKPTR2Dm1(vv, vplus, d)
      MKPTR2Dm1(wz, wplus, d)
      MKPTR2Dm1(zcqce, cqce, v)
      MKPTR2Dm1(zcqe, cqe, v)
      MKPTR2Dm1(zcte, cte, v)
      MKPTR2Dm1(zen, en, f)
      MKPTR2Dm1(zfbl, fbl, f)
      MKPTR2Dm1(zfdc, fdc, f)
      MKPTR2Dm1(zfmc, fmc, f)
      MKPTR2Dm1(zfsc, fsc, f)
      MKPTR2Dm1(zgztherm, gztherm, v)
      MKPTR2Dm1(zhufcp, hufcp, f)
      MKPTR2Dm1(zmqe, mqe, v)
      MKPTR2Dm1(zhupostshal, hupostshal, f)
      MKPTR2Dm1(zhushal, hushal, v)
      MKPTR2Dm1(zkt, kt, v)
      MKPTR2Dm1(zprcten, prcten, f)
      MKPTR2Dm1(zprctnm, prctnm, v)
      MKPTR2Dm1(zprctns, prctns, v)
      MKPTR2Dm1(zpriten, priten, f)
      MKPTR2Dm1(zpritnm, pritnm, v)
      MKPTR2Dm1(zpritns, pritns, v)
      MKPTR2Dm1(zqckfc, qckfc, f)
      MKPTR2Dm1(zmqce, mqce, v)
      MKPTR2Dm1(zqcz, qcz, v)
      MKPTR2Dm1(zqdifv, qdifv, v)
      MKPTR2Dm1(zqlsc, qlsc, v)
      MKPTR2Dm1(zqssc, qssc, v)
      MKPTR2Dm1(ztfcp, tfcp, f)
      MKPTR2Dm1(zmte, mte, v)
      MKPTR2Dm1(zrnflx, rnflx, f)
      MKPTR2Dm1(zsnoflx, snoflx, f)
      MKPTR2Dm1(ztpostshal, tpostshal, f)
      MKPTR2Dm1(ztshal, tshal, v)
      MKPTR2Dm1(ztusc, tusc, v)
      MKPTR2Dm1(ztvsc, tvsc, v)
      MKPTR2Dm1(zufcp, ufcp, f)
      MKPTR2Dm1(zumfs, umfs, v)
      MKPTR2Dm1(zumid, umid, v)
      MKPTR2Dm1(zvfcp, vfcp, f)
      MKPTR2Dm1(zwklclplus, wklclplus, d)
      MKPTR2Dm1(zvmid, vmid, v)
      MKPTR2Dm1(zareaup, areaup, f)
      MKPTR2Dm1(zdmfkfc, dmfkfc, f)
      MKPTR2Dm1(zkfcrf, kfcrf, f)
      MKPTR2Dm1(zkfcsf, kfcsf, f)
      MKPTR2Dm1(zkfmrf, kfmrf, v)
      MKPTR2Dm1(zkfmsf, kfmsf, v)
      MKPTR2Dm1(zqldi, qldi, f)
      MKPTR2Dm1(zqlmi, qlmi, v)
      MKPTR2Dm1(zqrkfc, qrkfc, f)
      MKPTR2Dm1(zqsdi, qsdi, f)
      MKPTR2Dm1(zqsmi, qsmi, v)
      MKPTR2Dm1(zumfkfc, umfkfc, f)

      MKPTR2DN(zmrk2, mrk2, ni, ens_nc2d, f)

      call init2nan(l_en0, l_en, l_pw0, l_pw, l_enr, l_pwr)
      call init2nan(zfm, dotp, dummy1, dummy2, geop)
      call init2nan(ppres, pudr, pddr, phsflx, purv, dmsedt)
      call init2nan(pch1, pch1ten)

      ! Shallow convective schemes called before deep (post-deep schemes called later in this routine)
      IF_KTRSNT: if (any(conv_shal == (/'KTRSNT_MG', 'KTRSNT   '/))) then

         ! Pre-scheme state for energy budget
         if (associated(zconesc)) then
            ier  = eb_en(l_en0,ttp,qqp,qcp,sigma,psp,nkm1)
            ier2 = eb_pw(l_pw0,qqp,qcp,sigma,psp,nkm1)
            if (ier /= EB_OK .or. ier2 /= EB_OK) then
               call physeterror('cnv_main', 'Problem computing preliminary energy budget for '//trim(conv_shal))
               return
            endif
         endif

         ! Shallow convective scheme (Kuo transient types)
         call shallconv5(ttp, qqm, sigma, zgztherm, psm, &
              ztshal, zhushal, ztlcs, ztscs, zqlsc, zqssc, zkshal, &
              zfsc, zqcz, zqdifv, dt, ni, nk, nkm1)
         if (phy_error_L) return

         ! Impose conservation by adjustment of moisture and temperature tendencies
         SHAL_EARLY_TENDENCY_ADJUSTMENT: if (shal_conserve == 'TEND') then

            ! Apply humidity tendency correction for total water conservation
            ier = eb_conserve_pw(zhushal,zhushal,ttp,qqp,sigma,psp,nkm1,F_rain=ztlcs*RAUW,F_snow=ztscs*RAUW)
            ! Apply temperature tendency correction for liquid water static energy conservation
            ier2 = eb_conserve_en(ztshal,ztshal,zhushal,ttp,qqp,qcp,sigma,psp,nkm1,F_rain=ztlcs*RAUW,F_snow=ztscs*RAUW)
            if (ier /= EB_OK .or. ier2 /= EB_OK) then
               call physeterror('cnv_main', 'Problem correcting for liquid water static energy conservation for '//trim(conv_shal))
               return
            endif

         endif SHAL_EARLY_TENDENCY_ADJUSTMENT

         ! Apply shallow convective tendencies to state variables
         call apply_tendencies(ttp, ztshal, ztdmask, ni, nk, nkm1)
         call apply_tendencies(qqp, zhushal, ztdmask, ni, nk, nkm1)

         ! Post-scheme energy budget analysis
         if (associated(zconesc)) then
            ! Compute post-scheme state
            ier  = eb_en(l_en,ttp,qqp,qcp,sigma,psp,nkm1)
            ier2 = eb_pw(l_pw,qqp,qcp,sigma,psp,nkm1)
            if (ier == EB_OK .and. ier2 == EB_OK) then
               ! Compute residuals
               ier  = eb_residual_en(l_enr,l_en0,l_en,ttp,qqp,qcp,delt,nkm1,F_rain=ztlcs*RAUW,F_snow=ztscs*RAUW)
               ier2 = eb_residual_pw(l_pwr,l_pw0,l_pw,ttp,delt,nkm1,F_rain=ztlcs*RAUW,F_snow=ztscs*RAUW)
            endif
            if (ier /= EB_OK .or. ier2 /= EB_OK) then
               call physeterror('cnv_main', 'Problem computing final energy budget for '//trim(conv_shal))
               return
            endif
            zconesc(:) = real(l_enr)
            zconqsc(:) = real(l_pwr)
         endif

      endif IF_KTRSNT

      ! Prepare for deep CPS calls
      cdt1  = dt
      rcdt1 = 1./cdt1
      geop  = zgztherm*GRAV

      ier  = 0
      ilab = 0
      zcte = 0.
      zcqe = 0.
      zcqce= 0.
      qc0  = 0.
      bett = 0.
      pch1 = 0.
      pch1ten = 0.

      t0 = ttp
      q0 = qqp
      where(qqp < 0.) qqp = 0.

      if (stcond /= 'NIL' ) then
         qc0 = qcp
         where(qcp < 0.) qcp = 0.
         where(qcm < 0.) qcm = 0.
      endif

      ! Run selected deep convective scheme
      DEEP_CONVECTION: if (convec == 'SEC') then            ! ajustement convectif sec
         call secajus(zcte, ttp, sigma, psp, niter, 0.1, cdt1, ni, nkm1)
         if (phy_error_L) return
         call apply_tendencies(ttp, zcte, ztdmask, ni, nk, nkm1)

      else if (convec == 'OLDKUO') then

         zfm(:,:) = qcp(:,:)
         call omega_from_w(dotp,wz,ttp,qqp,psp,sigma,ni,nkm1)
         call kuo6(zcte,zcqe,ztlc,ztsc, &
              ilab,zfdc,zfbl,dotp,zfm, &
              ttp,qqp,qqm, &
              geop,psp,psm, &
              sigma, cdt1, ni, ni, nkm1, &
              satuco, ccctim, cmelt)
         if (phy_error_L) return
         if (stcond /= 'NIL' ) &
              zcqce = (zfm -qcp)*rcdt1

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

      else if (any(convec == (/'KFC2', 'KFC3'/))) then

         call kfrpn3(ni,nkm1,zfcpflg,zkkfc,psp,ttp,qqp, &
              uu,vv,wz, &
              ztfcp,zhufcp,zufcp,zvfcp, &
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

      else if (convec == 'KUOSTD') then

         call omega_from_w(dotp,wz,ttp,qqp,psp,sigma,ni,nkm1)
         call lsctrol(ilab, dotp, sigma, ni, nkm1)
         call kuostd2(zcte,zcqe,ilab,zfdc,bett, &
              ttp,ttm,qqp,qqm,geop,psp, &
              sigma, cdt1, ni, nkm1)
         if (phy_error_L) return

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

         ! Pre-scheme state for energy budget
         if (associated(zconesc)) then
            ier  = eb_en(l_en0,ttp,qqp,qcp,sigma,psp,nkm1)
            ier2 = eb_pw(l_pw0,qqp,qcp,sigma,psp,nkm1)
            if (ier /= EB_OK .or. ier2 /= EB_OK) then
               call physeterror('cnv_main', 'Problem computing preliminary energy budget for '//trim(conv_shal))
               return
            endif
         endif

         ! Shallow convective scheme (Bechtold type)
         cnv_active = 0.  !shallow does not activate where deep or mid-level are active
         if (associated(zkkfc)) cnv_active = cnv_active + zkkfc
!!$         if (associated(zkmid)) cnv_active = cnv_active + zkmid

         call bkf_shallow6(ni, nkm1, dt, &
              ppres, zgztherm,                                     &
              ttp, qqp, dummy1, dummy2, uu, vv, wz, dmsedt,        &
              zfsc, zqlsc, zqssc,                      &
              ztshal, zhushal, zprctns,zpritns,ztusc, ztvsc,zumfs, &
              pch1, pch1ten,                &
              zkshal, psp, zwstar, zdxdy, zmrk2, cnv_active)
         if (phy_error_L) return

         ! Impose conservation by adjustment of moisture and temperature tendencies
         SHAL_LATE_TENDENCY_ADJUSTMENT: if (shal_conserve == 'TEND') then

            ! Apply humidity tendency correction for total water conservation
            ier = eb_conserve_pw(zhushal,zhushal,ttp,qqp,sigma,psp,nkm1,F_dqc=zprctns,F_dqi=zpritns)
            ! Apply temperature tendency correction for liquid water static energy conservation
            ier2 = eb_conserve_en(ztshal,ztshal,zhushal,ttp,qqp,qcp,sigma,psp,nkm1,F_dqc=zprctns,F_dqi=zpritns)
            if (ier /= EB_OK .or. ier2 /= EB_OK) then
               call physeterror('cnv_main', 'Problem correcting for liquid water static energy conservation for '//trim(conv_shal))
               return
            endif

         endif SHAL_LATE_TENDENCY_ADJUSTMENT

         ! Apply shallow convective tendencies
         call apply_tendencies(ttp, ztshal, ztdmask, ni, nk, nkm1)
         call apply_tendencies(qqp, zhushal, ztdmask, ni, nk, nkm1)
         call apply_tendencies(qcp, zprctns, ztdmask, ni, nk, nkm1)
         call apply_tendencies(qcp, zpritns, ztdmask, ni, nk, nkm1)
         if (bkf_lshalm) then
            call apply_tendencies(uu, ztusc, ztdmask, ni, nk, nkm1)
            call apply_tendencies(vv, ztvsc, ztdmask, ni, nk, nkm1)
         endif

         ! Post-scheme energy budget analysis
         if (associated(zconesc)) then
            ! Compute post-scheme state
            ier  = eb_en(l_en,ttp,qqp,qcp,sigma,psp,nkm1)
            ier2 = eb_pw(l_pw,qqp,qcp,sigma,psp,nkm1)
            if (ier == EB_OK .and. ier2 == EB_OK) then
               ! Compute residuals
               ier  = eb_residual_en(l_enr,l_en0,l_en,ttp,qqp,qcp,delt,nkm1)
               ier2 = eb_residual_pw(l_pwr,l_pw0,l_pw,ttp,delt,nkm1)
            endif
            if (ier /= EB_OK .or. ier2 /= EB_OK) then
               call physeterror('cnv_main', 'Problem computing final energy budget for '//trim(conv_shal))
               return
            endif
            zconesc(:) = real(l_enr)
            zconqsc(:) = real(l_pwr)
         endif

         ! Store temperature and moisture profiles for MSE tendency calculation at next step
         ztpostshal = ttp
         zhupostshal = qqp

      endif IF_CONV_SHAL
      if (phy_error_L) return

      ! Pre-deep CPS state for energy budget
      if (associated(zconedc)) then
         ier  = eb_en(l_en0,ttp,qqp,qcp,sigma,psp,nkm1)
         ier2 = eb_pw(l_pw0,qqp,qcp,sigma,psp,nkm1)
         if (ier /= EB_OK .or. ier2 /= EB_OK) then
            call physeterror('cnv_main', 'Problem computing preliminary energy budget for '//trim(convec))
            return
         endif
      endif

      ! Translate deep convective tendency names if required
      if (associated(ztfcp)) zcte = zcte + ztfcp
      if (associated(zhufcp)) zcqe = zcqe + zhufcp
      if (associated(zqckfc)) zcqce = zcqce + zqckfc

      ! Impose conservation by adjustment of moisture and temperature tendencies
      DEEP_TENDENCY_ADJUSTMENT: if (deep_conserve == 'TEND') then

         ! Apply humidity tendency correction for total water conservation
         ier  = eb_conserve_pw(zcqe,zcqe,ttp,qqp,sigma,psp,nkm1,F_dqc=zcqce,F_rain=ztlc,F_snow=ztsc)
         ! Apply temperature tendency correction for liquid water static energy conservation
         ier2 = eb_conserve_en(zcte,zcte,zcqe,ttp,qqp,qcp,sigma,psp,nkm1,F_dqc=zcqce,F_rain=ztlc,F_snow=ztsc)
         if (ier /= EB_OK .or. ier2 /= EB_OK) then
            call physeterror('cnv_main', 'Problem correcting for liquid water static energy conservation for '//trim(convec))
            return
         endif

      endif DEEP_TENDENCY_ADJUSTMENT

      ! Apply deep convective tendencies for non-condensed state variables
      call apply_tendencies(ttp, zcte, ztdmask, ni, nk, nkm1)
      call apply_tendencies(qqp, zcqe, ztdmask, ni, nk, nkm1)
      if (kfcmom) then
         call apply_tendencies(uu, zufcp, ztdmask, ni, nk, nkm1)
         call apply_tendencies(vv, zvfcp, ztdmask, ni, nk, nkm1)
      endif

      ! Apply deep convective tendencies for consdensed variables
      call conv_mp_tendencies1(zprcten, zpriten, ttp, qcp, ncp, qip, nip, qti1p, nti1p, ztdmask, ni, nk, nkm1)

      ! Post-deep CPS energy budget analysis
      if (associated(zconedc)) then
         ! Compute post-scheme state
         ier  = eb_en(l_en,ttp,qqp,qcp,sigma,psp,nkm1)
         ier2 = eb_pw(l_pw,qqp,qcp,sigma,psp,nkm1)
         if (ier == EB_OK .and. ier2 == EB_OK) then
            ! Compute residuals
            ier  = eb_residual_en(l_enr,l_en0,l_en,ttp,qqp,qcp,delt,nkm1,F_rain=ztlc,F_snow=ztsc)
            ier2 = eb_residual_pw(l_pwr,l_pw0,l_pw,ttp,delt,nkm1,F_rain=ztlc,F_snow=ztsc)
         endif
         if (ier /= EB_OK .or. ier2 /= EB_OK) then
             call physeterror('cnv_main', 'Problem computing final energy budget for '//trim(convec))
            return
         endif
         zconedc(:) = real(l_enr)
         zconqdc(:) = real(l_pwr)
      endif

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

      ! Pre-mid-level convection state for energy budget
      if (associated(zconemc)) then
         ier  = eb_en(l_en0,ttp,qqp,qcp,sigma,psp,nkm1)
         ier2 = eb_pw(l_pw0,qqp,qcp,sigma,psp,nkm1)
         if (ier /= EB_OK .or. ier2 /= EB_OK) then
            call physeterror('cnv_main', 'Problem computing preliminary energy budget for '//trim(conv_mid))
            return
         endif
      endif

      ! Impose conservation by adjustment of moisture and temperature tendencies
      MID_TENDENCY_ADJUSTMENT: if (mid_conserve == 'TEND') then

         ! Apply humidity tendency correction for total water conservation
         ier  = eb_conserve_pw(zmqe,zmqe,ttp,qqp,sigma,psp,nkm1,F_dqc=zmqce,F_rain=ztlcm)
         ! Apply temperature tendency correction for liquid water static energy conservation
         ier2 = eb_conserve_en(zmte,zmte,zmqe,ttp,qqp,qcp,sigma,psp,nkm1,F_dqc=zmqce,F_rain=ztlcm)
         if (ier /= EB_OK .or. ier2 /= EB_OK) then
            call physeterror('cnv_main', 'Problem correcting for liquid water static energy conservation for '//trim(conv_mid))
            return
         endif

      endif MID_TENDENCY_ADJUSTMENT

      ! Apply mid-level convective tendencies for non-condensed state variables
      call apply_tendencies(ttp, zmte, ztdmask, ni, nk, nkm1)
      call apply_tendencies(qqp, zmqe, ztdmask, ni, nk, nkm1)
      if (associated(zumid)) call apply_tendencies(uu, zumid, ztdmask, ni, nk, nkm1)
      if (associated(zvmid)) call apply_tendencies(vv, zvmid, ztdmask, ni, nk, nkm1)

      ! Apply mid-level convective tendencies for condensed variables
      call conv_mp_tendencies1(zprctnm, zpritnm, ttp, qcp, ncp, qip, nip, qti1p, nti1p, ztdmask, ni, nk, nkm1)

      ! Post-mid-level convection energy budget analysis
      if (associated(zconemc)) then
         ! Compute post-scheme state
         ier  = eb_en(l_en,ttp,qqp,qcp,sigma,psp,nkm1)
         ier2 = eb_pw(l_pw,qqp,qcp,sigma,psp,nkm1)
         if (ier == EB_OK .and. ier2 == EB_OK) then
            ! Compute residuals
            ier  = eb_residual_en(l_enr,l_en0,l_en,ttp,qqp,qcp,delt,nkm1,F_rain=ztlcm)
            ier2 = eb_residual_pw(l_pwr,l_pw0,l_pw,ttp,delt,nkm1,F_rain=ztlcm)
         endif
         if (ier /= EB_OK .or. ier2 /= EB_OK) then
             call physeterror('cnv_main', 'Problem computing final energy budget for '//trim(conv_mid))
            return
         endif
         zconemc(:) = real(l_enr)
         zconqmc(:) = real(l_pwr)
      endif

      call msg_toall(MSG_DEBUG, 'cnv_main [END]')
      !----------------------------------------------------------------
      return
   end subroutine cnv_main3

end module cnv_main
