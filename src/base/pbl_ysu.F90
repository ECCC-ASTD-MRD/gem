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

module pbl_ysu
  implicit none
  private
  public :: pbl_ysu1

contains

  !/@*
  function pbl_ysu1(dbus, fbus, vbus, dt, trnch, ni, nk, nkm1) result(status)
    use, intrinsic :: iso_fortran_env, only: INT64
    use debug_mod, only: init2nan
    use tdpack_const, only: CPD, DELTA, GRAV, KARMAN, CAPPA, RGASD, EPS1, &
         CHLC, RGASV
    use phy_options
    use phybus
    use sfclayer, only: sl_prelim, sl_sfclayer, SL_OK
    use module_bl_ysu, only: ysu2d
    implicit none
#include <rmnlib_basics.hf>
    !@Object Wrapper for Yonsei University (YSU) PBL scheme
    !@Arguments
    real, pointer, dimension(:), contiguous :: dbus, fbus, vbus !Bus data
    real, intent(in) :: dt                                      !Time step (sec)
    integer, intent(in) :: trnch                                !Slice number (tile j)
    integer, intent(in) :: ni, nk, nkm1                         !Slab dimensions
    integer :: status                                           !Return status (RMN_OK or RMN_ERR)

    !@Author Ron McTaggart-Cowan (Summer 2022)
    !*@/
#include "phymkptr.hf"
#include "surface.cdk"

    ! Internal parameters
    integer, parameter :: NDIFF=3
    integer, parameter :: MAX_DT=30.                            !Maximum step (s)

    ! Internal variables
    integer :: k, ki, istat, step, nstep
    integer, dimension(ni) :: kpbl
    real :: stepdt, idt
#ifdef __INTEL19_COMPILER
    real, pointer, dimension(:) :: zz0, zqsurf, zz0t, zbt, zfc, zfv
#else
    real, pointer, dimension(:), contiguous :: zz0, zqsurf, zz0t, zbt, zfc, zfv
#endif
    real, pointer, dimension(:), contiguous :: zp0, zmg, &
         zwstar, zudiag, zvdiag, zzusl, zztsl, ztsrad, &
         zdlat, zfcor, zalfat, zalfaq, zbm
    real, dimension(ni) :: dusfc, dvsfc, dtsfc, dqsfc, ldelta, zero, wspd, &
         lmg, lfsh, vmod, vdir, thair, ribsl, fmsl, fhsl, &
         lfh, lfv, rho, lfrv, lfrv0, hpbl
    real, pointer, dimension(:,:), contiguous :: zhuplus, zqcplus, ztplus, &
         zuplus, zvplus, zudifv, zvdifv, ztdifv, zqdifv, zqcdifv, &
         zsigm, zsigt, zgztherm, zkm, zkt, ztrad, zgzmom
    real, dimension(ni,nkm1) :: iuplust, uplust, ivplust, vplust, itplus, &
         iqplus, iqcplus, iudifv, udifvt, ivdifv, vdifvt, itdifv, &
         iqdifv, iqcdifv, thdifv, exner, iexner, &
         iprest, prest, presm, dzt, idzt, ikm, ikt, thrad, ithrad, &
         ustept, vstept, ustepm, vstepm, exnerrpn
    real, dimension(ni,nk) :: ipresm, ipresmo
    real, dimension(ni,nkm1*3) :: imoist, imoisttnd
    real, dimension(:,:,:), pointer, contiguous :: zvcoef

    ! Basic initializations
    status = RMN_ERR
    call init2nan(dusfc, dvsfc, dtsfc, dqsfc, ldelta, zero, wspd)
    call init2nan(lmg, lfsh, vmod, vdir, thair, ribsl, fmsl, fhsl)
    call init2nan(lfh, lfv, rho, lfrv, lfrv0, hpbl)
    call init2nan(iuplust, uplust, ivplust, vplust, itplus)
    call init2nan(iqplus, iqcplus, iudifv, udifvt, ivdifv, vdifvt, itdifv)
    call init2nan(iqdifv, iqcdifv, thdifv, exner, iexner)
    call init2nan(iprest, prest, presm, dzt, idzt, ikm, ikt, thrad, ithrad)
    call init2nan(ustept, vstept, ustepm, vstepm, exnerrpn)
    call init2nan(ipresm, ipresmo)

    ! Estract bus pointers
    MKPTR2Dm1(zhuplus, huplus, dbus)
    MKPTR1D(zp0, p0_plus, dbus)
    MKPTR2Dm1(zqcplus, qcplus, dbus)
    MKPTR2Dm1(zsigm, sigm, dbus)
    MKPTR2Dm1(zsigt, sigt, dbus)
    MKPTR2Dm1(ztplus, tplus, dbus)
    MKPTR2Dm1(zuplus, uplus, dbus)
    MKPTR2Dm1(zvplus, vplus, dbus)
    MKPTR1D(zdlat, dlat, fbus)
    MKPTR1D(zmg, mg, fbus)
    MKPTR1DK(zqsurf, qsurf, indx_agrege, fbus)
    MKPTR1D(ztsrad, tsrad, fbus)
    MKPTR1D(zudiag, udiag, fbus)
    MKPTR1D(zvdiag, vdiag, fbus)
    MKPTR1DK(zz0, z0, indx_agrege, fbus)
    MKPTR1DK(zz0t, z0t, indx_agrege, fbus)
    MKPTR1D(zalfaq, alfaq, vbus)
    MKPTR1D(zalfat, alfat, vbus)
    MKPTR1D(zbm, bm, vbus)
    MKPTR1DK(zbt, bt, indx_agrege, vbus)
    MKPTR1DK(zfc, fc, indx_agrege, vbus)
    MKPTR1DK(zfv, fv, indx_agrege, vbus)
    MKPTR1D(zfcor, fcor, vbus)
    MKPTR2Dm1(zgzmom, gzmom, vbus)
    MKPTR2Dm1(zgztherm, gztherm, vbus)
    MKPTR2D(zkm, km, vbus)
    MKPTR2D(zkt, kt, vbus)
    MKPTR2Dm1(zqcdifv, qcdifv, vbus)
    MKPTR2Dm1(zqdifv, qdifv, vbus)
    MKPTR2Dm1(ztdifv, tdifv, vbus)
    MKPTR2Dm1(ztrad, trad, vbus)
    MKPTR2Dm1(zudifv, udifv, vbus)
    MKPTR3D(zvcoef, vcoef, 2, vbus)
    MKPTR2Dm1(zvdifv, vdifv, vbus)
    MKPTR1D(zwstar, wstar, vbus)
    MKPTR1D(zztsl, ztsl, vbus)
    MKPTR1D(zzusl, zusl, vbus)

    ! Compute pressure-based cubes
    do k=1,nkm1
       presm(:,k) = zsigm(:,k) * zp0(:)
       prest(:,k) = zsigt(:,k) * zp0(:)
    enddo
    exner(:,:) = (prest(:,:)/1.e5)**CAPPA
    exnerrpn(:,:) = zsigt(:,:)**CAPPA

    ! Interpolate momentum-level variables to thermodynamic levels
    call vint_mom2thermo(uplust, zuplus, zvcoef, ni, nkm1)
    call vint_mom2thermo(vplust, zvplus, zvcoef, ni, nkm1)

    ! Diagnose lowest thermo-level winds from surface layer
    istat = sl_prelim(ztplus(:,nkm1), zhuplus(:,nkm1), zuplus(:,nkm1), &
         zvplus(:,nkm1), zp0, zzusl, spd_air=vmod, dir_air=vdir, &
         min_wind_speed=VAMIN)
    if (istat /= SL_OK) then
       call physeterror('pbl_ysu::pbl_ysu1', 'problem in sl_prelim()')
       return
    endif
    thair(:) = ztplus(:,nkm1) / exnerrpn(:,nkm1)
    istat = sl_sfclayer(thair, zhuplus(:,nkm1), vmod, vdir, zzusl, zztsl, &
         ztsrad, zqsurf, zz0, zz0t, zdlat, zfcor, hghtm_diag_row=zgztherm(:,nkm1), &
         u_diag=uplust(:,nkm1), v_diag=vplust(:,nkm1), &
         rib=ribsl, stabm=fmsl, stabt=fhsl)
    if (istat /= SL_OK) then
       call physeterror('pbl_ysu::pbl_ysu1', 'problem in sl_sfclayer()')
       return
    endif

    ! Compute static derived quantities
    do k=1,nkm1-1
       dzt(:,k) = zgzmom(:,k) - zgzmom(:,k+1)
    enddo
    dzt(:,nkm1) = zgzmom(:,nkm1)
    wspd(:) = sqrt(uplust(:,nkm1)**2 + vplust(:,nkm1)**2)
    thrad(:,:) = ztrad(:,:)/exner(:,:)
    zero(:) = 0.
    lmg(:) = 1
    where (zmg(:) .lt. critmask) lmg(:) = 2
    
    ! Flip all arrays for input
    ipresm(:,1) = zp0(:)
    do k=1,nkm1
       ki = nkm1 - (k-1)
       iuplust(:,ki) = uplust(:,k)
       ivplust(:,ki) = vplust(:,k)
       itplus(:,ki) = ztplus(:,k)
       iqplus(:,ki) = zhuplus(:,k)
       iexner(:,ki) = exner(:,k)
       ipresm(:,ki+1) = presm(:,k)
       iprest(:,ki) = prest(:,k)
       idzt(:,ki) = dzt(:,k)
       ithrad(:,ki) = thrad(:,k)
    enddo
    if (associated(zqcplus)) then
       do k=1,nkm1
          ki = nkm1 - (k-1)
          iqcplus(:,ki) = zqcplus(:,k)
       enddo
    else
       iqcplus(:,:) = 0.
    endif
    ipresmo(:,:) = ipresm(:,:)
    ikm(:,:) = 0.
    ikt(:,:) = 0.
    
    ! Time-splitting for stability
    if (pbl_ysu_rpnsolve) then
       nstep = 1
    else
       nstep = max(ceiling(dt / MAX_DT), 1)
    endif
    stepdt = dt / nstep
    TIME_STEPS: do step=1,nstep

       ! Compute time-varying surface layer properties
       istat = sl_prelim(itplus(:,1), iqplus(:,1), iuplust(:,1), &
            ivplust(:,1), zp0, zztsl, spd_air=vmod, dir_air=vdir, &
            min_wind_speed=VAMIN)
       if (istat /= SL_OK) then
          call physeterror('pbl_ysu::pbl_ysu1', 'problem in sl_prelim() substep')
          return
       endif
       thair(:) = itplus(:,1) / exnerrpn(:,nkm1)
       istat = sl_sfclayer(thair, iqplus(:,1), vmod, vdir, zztsl, zztsl, &
            ztsrad, zqsurf, zz0, zz0t, zdlat, zfcor, &
            rib=ribsl, stabm=fmsl, stabt=fhsl)
       if (istat /= SL_OK) then
          call physeterror('pbl_ysu::pbl_ysu1', 'problem in sl_sfclayer() substep')
          return
       endif     

       ! Compute time-varying surface fluxes assuming fixed exchange coefficients
       rho(:) = iprest(:,1) / (RGASD * (itplus(:,1) * (1. + DELTA*iqplus(:,1))))
       wspd(:) = sqrt(iuplust(:,1)**2 + ivplust(:,1)**2) + 1.e-9
       do k=1,nkm1
          ki = nkm1 - (k-1)
          ustept(:,k) = iuplust(:,ki)
          vstept(:,k) = ivplust(:,ki)
       enddo
       call vint_thermo2mom(ustepm, ustept, zvcoef, ni, nkm1)
       call vint_thermo2mom(vstepm, vstept, zvcoef, ni, nkm1)
       lfrv(:) = max(sqrt(zbm(:) * (sqrt(ustepm(:,nkm1)**2 + vstepm(:,nkm1)**2))), 0.01)
       if (step == 1) lfrv0(:) = lfrv(:)
       lfh(:) = -rho(:)*CPD * (zalfat(:) + zbt(:)*thair(:)) * (lfrv(:)/lfrv0(:))
       lfsh(:) = -rho(:)*(zalfaq(:) + zbt(:)*iqplus(:,1)) * (lfrv(:)/lfrv0(:))       
       
       ! Pack moisture variables together
       imoist(:,1:nkm1) = iqplus(:,:)
       imoist(:,nkm1+1:2*nkm1) = iqcplus(:,:)
       imoist(:,2*nkm1+1:3*nkm1) = 0.

       ! Call YSU scheme to compute PBL tendencies
       call ysu2d(j=trnch, ux=iuplust, vx=ivplust, tx=itplus, qx=imoist, &
            p2d=iprest, p2di=ipresm, pi2d=iexner, &
            utnp=iudifv, vtnp=ivdifv, ttnp=itdifv, qtnp=imoisttnd, ndiff=NDIFF, &
            cp=CPD, g=GRAV, rovcp=CAPPA, rd=RGASD, rovg=RGASD/GRAV, &
            xlv=CHLC, rv=RGASV, ep1=DELTA, ep2=EPS1, karman=KARMAN, &
            dz8w2d=idzt, psfcpa=zp0, znt=zz0, ust=lfrv, &
            hpbl=hpbl, psim=fmsl, psih=fhsl, xland=lmg, &
            hfx=lfh, qfx=lfsh, wspd=wspd, br=ribsl, &
            dusfc=dusfc, dvsfc=dvsfc, dtsfc=dtsfc, dqsfc=dqsfc, &
            dt=stepdt, rcl=1., kpbl1d=kpbl, exch_hx=ikt, exch_mx=ikm, &
            wstar=zwstar, delta=ldelta, u10=zudiag, v10=zvdiag, &
            uox=zero, vox=zero, rthraten=ithrad, p2diORG=ipresmo, &
            ysu_topdown_pblmix=0, &
            ids=1, ide=ni, jds=1, jde=1, kds=1, kde=nkm1, &
            ims=1, ime=ni, jms=1, jme=1, kms=1, kme=nkm1, &
            its=1, ite=ni, jts=1, jte=1, kts=1, kte=nkm1)

       ! Unpack moisture variables
       iqdifv(:,:) = imoisttnd(:,1:nkm1)
       if (NDIFF > 1) iqcdifv(:,:) = imoisttnd(:,nkm1+1:2*nkm1)

       ! Apply substep tendencies to state variables
       iuplust(:,:) = iuplust(:,:) + iudifv(:,:) * stepdt
       ivplust(:,:) = ivplust(:,:) + ivdifv(:,:) * stepdt
       itplus(:,:) = itplus(:,:) + itdifv(:,:) * stepdt
       iqplus(:,:) = iqplus(:,:) + iqdifv(:,:) * stepdt
       if (NDIFF > 1) iqcplus(:,:) = iqcplus(:,:) + iqcdifv(:,:) * stepdt
       
    enddo TIME_STEPS

    ! Compute full-step tendencies as arrays are flipped
    idt = 1./dt
    do k=1,nkm1
       ki = nkm1 - (k-1)
       udifvt(:,k) = (iuplust(:,ki) - uplust(:,k)) * idt
       vdifvt(:,k) = (ivplust(:,ki) - vplust(:,k)) * idt
       ztdifv(:,k) = (itplus(:,ki) - ztplus(:,k)) * idt
       zqdifv(:,k) = (iqplus(:,ki) - zhuplus(:,k)) * idt
       zkm(:,k) = ikm(:,min(ki+1,nkm1))
       zkt(:,k) = ikt(:,min(ki+1,nkm1))
    enddo
    if (NDIFF > 1) then
       if (associated(zqcplus)) then
          do k=1,nkm1
             ki = nkm1 - (k-1)
             zqcdifv(:,k) = (iqcplus(:,ki) - zqcplus(:,k)) * idt
          enddo
       else
          zqcdifv(:,:) = 0.
       endif
    endif
    
    ! Interpolate momentum-level variables to thermodynamic levels
    call vint_thermo2mom(zudifv, udifvt, zvcoef, ni, nkm1)
    call vint_thermo2mom(zvdifv, vdifvt, zvcoef, ni, nkm1)

    ! Estimate near-surface diffusion coefficients
    zkm(:,nk) = min(KARMAN * lfrv(:) * zz0(:), zkm(:,nkm1))
    zkt(:,nk) = zkm(:,nk) / beta

    ! Successful completion
    status = RMN_OK

    return
  end function pbl_ysu1

end module pbl_ysu
