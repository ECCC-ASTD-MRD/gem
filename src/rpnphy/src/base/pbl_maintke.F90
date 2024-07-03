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

module pbl_maintke
  implicit none
  private
  public :: maintke

contains

  !/@*
  subroutine maintke(pvars, kount, ni, nk, nkm1, trnch)
    use, intrinsic :: iso_fortran_env, only: INT64
    use debug_mod, only: init2nan
    use tdpack_const, only: CPD, DELTA, GRAV, KARMAN
    use series_mod, only: series_xst
    use phy_options
    use phy_status, only: phy_error_L
    use phybusidx
    use phymem, only: phyvar
    use vintphy, only: vint_thermo2mom1
    use pbl_utils, only: WSTAR_MIN
    use pbl_ri, only: rpnint
    use pbl_ri_utils, only: PBL_RI_CK
    use pbl_mtke, only: moistke
    use pbl_mtke_utils, only: BLCONST_CU, BLCONST_CK
    implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
    !@Object interface for turbulent kinetic energy calculations
    !@Arguments
    !          - Input/Output -
    ! pvars    list of all phy vars (meta + slab data)
    !          - Input -
    ! kount    index of timestep
    ! trnch    number of the slice
    ! ni       horizontal dimension
    ! nk       vertical dimension
    ! nkm1     vertical scope of the operator

    type(phyvar), pointer, contiguous :: pvars(:)
    integer, intent(in) :: kount, trnch, ni, nk, nkm1

    !@Author J. Mailhot and B. Bilodeau (Jan 1999)
    !*@/
#include "phymkptr.hf"
    include "surface.cdk"
    include "tables.cdk"
    include "phyinput.inc"

    integer, external :: neark

    real, save :: tfilt = 0.1

    integer :: i, k, stat, ksl(ni)
    real :: cf1, cf2, eturbtau, uet, ilmot, &
         fhz, fim, fit, hst_local, &
         b(ni,nkm1*4), xb(ni), xh(ni), fbsurf(ni), zns

    real, pointer, dimension(:), contiguous :: zbt_ag, zfrv_ag, zftemp_ag, &
         zfvap_ag, zhst_ag, zilmo_ag, zz0_ag, ztsurf, zqsurf, zz0t_ag, &
         zdlat, zfcor
    real, pointer, dimension(:), contiguous   :: zalfat, zztsl, zh, zlh, &
         zwstar, zudiag, zvdiag, zhpar, zpmoins, zsdtswd, &
         zsdtsws, zwge, zwgmax, zwgmin, zdxdy, ztstar, zqstar

    real, pointer, dimension(:,:), contiguous :: tke, zenmoins, zkm, zkt, ztmoins, &
         ztve, zze, zzn, zwtng, zwqng, zuwng, zvwng, zqcmoins, zhumoins, &
         zsigm, zsigt, zsigw, zumoins, zvmoins, zfbl, zfblgauss, zfblnonloc, zfnn, &
         zftot, zlwc, ziwc, zqtbl, zturbreg, zzd, zgq, zgql, zgte, zgzmom, zrif, &
         zrig, zshear2, zvcoef, zc1pbl, zbuoy, zmrk2, zdiffen, zdissen, &
         zpri, zsige, zbuoyen, zshren, zgztherm, zfxp, zqwvar, zqwvarmoins, &
         zqwvarplus

    real, dimension(ni,nkm1) :: c, x, wk2d, enold, tmom, qmom, te, qce
    !----------------------------------------------------------------

    call init2nan(xb, xh, fbsurf)
    call init2nan(b, c, x, wk2d, enold, tmom, qmom, te, qce)

    MKPTR1D(zalfat, alfat, pvars)
    MKPTR1DK(zbt_ag, bt, indx_agrege, pvars)
    MKPTR1D(zdxdy, dxdy, pvars)
    MKPTR1D(zdlat, dlat, pvars)
    MKPTR1D(zfcor, fcor, pvars)
    MKPTR1DK(zfrv_ag, frv, indx_agrege, pvars)
    MKPTR1DK(zftemp_ag, ftemp, indx_agrege, pvars)
    MKPTR1DK(zfvap_ag, fvap, indx_agrege, pvars)
    MKPTR1D(zh, h, pvars)
    MKPTR1D(zhpar, hpar, pvars)
    MKPTR1DK(zhst_ag, hst, indx_agrege, pvars)
    MKPTR1DK(zilmo_ag, ilmo, indx_agrege, pvars)
    MKPTR1D(zlh, lhtg, pvars)

    MKPTR1D(zpmoins, pmoins, pvars)
    MKPTR1D(zqstar, qstar, pvars)
    MKPTR1DK(zqsurf, qsurf, indx_agrege, pvars)
    MKPTR1D(zsdtswd, sdtswd, pvars)
    MKPTR1D(zsdtsws, sdtsws, pvars)
    MKPTR1D(ztstar, tstar, pvars)
    MKPTR1DK(ztsurf, tsurf, indx_agrege, pvars)
    MKPTR1D(zudiag, udiag, pvars)
    MKPTR1D(zvdiag, vdiag, pvars)
    MKPTR1D(zwge, wge, pvars)
    MKPTR1D(zwgmax, wgmax, pvars)
    MKPTR1D(zwgmin, wgmin, pvars)
    MKPTR1D(zwstar, wstar, pvars)
    MKPTR1DK(zz0_ag, z0, indx_agrege, pvars)
    MKPTR1DK(zz0t_ag, z0t, indx_agrege, pvars)
    MKPTR1D(zztsl, ztsl, pvars)

    MKPTR2D(zbuoy, buoy, pvars)
    MKPTR2D(zbuoyen, buoyen, pvars)
    MKPTR2D(zc1pbl, c1pbl, pvars)
    MKPTR2D(zdiffen, diffen, pvars)
    MKPTR2D(zdissen, dissen, pvars)
    MKPTR2D(zfbl, fbl, pvars)
    MKPTR2D(zfblgauss, fblgauss, pvars)
    MKPTR2D(zfblnonloc, fblnonloc, pvars)
    MKPTR2D(zfnn, fnn, pvars)
    MKPTR2D(zftot, ftot, pvars)
    MKPTR2D(zfxp, fxp, pvars)
    MKPTR2D(zgq, gq, pvars)
    MKPTR2D(zgql, gql, pvars)
    MKPTR2D(zgte, gte, pvars)
    MKPTR2D(zgzmom, gzmom, pvars)
    MKPTR2D(zgztherm, gztherm, pvars)
    MKPTR2D(zhumoins, humoins, pvars)
    MKPTR2D(ziwc, iwc, pvars)
    MKPTR2D(zkm, km, pvars)
    MKPTR2D(zkt, kt, pvars)
    MKPTR2D(zlwc, lwc, pvars)
    MKPTR2D(zpri, pri, pvars)
    MKPTR2D(zqcmoins, qcmoins, pvars)
    MKPTR2D(zqtbl, qtbl, pvars)
    MKPTR2D(zrif, rif, pvars)
    MKPTR2D(zrig, rig, pvars)
    MKPTR2D(zshear2, shear2, pvars)
    MKPTR2D(zshren, shren, pvars)
    MKPTR2D(zsige, sige, pvars)
    MKPTR2D(zsigm, sigm, pvars)
    MKPTR2D(zsigt, sigt, pvars)
    MKPTR2D(zsigw, sigw, pvars)
    MKPTR2D(ztmoins, tmoins, pvars)
    MKPTR2D(zturbreg, turbreg, pvars)
    MKPTR2D(ztve, tve, pvars)
    MKPTR2D(zumoins, umoins, pvars)
    MKPTR2D(zuwng, uwng, pvars)
    MKPTR2D(zvcoef, vcoef, pvars)
    MKPTR2D(zvmoins, vmoins, pvars)
    MKPTR2D(zvwng, vwng, pvars)
    MKPTR2D(zwqng, wqng, pvars)
    MKPTR2D(zwtng, wtng, pvars)
    MKPTR2D(zzd, zd, pvars)
    MKPTR2D(zze, ze, pvars)

    MKPTR2D(zmrk2, mrk2, pvars)

    if (advectke) then
       MKPTR2D(zenmoins, enmoins, pvars)
       MKPTR2D(tke, enplus, pvars)
       MKPTR2D(zqwvarmoins, qwvarmoins, pvars)
       MKPTR2D(zqwvar, qwvarplus, pvars)
    else
       nullify(zenmoins, zqwvarmoins)
       MKPTR2D(tke, en, pvars)
       MKPTR2D(zqwvar, qwvar, pvars)
    endif
    if (znplus > 0) then
       MKPTR2D(zzn, znplus, pvars)
    else
       MKPTR2D(zzn, zn, pvars)
    endif

    eturbtau = delt

    !  filtre temporel

    cf1=1.0-tfilt
    cf2=tfilt

    !  initialiser e avec ea,z et h

    if (kount == 0) then
       
       if (.not.any([(any((/'tr/qtbl:m','qtbl     '/) == phyinread_list_s(i)),i=1,phyinread_n)]) .or. &
            fluvert == 'RPNINT') zqtbl = 0.

       INIT_TKE: if (any('en'==phyinread_list_s(1:phyinread_n))) then
          tke = max(tke,0.)
       else
          do k=1,nkm1
             !vdir nodep
             do i=1,ni
                tke(i,k)= max( etrmin2, blconst_cu*zfrv_ag(i)**2 * &
                     exp(-(zze(i,k)- zze(i,nkm1))/zh(i)) )
             end do
          end do
       endif INIT_TKE

       if (advectke) then
          do k=1,nkm1
             !vdir nodep
             do i=1,ni
                zenmoins(i,k) = tke(i,k)
             end do
          end do
       endif

    endif

    if (kount.gt.0) then
       call series_xst(zzn, 'LM', trnch)
       call series_xst(tke, 'EN', trnch)
    endif

    do k=1,nkm1
       !vdir nodep
       do i=1,ni
          if (advectke) then
             enold(i,k) = max(zenmoins(i,k),etrmin2)
             tke(i,k)   = max(tke(i,k),etrmin2)
          else
             enold(i,k) = tke(i,k)
          endif
       end do
    end do

    ! Compute surface layer scales
    stat = neark(zsigt, zpmoins, 1000., ni, nkm1, ksl)
    do i=1,ni
       xb(i)=1.0+DELTA*zhumoins(i,ksl(i))
       xh(i)=(GRAV/(xb(i)*ztve(i,ksl(i)))) * ( xb(i)*zftemp_ag(i) &
            + DELTA*ztve(i,ksl(i))*zfvap_ag(i) )
       fbsurf(i)=xh(i)
       xh(i)=max(0.0,xh(i))
       xh(i) = (zh(i)*xh(i))**(1./3.)
       zwstar(i) = xh(i)
       ztstar(i) = max(zftemp_ag(i)/max(zwstar(i),WSTAR_MIN),0.)
       zqstar(i) = max(zfvap_ag(i)/max(zwstar(i),WSTAR_MIN),0.)
    enddo

    ! Compute diffusion coefficients via TKE and mixing length
    if (fluvert == 'RPNINT') then

       call rpnint(tke, zkm, zkt, zpri, zrif, zrig, zbuoy, zshear2, &
            zbuoyen, zshren, zdiffen, zdissen, zqwvar, zzn, zzd, &
            enold, zumoins, zvmoins, ztmoins, zhumoins, zlwc, ziwc, zfxp, &
            zh, zlh, zpmoins, zz0_ag, zz0t_ag, zfrv_ag, zwstar, zqstar, &
            zhpar, ztsurf, zqsurf, zdlat, zfcor, &
            zsigm, zsigt, zgzmom, zgztherm, zdxdy, zmrk2, zvcoef, &
            eturbtau, kount, ni, nkm1)
       if (phy_error_L) return
       
    elseif (fluvert == 'MOISTKE') then

       ! Initialize non-gradient fluxes to zero
       zwtng = 0.
       zwqng = 0.
       zuwng = 0.
       zvwng = 0.
       
       call moistke(tke,enold,zzn,zzd,zrif,zrig,zbuoy,zshear2,zpri,zqtbl,zc1pbl,zfnn, &
            zfblgauss,zfblnonloc,zgte,zgq,zgql,zh,zlh,zhpar,zwtng,zwqng,zuwng,zvwng,&
            zumoins,zvmoins,ztmoins,ztve,zhumoins,zhumoins,zpmoins,zsigm,zsige,zsigw, &
            zze,zz0_ag,zgzmom,zfrv_ag,zwstar,fbsurf,zturbreg, &
            zmrk2,zvcoef,zdxdy,eturbtau,kount,trnch,ni,nkm1)
       if (phy_error_L) return
       
       do k=1,nkm1-1
          do i=1,ni
             zkm(i,k) = blconst_ck*zzn(i,k)*sqrt(tke(i,k))
             zkt(i,k) = zkm(i,k)*zpri(i,k)
             zbuoy(i,k) = -zkt(i,k)*zbuoy(i,k)
          enddo
       enddo

    endif

    ! Diagnose variables for turbulent wind (gusts and standard deviations)
    if (associated(zwge)) then
       call twind(zwge, zwgmax, zwgmin, zsdtsws, zsdtswd, &
            ztve, enold, zumoins, zvmoins, zudiag, zvdiag, &
            zsige, zze, zh, zfrv_ag, &
            zwstar, ni, nkm1)
       call series_xst(zwge   , 'WGE' , trnch)
       call series_xst(zwgmax , 'WGX' , trnch)
       call series_xst(zwgmin , 'WGN' , trnch)
       call series_xst(zsdtsws, 'SDWS', trnch)
       call series_xst(zsdtswd, 'SDWD', trnch)
    endif

    if (kount.eq.0) then
       call series_xst(zzn, 'LM', trnch)
       call series_xst(tke, 'EN', trnch)
    endif

    ! Set near-surface values (lowest thermo and diagnostic), using the Delage
    ! class of stability functions
    do i=1,ni
       if (fluvert /= 'MOISTKE' .and. fluvert /= 'RPNINT') tke(i,nkm1) = tke(i,nkm1-1)
       tke(i,nk) = tke(i,nkm1)
       uet = zfrv_ag(i)
       ilmot = zilmo_ag(i)
       if (ilmot > 0.) then
          !# hst_local is used to avoid division by zero
          hst_local = max( zhst_ag(i), zztsl(i)+1.)
          fhz = 1-zztsl(i)/hst_local
          fim = 0.5*(1+sqrt(1+4*d97_as*zztsl(i)*ilmot*beta/fhz))
          fit = beta*fim
       else
          fim=(1-dg92_ci*zztsl(i)*ilmot)**(-.16666666)
          fit=beta*fim**2
          fhz=1
       endif
       if (fluvert == 'RPNINT') then
          zns = KARMAN*(zztsl(i)+zz0_ag(i))/fim
          zkm(i,nkm1) = uet*zns*fhz
       else
          zzn(i,nkm1)   = KARMAN*(zztsl(i)+zz0_ag(i))/fim
          zkm(i,nkm1)   = uet*zzn(i,nkm1)*fhz
          zkt(i,nkm1)   = zkm(i,nkm1)*fim/fit
       endif
       zzn(i,nk) = KARMAN*zz0_ag(i)
       zkm(i,nk) = uet*zzn(i,nk)
       zkt(i,nk) = zkm(i,nk)/beta
    enddo

    return
  end subroutine maintke

end module pbl_maintke
