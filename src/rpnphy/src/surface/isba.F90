
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

!/@*
subroutine isba4(BUS, BUSSIZ, PTSURF, PTSURFSIZ, DT, KOUNT, N, M, NK)
   use tdpack_const, only: PI
   use phy_status, only: phy_error_L
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use sfc_options, only: atm_external, atm_tplus, radslope, vamin, sl_Lmin_soil, &
        zu, zt, impflx, thermal_stress, z0tevol
   use sfcbus_mod
   implicit none
!!!#include <arch_specific.hf>
   !@Object Multitasking of the surface scheme ISBA
   !@Arguments
   !               - Input/Output -
   ! BUS           bus of surface variables
   !               - Input -
   ! BUSSIZ        size of the surface bus
   ! PTSURF        surface pointers
   ! PTSURFSIZ     dimension of ptsurf
   ! KOUNT         number of timestep
   ! DT            timestep
   ! N             running length
   ! M             horizontal dimension
   ! NK            vertical dimension

   integer BUSSIZ, N, NK, KOUNT
   real DT
   real,target :: bus(bussiz)
   integer PTSURFSIZ
   integer PTSURF(PTSURFSIZ)

   !@Author S. Belair (January 1997)
   !*@/

   integer SURFLEN
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   integer, parameter :: INDX_SFC = INDX_SOIL
   logical, parameter :: ISBA_TDIAGLIM = .false.
   real, parameter :: HGHTM_DIAG0 = 10.
   real, parameter :: HGHTT_DIAG0 = 1.5

   integer I,m,zopt

   real,dimension(n) :: alphast, cd,      cg,     ch,     ct,     del
   real,dimension(n) :: dsnowdt, dwaterdt,freezs, hrsurf, leff,   psnz0
   real,dimension(n) :: rhoa,    rhomax,  rhost,  rsnow
   real,dimension(n) :: t2t,     tst,     tva,    vdir,   w2t 
   real,dimension(n) :: wft,     wgeq,    wgt,    wlt,    wrt,    wst
   real,dimension(n) :: z0tot,   zc1,     zc2,    zcs
   real,dimension(n) :: my_ta,   my_ua,   my_va,  vmod,   vmod0
   real,dimension(n) :: zref_sw_surf, zemit_lw_surf
   real,dimension(n) :: zu10, zusr
   real,dimension(n) :: zusurfzt, zvsurfzt, zqd
   real,dimension(n) :: zzenith

   real,pointer,dimension(:) :: cmu
   real,pointer,dimension(:) :: ctu
   real,pointer,dimension(:) :: hu
   real,pointer,dimension(:) :: ps
   real,pointer,dimension(:) :: tt
   real,pointer,dimension(:) :: uu
   real,pointer,dimension(:) :: vv
   real,pointer,dimension(:) :: z0h
   real,pointer,dimension(:) :: z0m
   real,pointer,dimension(:) :: zacoef
   real,pointer,dimension(:) :: zalfaq
   real,pointer,dimension(:) :: zalfat
   real,pointer,dimension(:) :: zalveg
   real,pointer,dimension(:) :: zalvis
   real,pointer,dimension(:) :: zbcoef
   real,pointer,dimension(:) :: zc1sat
   real,pointer,dimension(:) :: zc2ref
   real,pointer,dimension(:) :: zc3ref
   real,pointer,dimension(:) :: zcgsat
   real,pointer,dimension(:) :: zcveg
   real,pointer,dimension(:) :: zdlat
   real,pointer,dimension(:) :: zdrain
   real,pointer,dimension(:) :: zeflux
   real,pointer,dimension(:) :: zemisr
   real,pointer,dimension(:) :: zemsvc
   real,pointer,dimension(:) :: zfc
   real,pointer,dimension(:) :: zfcor
   real,pointer,dimension(:) :: zfdsi
   real,pointer,dimension(:) :: zfl
   real,pointer,dimension(:) :: zfq
   real,pointer,dimension(:) :: zfrv
   real,pointer,dimension(:) :: zfsolis
   real,pointer,dimension(:) :: zftemp
   real,pointer,dimension(:) :: zfv
   real,pointer,dimension(:) :: zfvap
   real,pointer,dimension(:) :: zgamveg
   real,pointer,dimension(:) :: zhst
   real,pointer,dimension(:) :: zhusurf
   real,pointer,dimension(:) :: zhv
   real,pointer,dimension(:) :: zilmo
   real,pointer,dimension(:) :: zisoil
   real,pointer,dimension(:) :: zlai
   real,pointer,dimension(:) :: zleg
   real,pointer,dimension(:) :: zler
   real,pointer,dimension(:) :: zles
   real,pointer,dimension(:) :: zletr
   real,pointer,dimension(:) :: zlev
   real,pointer,dimension(:) :: zmelts
   real,pointer,dimension(:) :: zmeltsr
   real,pointer,dimension(:) :: zoverfl
   real,pointer,dimension(:) :: zpcoef
   real,pointer,dimension(:) :: zpsn
   real,pointer,dimension(:) :: zpsng
   real,pointer,dimension(:) :: zpsnv
   real,pointer,dimension(:) :: zqdiag
   real,pointer,dimension(:) :: zqdiagtyp
   real,pointer,dimension(:) :: zqdiagtypv
   real,pointer,dimension(:) :: zqsurf
   real,pointer,dimension(:) :: zrainrate
   real,pointer,dimension(:) :: zresa
   real,pointer,dimension(:) :: zrgl
   real,pointer,dimension(:) :: zrnet_s
   real,pointer,dimension(:) :: zrootdp
   real,pointer,dimension(:) :: zrst
   real,pointer,dimension(:) :: zrunofftot
   real,pointer,dimension(:) :: zsnoal
   real,pointer,dimension(:) :: zsnoden
   real,pointer,dimension(:) :: zsnodp
   real,pointer,dimension(:) :: zsnoma
   real,pointer,dimension(:) :: zsnoro
   real,pointer,dimension(:) :: zsnowrate
   real,pointer,dimension(:) :: zstomr
   real,pointer,dimension(:) :: ztdiag
   real,pointer,dimension(:) :: ztdiagtyp
   real,pointer,dimension(:) :: ztdiagtypv
   real,pointer,dimension(:) :: zthetaa
   real,pointer,dimension(:) :: ztsoil1
   real,pointer,dimension(:) :: ztsoil2
   real,pointer,dimension(:) :: ztsrad
   real,pointer,dimension(:) :: ztsurf
   real,pointer,dimension(:) :: zudiag
   real,pointer,dimension(:) :: zudiagtyp
   real,pointer,dimension(:) :: zudiagtypv
   real,pointer,dimension(:) :: zvdiag
   real,pointer,dimension(:) :: zvdiagtyp
   real,pointer,dimension(:) :: zvdiagtypv
   real,pointer,dimension(:) :: zvegfrac
   real,pointer,dimension(:) :: zwfc
   real,pointer,dimension(:) :: zwflux
   real,pointer,dimension(:) :: zwsat
   real,pointer,dimension(:) :: zwsnow
   real,pointer,dimension(:) :: zwsoil1
   real,pointer,dimension(:) :: zwsoil2
   real,pointer,dimension(:) :: zwveg
   real,pointer,dimension(:) :: zwwilt
   real,pointer,dimension(:) :: zz0veg
   real,pointer,dimension(:) :: zz0tveg
   real,pointer,dimension(:) :: zza
   real,pointer,dimension(:) :: zzusl
   real,pointer,dimension(:) :: zztsl
   !      real,pointer,dimension(:,:) :: ts

   real,pointer,dimension(:) :: zfsd
   real,pointer,dimension(:) :: zfsf
   real,pointer,dimension(:) :: zcoszeni
   real,pointer,dimension(:) :: zutcisun
   real,pointer,dimension(:) :: zutcishade
   real,pointer,dimension(:) :: zradsun
   real,pointer,dimension(:) :: zradshade
   real,pointer,dimension(:) :: zwbgtsun
   real,pointer,dimension(:) :: zwbgtshade
   real,pointer,dimension(:) :: ztglbsun
   real,pointer,dimension(:) :: ztglbshade
   real,pointer,dimension(:) :: ztwetb
   real,pointer,dimension(:) :: zq1
   real,pointer,dimension(:) :: zq2
   real,pointer,dimension(:) :: zq3
   real,pointer,dimension(:) :: zq4
   real,pointer,dimension(:) :: zq5
   real,pointer,dimension(:) :: zq6
   real,pointer,dimension(:) :: zq7
   !------------------------------------------------------------------------

   !# In offline mode the t-step 0 is (correctly) not performed
   if (atm_external .and. kount == 0) return

   SURFLEN = M

   cmu      (1:n) => bus( x(bm,1,1)           : )
   ctu      (1:n) => bus( x(bt,1,indx_sfc)    : )
   z0h      (1:n) => bus( x(z0t,1,indx_sfc)   : )
   z0m      (1:n) => bus( x(z0,1,indx_sfc)    : )
   zacoef   (1:n) => bus( x(acoef,1,1)        : )
   zalfaq   (1:n) => bus( x(alfaq,1,1)        : )
   zalfat   (1:n) => bus( x(alfat,1,1)        : )
   zalveg   (1:n) => bus( x(alveg,1,1)        : )
   zalvis   (1:n) => bus( x(alvis,1,indx_sfc) : )
   zbcoef   (1:n) => bus( x(bcoef,1,1)        : )
   zc1sat   (1:n) => bus( x(c1sat,1,1)        : )
   zc2ref   (1:n) => bus( x(c2ref,1,1)        : )
   zc3ref   (1:n) => bus( x(c3ref,1,1)        : )
   zcgsat   (1:n) => bus( x(cgsat,1,1)        : )
   zcveg    (1:n) => bus( x(cveg,1,1)         : )
   zdlat    (1:n) => bus( x(dlat,1,1)         : )
   zdrain   (1:n) => bus( x(drain,1,1)        : )
   zeflux   (1:n) => bus( x(eflux,1,1)        : )
   zemisr   (1:n) => bus( x(emisr,1,1)        : )
   zemsvc  (1:n) => bus( x(emsvc,1,1)        : )
   zfc      (1:n) => bus( x(fc,1,indx_sfc)    : )
   zfcor    (1:n) => bus( x(fcor,1,1)         : )
   zfdsi    (1:n) => bus( x(fdsi,1,1)         : )
   zfl      (1:n) => bus( x(fl,1,1)           : )
   zfq      (1:n) => bus( x(fq,1,1)           : )
   zfrv     (1:n) => bus( x(frv,1,indx_sfc)   : )
   zftemp   (1:n) => bus( x(ftemp,1,indx_sfc) : )
   zfv      (1:n) => bus( x(fv,1,indx_sfc)    : )
   zfvap    (1:n) => bus( x(fvap,1,indx_sfc)  : )
   zgamveg  (1:n) => bus( x(gamveg,1,1)       : )
   zhst     (1:n) => bus( x(hst,1,indx_sfc)   : )
   zhusurf  (1:n) => bus( x(husurf,1,1)       : )
   zhv      (1:n) => bus( x(hv,1,1)           : )
   zilmo    (1:n) => bus( x(ilmo,1,indx_sfc)  : )
   zisoil   (1:n) => bus( x(isoil,1,1)        : )
   zlai     (1:n) => bus( x(lai,1,1)          : )
   zleg     (1:n) => bus( x(leg,1,1)          : )
   zler     (1:n) => bus( x(ler,1,1)          : )
   zles     (1:n) => bus( x(les,1,1)          : )
   zletr    (1:n) => bus( x(letr,1,1)         : )
   zlev     (1:n) => bus( x(lev,1,1)          : )
   zmelts   (1:n) => bus( x(melts,1,1)        : )
   zmeltsr  (1:n) => bus( x(meltsr,1,1)       : )
   zoverfl  (1:n) => bus( x(overfl,1,1)       : )
   zpcoef   (1:n) => bus( x(pcoef,1,1)        : )
   zpsn     (1:n) => bus( x(psn,1,1)          : )
   zpsng    (1:n) => bus( x(psng,1,1)         : )
   zpsnv    (1:n) => bus( x(psnv,1,1)         : )
   zqdiag   (1:n) => bus( x(qdiag,1,1)        : )
   zqdiagtyp(1:n) => bus( x(qdiagtyp,1,indx_sfc) : )
   zqdiagtypv(1:n) => bus( x(qdiagtypv,1,indx_sfc) : )
   zqsurf   (1:n) => bus( x(qsurf,1,indx_sfc) : )
   zrainrate(1:n) => bus( x(rainrate,1,1)     : )
   zresa    (1:n) => bus( x(resa,1,1)         : )
   zrgl     (1:n) => bus( x(rgl,1,1)          : )
   zrnet_s  (1:n) => bus( x(rnet_s,1,1)       : )
   zrootdp  (1:n) => bus( x(rootdp,1,1)       : )
   zrst     (1:n) => bus( x(rst,1,1)          : )
   zrunofftot(1:n) => bus( x(runofftot,1,indx_sfc) : )
   zsnoal   (1:n) => bus( x(snoal,1,1)        : )
   zsnoden  (1:n) => bus( x(snoden,1,1)       : )
   zsnodp   (1:n) => bus( x(snodp,1,indx_sfc) : )
   zsnoma   (1:n) => bus( x(snoma,1,1)        : )
   zsnoro   (1:n) => bus( x(snoro,1,1)        : )
   zsnowrate(1:n) => bus( x(snowrate,1,1)     : )
   zstomr   (1:n) => bus( x(stomr,1,1)        : )
   ztdiag   (1:n) => bus( x(tdiag,1,1)        : )
   ztdiagtyp(1:n) => bus( x(tdiagtyp,1,indx_sfc) : )
   ztdiagtypv(1:n) => bus( x(tdiagtypv,1,indx_sfc) : )
   ztsoil1  (1:n) => bus( x(tsoil,1,1)        : )
   ztsoil2  (1:n) => bus( x(tsoil,1,2)        : )
   ztsrad   (1:n) => bus( x(tsrad,1,1)        : )
   ztsurf   (1:n) => bus( x(tsurf,1,indx_sfc) : )
   zudiag   (1:n) => bus( x(udiag,1,1)        : )
   zudiagtyp(1:n) => bus( x(udiagtyp,1,indx_sfc) : )
   zudiagtypv(1:n) => bus( x(udiagtypv,1,indx_sfc) : )
   zvdiag   (1:n) => bus( x(vdiag,1,1)        : )
   zvdiagtyp(1:n) => bus( x(vdiagtyp,1,indx_sfc) : )
   zvdiagtypv(1:n) => bus( x(vdiagtypv,1,indx_sfc) : )
   zvegfrac (1:n) => bus( x(vegfrac,1,1)      : )
   zwfc     (1:n) => bus( x(wfc,1,1)          : )
   zwflux   (1:n) => bus( x(wflux,1,1)        : )
   zwsat    (1:n) => bus( x(wsat,1,1)         : )
   zwsnow   (1:n) => bus( x(wsnow,1,1)        : )
   zwsoil1  (1:n) => bus( x(wsoil,1,1)        : )
   zwsoil2  (1:n) => bus( x(wsoil,1,2)        : )
   zwveg    (1:n) => bus( x(wveg,1,1)         : )
   zwwilt   (1:n) => bus( x(wwilt,1,1)        : )
   zz0tveg  (1:n) => bus( x(z0tveg,1,1)       : )
   zz0veg   (1:n) => bus( x(z0veg,1,1)        : )
   zza      (1:n) => bus( x(za,1,1)           : )
   zzusl    (1:n) => bus( x(zusl,1,1)         : )
   zztsl    (1:n) => bus( x(ztsl,1,1)         : )

   zfsd     (1:n) => bus( x(fsd,1,1)          : )
   zfsf     (1:n) => bus( x(fsf,1,1)          : )
   zcoszeni (1:n) => bus( x(cang,1,1)         : )

   zutcisun     (1:n) => bus( x(yutcisun     , 1, indx_sfc) : )
   zutcishade   (1:n) => bus( x(yutcishade   , 1, indx_sfc) : )
   zwbgtsun     (1:n) => bus( x(ywbgtsun     , 1, indx_sfc) : )
   zwbgtshade   (1:n) => bus( x(ywbgtshade   , 1, indx_sfc) : )
   ztglbsun     (1:n) => bus( x(ytglbsun     , 1, indx_sfc) : )
   ztglbshade   (1:n) => bus( x(ytglbshade   , 1, indx_sfc) : )
   zradsun      (1:n) => bus( x(yradsun      , 1, indx_sfc) : )
   zradshade    (1:n) => bus( x(yradshade    , 1, indx_sfc) : )
   ztwetb       (1:n) => bus( x(ytwetb       , 1, indx_sfc) : )
   zq1          (1:n) => bus( x(yq1          , 1, indx_sfc) : )
   zq2          (1:n) => bus( x(yq2          , 1, indx_sfc) : )
   zq3          (1:n) => bus( x(yq3          , 1, indx_sfc) : )
   zq4          (1:n) => bus( x(yq4          , 1, indx_sfc) : )
   zq5          (1:n) => bus( x(yq5          , 1, indx_sfc) : )
   zq6          (1:n) => bus( x(yq6          , 1, indx_sfc) : )
   zq7          (1:n) => bus( x(yq7          , 1, indx_sfc) : )

   if (atm_tplus) then
      hu       (1:n) => bus( x(huplus,1,nk)      : )
      ps       (1:n) => bus( x(pplus,1,1)        : )
      tt       (1:n) => bus( x(tplus,1,nk)       : )
      zthetaa  (1:n) => bus( x(thetaap,1,1)      : )
      uu       (1:n) => bus( x(uplus,1,nk)       : )
      vv       (1:n) => bus( x(vplus,1,nk)       : )
   else
      hu       (1:n) => bus( x(humoins,1,nk)     : )
      ps       (1:n) => bus( x(pmoins,1,1)       : )
      zthetaa  (1:n) => bus( x(thetaa,1,1)       : )
      tt       (1:n) => bus( x(tmoins,1,nk)      : )
      uu       (1:n) => bus( x(umoins,1,nk)      : )
      vv       (1:n) => bus( x(vmoins,1,nk)      : )
   endif

   if (RADSLOPE) then
      zFSOLIS(1:n)   => bus( x(fluslop,1,1)      : )
   else
      zFSOLIS(1:n)   => bus( x(flusolis,1,1)     : )
   endif
 
   i = sl_prelim(tt,hu,uu,vv,ps,zzusl,VMOD0,VDIR,TVA,RHOA,min_wind_speed=VAMIN)
   if (sl_Lmin_soil > 0.) then
      VMOD(:) = VMOD0(:)
   elseif (i == SL_OK) then
      i = sl_prelim(tt,hu,uu,vv,ps,zzusl,VMOD,VDIR,TVA,RHOA, &
           min_wind_speed=2.5,min_wind_reduc='linear')
   endif
   if (i /= SL_OK) then
      call physeterror('isba', 'error returned by sl_prelim()')
      return
   endif

   call SOILI2(zTSOIL1, zWSOIL1, &
        zWSOIL2, zISOIL, &
        zSNOMA, zSNORO, &
        zVEGFRAC, zCGSAT, &
        zWSAT, zWWILT, &
        zBCOEF, zC1SAT, &
        zC2REF, zACOEF, &
        zPCOEF, zCVEG, &
        z0m, zz0veg, &
        CG, ZC1, ZC2, WGEQ, CT, ZCS, &
        zPSN, &
        zPSNG, zPSNV, &
        PSNZ0, Z0TOT, z0h, zz0tveg, N)

   call VEGI(zfsolis, &
        tt, ztsoil1, &
        hu, ps, &
        zWSOIL2, zRGL, &
        zLAI, zSTOMR, &
        zGAMVEG, zWWILT, &
        zWFC, zRST, &
        N)

   call drag7(ztsoil1, zWSOIL1, &
        zWVEG, zthetaa, &
        VMOD, VDIR, hu, &
        ps, zRST, &
        zVEGFRAC, &
        z0h, Z0TOT, zWFC, &
        zPSNG, zPSNV, &
        zLAI, zzusl, zztsl, &
        zdLAT, zFCOR, zRESA, &
        zILMO, zHST, &
        zFRV, zFTEMP, &
        zFVAP, &
        CH, CD, HRSURF, zHUSURF, &
        zHV, &
        DEL, zqsurf, &
        ctu, &
        N)
   if (phy_error_L) return

   call ebudget4(tt, &
        ztsoil1, ztsoil2, &
        zWSOIL2, zISOIL, &
        zWSNOW, zSNOMA, DT, &
        zSNOAL, zRAINRATE, &
        zfsolis, &
        zALVEG, zALVIS, zemisr, zemsvc,&
        zFDSI, zthetaa, &
        hu, ps, RHOA, &
        zVEGFRAC, HRSURF, &
        zHV, DEL, &
        zRESA, zRST, &
        CT, CG, ZCS, &
        zPSN, zPSNV, &
        zPSNG, zWSAT, &
        zROOTDP, zSNODP, &
        TST, T2T, &
        zRNET_S, &
        zfc, &
        zfv, &
        zLEG, zLEV, &
        zLES, zLER, &
        zLETR, zFL, &
        zEFLUX, &
        LEFF, DWATERDT, DSNOWDT, FREEZS, RHOMAX, &
        zMELTS, zMELTSR, &
        zFTEMP, zFVAP, &
        N)
 
   ! Local calculation of temperature and winds at chosen levels
   ! (for now, screen level heights are chosen, i.e. 1.5 and 10m)
   ! which are used in the calculation of "density of falling snow"
   ! in hydro1

   ! option for z0t calculation over vegetation
    if (z0tevol == 'ZILI95') then
      zopt = 9
   elseif (z0tevol == 'FIXED') then
      zopt = 0
   else
      call physeterror('isba', 'unknown option for z0tevol='//trim(z0tevol))
      return
   endif

   if (kount.eq.0) then
      my_ta = ztdiag
      my_ua = zudiag
      my_va = zvdiag
   else
      i = sl_sfclayer(zthetaa,hu,vmod,vdir,zzusl,zztsl,ztsoil1,zqsurf, &
           z0m,z0h,zdlat,zfcor,optz0=zopt,L_min=sl_Lmin_soil, &
           hghtm_diag=HGHTM_DIAG0,hghtt_diag=HGHTT_DIAG0,t_diag=my_ta,u_diag=my_ua, &
           v_diag=my_va,tdiaglim=ISBA_TDIAGLIM)
      if (i /= SL_OK) then
         call physeterror('isba', 'error returned by sl_sfclayer()')
         return
      endif
   endif

   call HYDRO3(DT, tt, ztsoil1, &
        zWSOIL1, &
        zWSOIL2, zISOIL, &
        zWSNOW, zWVEG, &
        zSNOMA, zSNOAL, &
        zSNORO, zSNODP, &
        zRAINRATE, zSNOWRATE, &
        my_ta, &
        my_ua, &
        my_va, &
        zLEV, zLETR, &
        zLEG, zLES, &
        ZC1, ZC2, zC3REF, &
        WGEQ, CT, zLAI, zVEGFRAC, &
        zROOTDP, zWSAT, &
        zWFC, zPSN, &
        zPSNG, zPSNV, &
        LEFF, DWATERDT, DSNOWDT, FREEZS, RHOMAX, &
        TST, &
        WGT, W2T, WFT, WLT, WRT, WST, ALPHAST, RHOST, &
        zDRAIN, zOVERFL, &
        RSNOW, N )

   call UPDATE4(ztsoil1, ztsoil2, &
        zwsoil1, zWSOIL2, &
        zISOIL, ZWSNOW, &
        zWVEG, ZSNOMA, &
        ZSNOAL, ZSNORO, &
        ZSNODEN, VMOD, CD, RHOA, &
        zfc, &
        zEFLUX, &
        TST, T2T, WGT, W2T, WFT, WLT, WRT, WST, &
        ALPHAST, RHOST, &
        cmu, zFQ, &
        zALFAT, &
        zALFAQ, &
        ctu, &
        zthetaa, &
        hu, &
        N)

   !# Compute values at the diagnostic level
   i = sl_sfclayer(zthetaa,hu,vmod,vdir,zzusl,zztsl,ztsoil1,zqsurf, &
        z0m,z0h,zdlat,zfcor,optz0=zopt,L_min=sl_Lmin_soil,spdlim=vmod, &
        hghtm_diag=zu,hghtt_diag=zt,t_diag=ztdiag,q_diag=zqdiag, &
        u_diag=zudiag,v_diag=zvdiag,tdiaglim=ISBA_TDIAGLIM)
   if (i /= SL_OK) then
      call physeterror('isba', 'error 2 returned by sl_sfclayer()')
      return
   endif
   zudiag = zudiag * vmod0 / vmod
   zvdiag = zvdiag * vmod0 / vmod

   !# Fill surface type-specific diagnostic values
   zqdiagtyp = zqdiag
   ztdiagtyp = ztdiag
   zudiagtyp = zudiag
   zvdiagtyp = zvdiag

   !# Compute diagnostics using vegetation-only roughness
   i = sl_sfclayer(zthetaa,hu,vmod,vdir,zzusl,zztsl,ztsoil1,zqsurf, &
        zz0veg,zz0tveg,zdlat,zfcor,optz0=zopt,   &
        L_min=sl_Lmin_soil,spdlim=vmod, &
        hghtm_diag=zu,hghtt_diag=zt,t_diag=ztdiagtypv,q_diag=zqdiagtypv, &
        u_diag=zudiagtypv,v_diag=zvdiagtypv,tdiaglim=ISBA_TDIAGLIM) 
   if (i /= SL_OK) then
      call physeterror('isba', 'error 3 returned by sl_sfclayer()')
      return
   endif
   zudiagtypv = zudiagtypv * vmod0 / vmod
   zvdiagtypv = zvdiagtypv * vmod0 / vmod

!VDIR NODEP
   do i=1,n

      ZTSURF   (I) = ZTSOIL1 (I)
      ZTSRAD   (I) = ZTSOIL1 (I)

      !# CALCULATE LAND-ATMOSPHERE OUTCOMING WATER FLUX
      zWFLUX(I) = RHOA(I)*zEFLUX(I)

      if (.not.IMPFLX) CTU (I) = 0.
      ! Fill runoff variable
      zRUNOFFTOT(i) = zOVERFL(i)

   end do
!-----------------------------------------------------
   !#TODO: at least 4 times identical code in surface... separeted s/r to call
   IF_THERMAL_STRESS: if (thermal_stress) then

      ! Compute wind at the globe sensor level
   i = sl_sfclayer(zthetaa,hu,vmod,vdir,zzusl,zztsl,ztsoil1,zqsurf, &
        zz0veg,zz0tveg,zdlat,zfcor,optz0=zopt,   &
        L_min=sl_Lmin_soil,spdlim=vmod, &
        hghtm_diag=zt,hghtt_diag=zt, &
        u_diag=zusurfzt,v_diag=zvsurfzt,tdiaglim=ISBA_TDIAGLIM) 

  if (i /= SL_OK) then
      call physeterror('isba', 'error 3 returned by sl_sfclayer()')
      return
   endif
   zusurfzt = zusurfzt * vmod0 / vmod
   zvsurfzt = zvsurfzt * vmod0 / vmod

      do i=1,N

         if (abs(zzusl(i)-zu) <= 2.0) then
            zu10(i) = sqrt(uu(i)**2+vv(i)**2)
         else
            zu10(i) = sqrt(zudiag(i)**2+zvdiag(i)**2)
         endif

         ! wind  at SensoR level zubos at z=zt
         if( (abs(zusurfzt(i)) >= 0.1) .and. (abs(zvsurfzt(i)) >= 0.1)) then
         zusr(i) = sqrt( zusurfzt(i)**2 + zvsurfzt(i)**2)
         else
         zusr(i) = zu10(i)
         endif

            zqd(i) = max( ZQDIAG(i) , 1.e-6) 

         zref_sw_surf (i)  = zalvis(i) * zfsolis(i)
         zemit_lw_surf(i)  = zfsolis(i) - zref_sw_surf (i) + zfdsi(i) - zrnet_s(i)

         zzenith(i) = acos(zcoszeni(i))      
         if (zfsolis(i) > 0.0) then
            zzenith(i) = min(zzenith(i), pi/2.)
         else
            zzenith(i) = max(zzenith(i), pi/2.)
         endif

      end do

      call SURF_THERMAL_STRESS(ZTDIAG, zqd,            &
           ZU10, zusr, ps,                             &
           ZFSD, ZFSF, ZFDSI, ZZENITH,                 &
           ZREF_SW_SURF,ZEMIT_LW_SURF,                 &
           Zutcisun ,Zutcishade,                       &
           zwbgtsun, zwbgtshade,                       &
           zradsun, zradshade,                         &
           ztglbsun, ztglbshade, ztwetb,               &
           ZQ1, ZQ2, ZQ3, ZQ4, ZQ5,                    &
           ZQ6, ZQ7, N)
   endif IF_THERMAL_STRESS

   !# FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
   call FILLAGG(BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_SOIL, &
        SURFLEN)

   return
end subroutine isba4
