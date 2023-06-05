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

module kfrpn
   implicit none
   private
   public :: kfrpn4

contains

  subroutine kfrpn4(ix, kx, flagconv, kkfc, ps, tp1, qp1, &
       ub, vb, scr3, &
       dtdt, dqdt, dudt, dvdt, &
       dudt1, dvdt1, dudt2, dvdt2, dudt3, dvdt3, sumdudt, sumdvdt, &
       dqcdt, dqidt, dqrdt, &
       sigt, dxdy, zcrr, gzm, &
       abeout, capeout, tauc, cinout, &
       wklclp, wklclout, areaup, clouds, dmfout, peffout, &
       umfout,  zbaseout,  ztopout, &
       wumaxout, rliqout, riceout, &
       rliq_int, rice_int, &
       rnflx, snoflx, &
       kount, xlat, mg, mlac, wstar, tstar, tke, kt, &
       coadvu, coadvv, coage, cowlcl, cozlcl, mrk2, critmask, delt)
    use, intrinsic :: iso_fortran_env, only: INT64
    use tdpack_const, only: CHLF, CPD, GRAV, PI, RGASD, TRPL
    use cnv_options
    use phy_options, only: dyninread_list_S, cmt_comp_diag
    use debug_mod, only: init2nan
    use tpdd, only: tpdd1
    use ens_perturb, only: ens_spp_get
    implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

    integer, intent(in) :: ix, kx
    real, intent(inout) :: flagconv(ix) ,kkfc(ix)
    real, intent(in)    :: ps(ix), tp1(ix,kx), qp1(ix,kx)
    real, intent(in)    :: ub(ix,kx), vb(ix,kx), scr3(ix,kx)
    real, intent(inout) :: dtdt(ix,kx), dqdt(ix,kx), dudt(ix,kx), dvdt(ix,kx)
    real, pointer, dimension(:,:) :: dudt1, dvdt1, dudt2, dvdt2, dudt3, dvdt3, sumdudt, sumdvdt
    real, intent(inout) :: dqcdt(ix,kx), dqidt(ix,kx), dqrdt(ix,kx)
    real, intent(in)    :: sigt(ix,kx), dxdy(ix), gzm(ix,kx)
    real, intent(out)   :: zcrr(ix), abeout(ix), capeout(ix), tauc(ix), cinout(ix)
    real, pointer       :: wklclp(:,:)
    real, intent(out)   :: wklclout(ix), areaup(ix,kx), clouds(ix,kx), dmfout(ix,kx)
    real, intent(out)   :: peffout(ix), umfout(ix,kx), zbaseout(ix), ztopout(ix)
    real, intent(out)   :: wumaxout(ix), rliqout(ix,kx), riceout(ix,kx)
    real, intent(out)   :: rliq_int(ix), rice_int(ix)
    real, intent(out)   :: rnflx(ix,kx), snoflx(ix,kx)
    integer, intent(in) ::  kount
    real, intent(in) :: xlat(ix), mg(ix), mlac(ix), wstar(ix), tstar(ix)
    real, intent(in) :: tke(ix,kx), kt(ix,kx)
    real, intent(inout) :: coadvu(ix), coadvv(ix), coage(ix), cowlcl(ix), cozlcl(ix)
    real, pointer :: mrk2(:,:)
    real, intent(in) :: critmask, delt

    !@Author Jack Kain and JM Fritsch (Oct 14,1990)

    !@Revision
    ! 001      Stephane Belair and Zonghui Huo (Dec.1994)
    !               inclusion into the RPN physics package
    ! 002      Stephane Belair (Nov. 1998) documentation of the code
    ! 003      Gerard Pellerin (Dec.1998) optimisation trigger
    ! 004      Richard Moffet  (Fev. 1999) RPN/CMC thermodynamic functions
    ! 005      Andre Methot correct bugs with NCA update and < 0
    ! 006      R. Moffet and S. Belair (Jan. 2000)
    !                   Cleaning of the code
    ! 007      S. Belair (Nov. 2000)  new output diagnostics
    ! 008      A.-M. Leduc (Jan 2001) Automatic arrays
    ! 009      S. Belair (Apr. 2001)  force detrainment in the upper portion
    !                                 of the updraft
    ! 010      S.Belair, G.Lemay, A-M Leduc (Apr 2001) Debugging
    ! 011      S.Belair (July 2001) convective transport of momentum
    ! 012      A-M. Leduc (Nov 2001) kfcp1 -->kfcp2 (changed arguments in call+
    !                                NCA becomes FLAGCONV)
    ! 013      A-M. Leduc (Feb 2002) reactivation of cloud fraction.
    ! 014      S. Belair, A-M. Leduc (nov 2002) add convective 1D counter kkfc.
    !                                argument kkfc added, kfcp2-->kfcp3.
    ! 015      K. Winger (nov 2002) conservation of water.
    ! 016      A-M. Leduc (Dec 2002) add switch ikfcpcp
    ! 017      L. Spacek (Sep 2003) do loop 107 changes LFS to LFS-1
    ! 018      A-M. Leduc (Sep 2003) Initialize LDB and minimize value of wabs.
    ! 019      L.Spacek (Dec 2003) loops 260,261 D/UMFOUT(kx-nk+1)=D/UMF(nk)
    ! 020      B. Bilodeau (Feb 2004) NUP is declared NK instead of IX
    ! 021      S. Belair (Fall 2003) Output precipitation fluxes
    ! 022      B. Bilodeau (Mar 2004) Add rliq_int and rice_int
    ! 023      B. Bilodeau (Mar 2004) Add "ramp" for wklcl
    ! 024      L.Spacek (May 2004)    In convective loop test
    !                                 the vertical loop to ITOP instead LMAX
    !                                 JKP replaced by LC+1
    ! 025      L.Spacek (Nov 2004)    cloud clean-up, multiplication of RLIQOUT
    !                                 RICEOUT by CLOUDS
    ! 026      B. Bilodeau (Dec 2004) Correct KKFC bug
    ! 027      B. Bilodeau (Jan 2005) Move the normalization of pcpn fluxes
    !                                 after ZCRR is recalculated
    ! 028      D. Figueras Nieto (Fev 2005) add D/UMFOUT(kx-nk+1)=D/UMF(nk)
    !                                       loop 170
    ! 029      D. Talbot (Spring 2005) wklcl = 0.02 at end of trigger function
    !                    to make the scheme more active
    ! 030      B.Dugas (Jun 2005)     Add surface modulation of WKLCL via KFCTRIGA
    ! 031      B. Bilodeau and S. Belair (Aug 2005) -
    !                       Correct RLIQOUT and RICEOUT bug.
    !                       Limit cloud fraction to values between (0.,1.)
    ! 032      R. McTaggart-Cowan and M. Desgagne (Jul 2006) -
    !                       Revert back to a more vectorizable form
    !                       for Vector processors
    ! 033      K. Winger (Sep 2004)   Correct normalization of precipitation fluxes

    ! 034      A-M. Leduc (Feb 2009) Add latitudinal variation to WKLCL over the ocean
    !                      using logical key kfctriglat. Add arguments MG ,  Mlac and
    !                      xlat in call. kfcp4 becomes kfcp5.
    ! 035      L. Spacek (Feb 2009)  Add UD=0, VD=0 in loop 117
    ! 036      L. Spacek    (Sep 2011)   Eliminate obsolete options
    ! 037      P. Vaillancourt,A.Zadra,M.Roch  (Dec 2012)
    !              Correction to avoid division by zero problem
    !              Internal report available(~armppva/RapportsInternes.html)
    ! 038      S. Gravel (Oct 2013)
    !              Introduce logical key kfcprod to calculate production term used
    !              by aqueous phase chemistry
    ! 039      P. Vaillancourt,V. Lee  (May 2014)
    !              Correction to avoid division by zero problem when TU10-TVQU(NK1)=0
    ! 040      P. Vaillancourt (Dec 2014)
    !              Output solid and liq precip fluxes independently; change units of
    !              updraft and downdraft mass fluxes for output
    !Object
    !          compute the effects of deep convection using
    !          the Kain-Fritsch convective paramerization scheme.(MKS units)

    !Arguments

    !          - Input -
    ! IX       X dimension of the model grid (NI)
    ! KX       Z dimension of the model grid (NK)

    !          - Input/Output -
    ! FLAGCONV   counter for whether convection is activated?

    !          - Input -
    ! PS       surface pressure
    ! TP1      temperature at time (T+1)
    ! QP1      specific humidity at time (T+1)
    ! UB       wind in X direction at time (T+1)
    ! VB       wind in Y direction at time (T+1)
    ! SCR3     vertical velocity at time (T+1)

    !          - Input/Output
    ! DTDT     convective effects of heating
    ! DQDT     convective effects of moistening
    ! DUDT     convective effects on u-momentum
    ! DVDT     convective effects of v-momentum
    ! DQCDT    cloud water due to the Kain and Fritsch scheme
    ! DQIDT    cloud ice due to the Kain and Fritsch scheme
    ! DQRDT    rain water/ snow due to the Kain and Fritsch scheme
    ! COADVU   cloud object advecting x-direction wind (m/s)
    ! COADVV   cloud object advecting y-direction wind (m/s)
    ! COAGE    cloud object age (s)
    ! COWLCL   cloud object updraft vertical motion at the LCL (m/s)
    ! COZLCL   cloud object cloud base height (m)
    ! WKLCLP   threshold vertical motion for triggering at time-plus

    !          - Input -
    ! GZM      geopotential
    ! SIGT     thermo sigma levels of the model
    ! DXDY     area of each tile of the grid
    ! XLAT     latitude(radians)
    ! MG       land-sea mask
    ! MLAC     fraction of lakes (mask)

    !          - Output -
    ! ZCRR     convective rainfall rate


    !          - Input -
    ! RAD      radius of the convective updraft at cloud base
    ! CDEPTH   minimum cloud depth

    !          - Output (diagnostics) -
    ! ABEOUT   available buoyant energy
    ! CAPEOUT  convective available potential energy
    ! TAUC     convection adjustment time
    ! CINOUT   convective inhibition
    ! AREAUP   cloud coverage area (m^2)
    ! CLOUDS   cloud fractional coverage area (fraction between 0 and 1)
    ! DMFOUT   downdraft mass flux
    ! PEFFOUT  precipitation efficiency
    ! UMFOUT   updraft mass flux
    ! ZBASEOUT cloud base height (i.e., LCL height)
    ! ZTOPOUT  cloud base top
    ! WUMAXOUT maximum velocity in the convective updrafts (m/s)
    ! RLIQOUT  mixing ratio of liquid water in the updrafts (kg/kg)
    ! RICEOUT  mixing ratio of ice in the updrafts (kg/kg)
    ! RLIQ_INT vertical integral of RLIQOUT
    ! RICE_INT vertical integral of RICEOUT


    !  References:
    !  Fritsch and Chappell (1980), J. Atmos. Sci., 1722-1733.
    !  Zhang and Fritsch (1986), J. Atmos. Sci., 1913-1943.
    !  Kain and Fritsch (1990), J. Atmos. Sci., 2784-2802.

    !Notes:
    !
    ! Suggested values for lamda and GKI
    ! 1. cmt as implemented in phase 2 (with bug)
    !      cmt_ecmwf_lambda = 2.  
    !      cmt_gki_cup = 0.0
    !      cmt_gki_cdown = 0.0
    ! 2. like opsph2 but debugged and applied over whole cloud
    !      cmt_ecmwf_lambda = 2.  
    !      cmt_gki_cup = 0.0
    !      cmt_gki_cdown = 0.0
    ! 3. zero-drag from Romps - results in max tendencies
    !      cmt_ecmwf_lambda = 0.   
    !      cmt_gki_cup = 0.0
    !      cmt_gki_cdown = 0.0
    ! 4. standard GKI; not including cmt in downdraft 
    !      cmt_ecmwf_lambda = 0.   
    !      cmt_gki_cup = 0.7
    !      cmt_gki_cdown = 0.0
    ! 5. GKI; stronger pressure effect - smaller tendencies
    !      cmt_ecmwf_lambda = 0.   
    !      cmt_gki_cup = 0.9
    !      cmt_gki_cdown = 0.0
    ! 6. GKI; lower pressure effect - stronger tendencies
    !      cmt_ecmwf_lambda = 0.   
    !      cmt_gki_cup = 0.3
    !      cmt_gki_cdown = 0.0
    ! 7. Tendencies=0. with following parameters (analytically)
    !      cmt_ecmwf_lambda = 0.   
    !      cmt_gki_cup = 1.0
    !      cmt_gki_cdown = 1.0
    ! 8. for testing, e.g. GKI-orig with pressure term in downdraft
    !      cmt_ecmwf_lambda = 0.   
    !      cmt_gki_cup = 0.7
    !      cmt_gki_cdown = 0.7


    !THINGS THAT REMAINS TO BE DONE IN THIS SCHEME:
    !=============================================

    ! 1) Verify the sensitivity to RAD (=3000 for the moment)

    ! 2) Convert to "mixing ratios" at the beginning of the
    !    subroutine

    ! 3) Allow detrainment of ice nucleis, as well as cloud water

    !*
    real, parameter :: WU_MIN=0.

    integer i,k,ifexfb,iflag
    integer klcl,klclm1,kmin,kmix,kfrz,kpbl,lc,lcl,ldt
    integer ldb,let,lmax,low,ltop,ltop1,ltopm1,lvf,lfs,ml
    integer nd,nd1,ndk,ncount,nic,nj,nk,nk1,nm
    integer nstep,ntc,nlayrs,nupnk
    integer ktop

    real a1,abe,abeg,aice,ainc,aincmx,aliq,au0
    real be,bice,bliq,boterm,cldhgt
    real c5,cbh,cice,cliq,clvf,cndtnf,cpm,cporq,cpr,cv
    real dabe,dq,dlp,dpt,dtt,dtt1,dzz,dliq,dice
    real devdmf,dmfmin,dmflfs,dptt
    real zlcl
    real betaw,betai,gamw,gami
    real dpdd,dpddmx,dpptdf,ddinc,dqsdt
    real dtmp,dtime,dtmltd,dtlcl
    real dumfdp
    real es,ee,ee1,ee2,emix,effq,enterm
    real f1,f2,fabe,frc,frc1
    real gdt,gdry,liqfrac,icefrac
    real p00,p165,pef,peff,pefmin,pefmax,pefcbh
    real pmid,ppr,pptflx,pptfl2,pptmlt
    real qenv,qese,qnewic,qnewlq,qnwfrz,qs,qsrh
    real r1,rl,rf,rei,rdd
    real rovg,rocpq,rate,rcbh,rced,rhbc,rhic,rtmp
    real stab,shsign,sumflx
    real tvlcl,tlog,tenv,dtenv,tven,tdpt,tbfrz,ktlcl
    real tvavg,tvbar,thata,thtfc,thtudl,thtmin,thttmp
    real ttmp,ttemp,tmpliq,tmpice,ttfrz,trppt,tsat
    real tder,tder2,t1rh
    real oneovg, maxzpar, dpmix, betaent
    real tu10,tu95,tudl
    real ud1,ud2,usr,udlbe,updinc,updin2
    real upold,upnew
    real vws,vmflcl
    real wlcl,wkl,wsigne,wtw
    real xls0,xls1,xlv0,xlv1
    real dpup , dpdown
    real dqdtdk,dqcdtdk
    real wklcld, trigger_energy
    real kfscale,cin,cape,wumean
    real denom,eps,tvdiff
    real intdudt,intdvdt,intdp,intu,intv
    real t1,t2,cmean,dcape,tp,tm,wlcl0
    real onemincu,onemincd,norm

    ! For the optimisation and thermodynamic functions

    integer :: l5,klm,kl,llfc,mxlayr

    real :: dxsq,rholcl,wabs,zmix,p300,thta

    real :: wsmax,wsmin

    logical, dimension(ix) :: activ

    integer, dimension(kx) :: nup

    real, dimension(ix) :: dpthmxg, pmixg,  tmixg, qmixg, &
         thmixg, tlclg, plclg, &
         wklcla, psb, thvmixg, tkemixg, &
         trigs,   trige, radmult, dpddmult, w_wsmin, w_wsmax, &
         acrate
    real, dimension(kx) :: ddr,     ddr2,   der,    der2,  detic, detic2, &
         detlq,   detlq2, dmf,    dmf2,  domgdp, &
         dtfm,    ems,    emsd,   eqfrc, exn, omgaeff,           &
         omga,    pptice, pptliq, qadv,  qd,    qdt,    &
         qg,      qicout, qlqout, qpa,   qu,     &
         ratio2,  rice,   rliq,   tg,    theted,thetee, &
         theteu,  thadv,  thpa,  thtad, thtag,  &
         thtau,   thta0,  thtes,  thtesg,tu,    tvd,    &
         tvg,     tvqu,   tvu,    tz,    udr,   udr2,   &
         uer,     uer2,   umf,    umf2,  wspd,   &
         wu,      uu,     vu,     ud,    vd,    ug,     &
         vg,      upa,    vpa,    uadv,  vadv,  qlpa,   &
         qipa,    qlg,    qig,    qtpa,  qtdt,  qtg,    &
         thpai,   qtpai,  qlpai,  qipai, ta,    tb,     &
         tc,      td_thta,td_qt,  td_ql, td_qi, tri_thta,&
         tri_qt,  tri_ql, tri_qi, workk, emf,   td_q,   &
         tri_q,   qpai,   td_u,   td_v,  tri_u, tri_v,  &
         upai,    vpai, cmtdenom, ct1, ct2, ct3, ct4, delpup
    real, dimension(kx+1)  :: omg
    real, dimension(ix,kx) :: tt0,  tv00,   q00,  u00,    v00,   ww0,    &
         dzp,  dpp,    qst1, pp0,    z0g,     &
         sigkfc, ql0,  qi0,  thv0, tke0,    &
         kt0,  wklcl0, idudt1, idvdt1, idudt2, idvdt2, idudt3, idvdt3
    character(len=64), dimension(ix) :: deeptrig

    external tpmix
    external condload_safe
    external envirtht
    external prof5

    ! Basic parameters
#include "clefcon.cdk"
    include "phyinput.inc"

    ! User-adjustable parameters
    real, parameter :: DETREG=0.5                             !Level of dynamic detrainment for cloud ensemble
    real, parameter :: DETTOT=0.                              !Total dynamic detrainment for cloud ensemble
    logical, parameter :: DOWNDRAFT=.true.                    !Include downdraft calculations
    logical, parameter :: TOTAL_WATER=.false.                 !Use total water and condensate in mass flux equations
    logical, parameter :: MIXING_RATIO=.false.                !Convert to/from mixing ratio for the scheme

    ! Local variables
    real :: rad, cdepth, timec, timer
    real :: nuer,nudr
    character(len=32), dimension(ix) :: scheme

    ! character(len=16) :: trigger_type='original' !original,kain04,scaled,tke

    call init2nan(dpthmxg, pmixg,  tmixg, qmixg)
    call init2nan(thmixg, tlclg, plclg)
    call init2nan(wklcla, psb, thvmixg,tkemixg)
    call init2nan(ddr,     ddr2,   der,    der2,  detic, detic2)
    call init2nan(detlq,   detlq2, dmf,    dmf2,  domgdp)
    call init2nan(dtfm,    ems,    emsd,   eqfrc, exn)
    call init2nan(omga,    pptice, pptliq, qadv,  qd,    qdt)
    call init2nan(qg,      qicout, qlqout, qpa,   qu)
    call init2nan(ratio2,  rice,   rliq,   tg,    theted,thetee)
    call init2nan(theteu,  thadv,  thpa,   thtad, thtag)
    call init2nan(thtau,   thta0,  thtes,  thtesg,tu,    tvd)
    call init2nan(tvg,     tvqu,   tvu,    tz,    udr,   udr2)
    call init2nan(uer,     uer2,   umf,    umf2,  wspd)
    call init2nan(wu,      uu,     vu,     ud,    vd,    ug)
    call init2nan(vg,      upa,    vpa,    uadv,  vadv,  qlpa)
    call init2nan(qipa,    qlg,    qig,    qtpa,  qtdt,  qtg)
    call init2nan(thpai,   qtpai,  qlpai,  qipai, ta,    tb)
    call init2nan(tc,      td_thta,td_qt,  td_ql, td_qi, tri_thta)
    call init2nan(tri_qt,  tri_ql, tri_qi, workk, emf,   td_q)
    call init2nan(tri_q,   qpai,   td_u,   td_v,  tri_u, tri_v)
    call init2nan(upai,    vpai, omg, radmult, dpddmult)
    call init2nan(omgaeff, cmtdenom)
    call init2nan(tt0,  tv00,   q00,  u00,    v00,   ww0)
    call init2nan(dzp,  dpp,    qst1, pp0,    z0g)
    call init2nan(sigkfc, ql0,  qi0,  thv0, tke0)
    call init2nan(kt0,  wklcl0)
    call init2nan(ct1, ct2, ct3, ct4, delpup)
    call init2nan(idudt1, idvdt1, idudt2, idvdt2, idudt3, idvdt3)

    cdepth = kfcdepth
    scheme(:) = ens_spp_get('deep', mrk2, default=deep)

    ! "RAMP" FOR WKLCL :
    ! ================

    !   WKLCL WILL INCREASE FROM KFCTRIG4(3) TO KFCTRIG4(4)
    !   BETWEEN TIMESTEPS KFCTRIG4(1) AND KFCTRIG4(2)

    trigs(:) = kfctrig4(3)
    trige(:) = ens_spp_get('kfctrig4', mrk2, default=kfctrig4(4)) 
    if (kount <= int(kfctrig4(1))) then
       wklcla(:) = trigs(:)
    else if (kount > int(kfctrig4(2))) then
       wklcla(:) = trige(:)
    else
       wklcla(:) = trigs(:) + (real(kount) - kfctrig4(1)) / &
            (kfctrig4(2) - kfctrig4(1)) * (trige(:) - trigs(:))
    endif


    ! Latitudinal ramp for WKLCL :
    ! ============================

    ! WKLCL will take on different values:
    ! over land and lakes: we kee the value set by the "ramp" above
    ! over sea water:
    !     for |lat| >= TRIGLAT(2) we keep value set by the "ramp" above
    !     for |lat| <= TRIGLAT(1) we use the new value KFCTRIGL
    !     and linear interpolation in between TRIGLAT(1) and TRIGLAT(2)

    if (KFCTRIGLAT) then
       do I=1,IX
          if (abs(XLAT(I)) .le. TRIGLAT(1).and. &
               MG(I) .le. critmask .and. &
               MLAC(I) .le. critmask) then
             WKLCLA(I)= KFCTRIGL
          else if (abs(XLAT(I)).gt.TRIGLAT(1) .and. &
               abs(XLAT(I)).lt.TRIGLAT(2) .and. &
               MG(I) .le. critmask .and. &
               MLAC(I) .le. critmask) then
             WKLCLA(I)= ( ((abs(XLAT(I))-TRIGLAT(1))/ &
                  (TRIGLAT(2)-TRIGLAT(1)))* &
                  (WKLCLA(I)-KFCTRIGL) ) + KFCTRIGL
          endif
       end do
    endif

    ! Convective velocity scale ramp for WKLCL :
    ! ============================================

    ! REPLACES the latitudinal ramp above.

    ! WKLCL will take on different values:
    ! over land and lakes: we keep the value set by the "ramp" above
    ! over sea water:
    !     for wstar <= KFCTRIGW(1) we use KFCTRIGW(3)
    !     for wstar >= KFCTRIGW(2) we use KFCTRIGW(4)
    !     and linear interpolation in between KFCTRIGW(1) and KFCTRIGW(2)

    if (any(abs(kfctrigw) > 0.)) then
       wsmin = kfctrigw(1)
       wsmax = kfctrigw(2)
       w_wsmin(:) = ens_spp_get('kfctrigwl', mrk2, default=kfctrigw(3))
       w_wsmax(:) = ens_spp_get('kfctrigwh', mrk2, default=kfctrigw(4))
       where (mg <= CRITMASK .and. mlac <= CRITMASK)
          WKLCLA = min(max(w_wsmin + (wstar - wsmin) / (wsmax - wsmin) * (w_wsmax - w_wsmin) , &
               min(w_wsmin,w_wsmax)), max(w_wsmin,w_wsmax))
       endwhere
    endif

    ! ===================================================

    if (KFCTRIGA.gt.0.0) then
       !  (B. Dugas, June 15 2005)
       !  IN WHAT FOLLOWS, WE ASSUME THAT
       !     KFCTRIG * RESOLUTION = CONSTANT,
       !   SO THAT
       !     KFCTRIG(RES2) = KFCTRIG(RES1) * RES1 / RES2
       !  THE 0.17 AND 0.01 LIMITS ARE THOSE THAT ARE DEEMED
       !  TO BE APPROPRIATE FOR 10KM AND 170KM, RESPECTIVELY.
       !  THESE LAST TWO VALUES DEFINE THE RESOLUTION INTERVAL
       !  AT WHICH THIS KF CONVECTION CODE WILL BE USED IN THE
       !  SGMIP AND/OR OPERATIONAL FRAMEWORK
       !  CONVERT KFCTRIGA TO METRES
       KFSCALE = 1000. * KFCTRIGA
       WKLCLA(:) = min( 0.17, &
            max( 0.01, &
            WKLCLA(:) * KFSCALE &
            / sqrt( DXDY(:) ) &
            ) &
            )
    endif

    ! Relax advected trigger velocity towards local value
    if (associated(wklclp) .and. kfctrigtau > 0.) then
       if (kount == 0) then
          if (.not.(any(dyninread_list_s == 'wklcl') .or. &
               any(phyinread_list_s(1:phyinread_n) == 'tr/wklcl:p'))) then
             do k=1,kx
                wklclp(:,k) = wklcla(:)
             enddo
          endif
       endif
       do k=1,kx
          wklclp(:,k) = (wklclp(:,k) + wklcla(:)*delt/kfctrigtau) / &
               (1. + delt/kfctrigtau)
       enddo
    endif

    ! 2.  USEFUL PARAMETERS
    ! =====================

    EPS = 1.E-8
    P00     = 1.E5
    CV      = 717.
    RHIC    = 1.
    RHBC    = 0.90
    TTFRZ =   268.16
    TBFRZ   = 248.16
    C5      = 1.0723E-3
    RATE    = 0.01
    MAXZPAR = 3.5E3
    DPMIX   = 60.E2
    BETAENT = 1.05
    PEFMIN  = 0.2
    PEFMAX = 0.9

    ROVG   = RGASD/GRAV
    GDRY   = -GRAV/CPD
    ONEOVG = 101.9368
    KL     = KX
    KLM    = KL-1

    BETAW = 6822.459384
    GAMW  = 5.13948
    BETAI = 6295.421
    GAMI  = 0.56313

    ONEMINCU = 1. - cmt_gki_cup
    ONEMINCD = 1. - cmt_gki_cdown

    ! Define constants for calculation of latent heating

    XLV0 = 3.147E+6
    XLV1 = 2369.
    XLS0 = 2.905E+6
    XLS1 = 259.532

    ! Define constants for calculation of saturation vapor
    ! pressure according to Buck (JAM, December 1981)

    ALIQ = 613.3
    BLIQ = 17.502
    CLIQ = 4780.8
    DLIQ = 32.19
    AICE = 613.2
    BICE = 22.452
    CICE = 6133.0
    DICE = 0.61

    ! Flag for feedback to explicit moisture
    ! IFEXFB = O  No feedback of the convective scheme on the grid-scale
    !             cloud and rainwater variables
    ! IFEXFB = 1  production term calculated for potential feedback calculation
    !PV note: feedback is not and should not be allowed for grid scale rainwater
    !         unless kfcp code is modified

    IFEXFB = 0
    if (kfcprod) IFEXFB = 1

    ! Retrieve stochastic parameter information on request
    radmult(:) = ens_spp_get('crad_mult', mrk2, default=1.)
    dpddmult(:) = ens_spp_get('dpdd_mult', mrk2, default=1.)
    deeptrig(:) = ens_spp_get('deeptrig', mrk2, default='ORIGINAL')
    acrate(:) = ens_spp_get('deeprate', mrk2, default=RATE)

    ! =============================================================


    ! 4. INPUT A VERTICAL SOUNDING (ENVIRONMENTAL)
    ! ============================================
    do I = 1, IX
       PSB(I)=PS(I)*1.E-3
    enddo
    do K=1,KX
       do I=1,IX
          SIGKFC(I,K) = SIGT(I,K)
       end do
    end do

    do K = 1, KX
       NK = KX-K+1
       do I = 1, IX

          ! The sounding
          ! In the RPN physics, the K indices goes
          ! from 1 at the top of the model to KX at the surface.
          ! In the Kain-Fritsch subroutine, we use the opposite (of course).
          ! The switch is done in this loop.

          PP0(I,K) = 1.E3*(SIGKFC(I,NK)*PSB(I))
          TT0(I,K) = TP1(I,NK)
          U00(I,K) = UB(I,NK)
          V00(I,K) = VB(I,NK)
          if (MIXING_RATIO) then
             Q00(I,K) = max(qp1(i,nk)/(1.-qp1(i,nk)),1.0E-10 ) !convert to mixing ratio
             QL0(I,K) = 0.
             QI0(I,K) = 0.
          else
             Q00(I,K) = max(qp1(i,nk),1.0E-10 )
             QL0(I,K) = 0.
             QI0(I,K) = 0.
          endif
          KT0(I,K) = KT(I,NK)
          WKLCL0(I,K) = 0.
          if (associated(wklclp)) WKLCL0(I,K) = WKLCLP(I,NK)
          ! add the definition of z0g (geopotential height) from gzm(geopotential)

          Z0G(I,K) = GZM(I,NK)/grav

          ! If Q0 is above saturation value, reduce it to saturation level

          ES=ALIQ*exp((BLIQ*TT0(I,K)-CLIQ)/(TT0(I,K)-DLIQ))
          ES=min( ES, 0.5*PP0(I,K) )
          QST1(I,K)=0.622*ES/(PP0(I,K)-ES)
          QST1(I,K) = max( min(QST1(I,K),0.050) , 1.E-6 )
          Q00(I,K)  = min(QST1(I,K),Q00(I,K))
          TV00(I,K) = TT0(I,K) * (1. + 0.608*Q00(I,K) )
          ROCPQ=0.2854*(1.-0.28*Q00(I,K))
          THV0(I,K) = TV00(I,K)*(1.E5/PP0(I,K))**ROCPQ
          WW0(I,K)  = SCR3(I,NK)
          TKE0(I,K) = TKE(I,NK)

       end do
    end do

    do I=1,IX
       DPP(I,1)   = ( 1. - 0.5*(SIGKFC(I,KX)+SIGKFC(I,KX-1))) * PSB(I)*1.E3    
       DPP(I,KX)  = 0.5 * ( SIGKFC(I,2) -SIGKFC(I,1)    ) * PSB(I)*1.E3        
    end do

    do K=2,KX-1
       NK = KX-K+1
       do I=1,IX
          DPUP      = ( SIGKFC(I,NK)-SIGKFC(I,NK-1) ) * PSB(I)*1.E3
          DPDOWN    = ( SIGKFC(I,NK+1)-SIGKFC(I,NK) ) * PSB(I)*1.E3
          DPP(I,K)  = 0.5 * (DPUP+DPDOWN)
       end do
    end do

    ! pass gzm as argument, so don't need to calculate z0g here,
    ! only need to calculate the dzp (thickness) from the inverted z0g.

    do K = 2, KX
       do I = 1, IX
          DZP(I,K-1) = Z0G(I,K)-Z0G(I,K-1)
       end do
    end do

    ! at the top: do not put the thickness to zero but instead
    ! put the level kx-1 in kx

    do I=1,IX
       DZP(I,KX) = DZP(I,KX-1)
    end do

    ! MAIN LOOP ACROSS I FOR KFC
    ! ================

    MAIN_I_LOOP: do I = 1, IX

       ! 5. COUNTER "FLAGCONV"
       ! ================


       ! If FLAGCONV > 0
       ! ==> convection is already activated
       !     the tendencies on T, Q, QC, and
       !     QR do not change.
       !     goto the next column

       ! If FLAGCONV =< 0
       ! ==> convective tendencies and diagnostics
       !     are put to "0"


       ACTIV(I) = .false.

       if (nint(FLAGCONV(I)) == 1) KKFC(I)=1.

       if (FLAGCONV(I).le.0) then

          RLIQ_INT(I) = 0.
          RICE_INT(I) = 0.
          if (cmt_comp_diag) then
             IDUDT1(I,:)  = 0.0 ;IDUDT2(I,:)  = 0.0 ;IDUDT3(I,:) = 0.0 
             IDVDT1(I,:)  = 0.0 ;IDVDT2(I,:)  = 0.0 ;IDVDT3(I,:)  = 0.0
             DUDT1(I,:)   = 0.0 ;DUDT2(I,:)  = 0.0 ;DUDT3(I,:)  = 0.0 
             DVDT1(I,:)   = 0.0 ;DVDT2(I,:)  = 0.0 ;DVDT3(I,:)  = 0.0 
             SUMDUDT(I,:) = 0.0 ;SUMDVDT(I,:)  = 0.0
          endif
          CMTDENOM(:)  = 0.0 

          do K=1,KX
             DTDT(I,K)  = 0.0
             DQDT(I,K)  = 0.0
             DQCDT(I,K) = 0.0
             DQIDT(I,K) = 0.0
             DUDT(I,K)  = 0.0
             DVDT(I,K)  = 0.0
             DQRDT(I,K) = 0.0

             UMFOUT(I,K) = 0.
             DMFOUT(I,K) = 0.

             RLIQOUT(I,K) = 0.
             RICEOUT(I,K) = 0.

             RNFLX(I,K)   = 0.
             SNOFLX(I,K)  = 0.

             AREAUP(I,K) = 0.
             CLOUDS(I,K) = 0.

          end do

          ZCRR(I) = 0.0

          ABEOUT(I)   = 0.
          CAPEOUT(I)  = 0.
          CINOUT(I)   = 0.
          TAUC(I)     = 0.
          PEFFOUT(I)  = 0.
          ZBASEOUT(I) = 0.
          ZTOPOUT(I)  = 0.
          WUMAXOUT(I) = 0.
          WKLCLOUT(I) = 0.
          COADVU(I)   = 0.
          COADVV(I)   = 0.

       else if (FLAGCONV(I).gt.0) then
          ACTIV(I) = .true.
          FLAGCONV(I) = FLAGCONV(I) - 1
          ZCRR(I) = 1000.*ZCRR(I)
          if (MIXING_RATIO) then
             do k=1,kx
                nk = kx-k+1
                dqdt(i,k) = dqdt(i,k) * (1.+Q00(i,nk))**2
                dqcdt(i,k) = dqcdt(i,k) * (1.+Q00(i,nk))**2
                dqidt(i,k) = dqidt(i,k) * (1.+Q00(i,nk))**2
             enddo
          endif
          goto 325 !Active point tendencies persist: no calculations done

       endif

       ! Establish local cloud properties
       rad = kfcrad
       if (mg(i) <= CRITMASK .or. mlac(i) <= CRITMASK) rad = kfcradw
       rad = rad * radmult(i)
       ! Highest TOP
       do K = 2, KX
          if(Z0G(I,K) - Z0G(I,1) .lt. MAXZPAR ) KTOP=K
       end do

       do K=1,KX
          RLIQ(K) = 0.
          RICE(K) = 0.
       end do


       DXSQ = DXDY(I)
10     P300 = 1000.*(PSB(I)*SIGKFC(I,KL)-30.)



       ML = 0
       !VDIR NOLSTVAL
       do K = 1, KX

          ! NK = KX-K+1

          PMID = 0.5 * ( 1000.*PSB(I) + 100.E2 )
          if (PP0(I,K).ge.PMID) L5 = K
          if (PP0(I,K).ge.P300) LLFC = K
          if (TT0(I,K).gt.TRPL) ML=K
       enddo


       ! 6. TRIGGER FUNCTION
       ! ===================

       KMIX=1
25     LOW = KMIX

       !  No levels (or parcels) in the lowest
       !  300 mb were found to be unstable
       !  ==> goto the next grid point


       if (LOW.gt.LLFC .or. Z0G(I,LOW)-Z0G(I,1) >= MAXZPAR) goto 325
       LC=LOW
       MXLAYR=0


       ! Find the number of layers needed to have a mixed layer of at least 60 mb

       ! Note: in a future version of the scheme,
       ! why not use results from the boundary-layer
       ! scheme concerning the well-mixed layer ?
       ! DPthmixg(i) = pressure depth of the mixed layer
       ! thmixg(i)   = potential temperature of the mixed layer
       ! NLAYRS  = number of model layers in the mixed layer

       NLAYRS=0
       DPTHMXG(I)=0.
       do NK = LC,KX
          DPTHMXG(I)=DPTHMXG(I)+DPP(I,NK)
          NLAYRS=NLAYRS+1
          if(DPTHMXG(I).gt.DPMIX)goto 64   !The layer is at least 60-mb deep
       enddo

       goto 325

64     KPBL=LC+NLAYRS-1
       KMIX=LC+1
       THMIXG(I)=0.
       QMIXG(I)=0.
       TKEMIXG(I)=0.
       ZMIX=0.
       PMIXG(I)=0.
       DPTHMXG(I)=0.


       ! Integrated properties

       do  NK = LC,KPBL
          DPTHMXG(I)=DPTHMXG(I)+DPP(I,NK)
          ROCPQ=0.2854*(1.-0.28*Q00(I,NK))
          THMIXG(I)=THMIXG(I)+DPP(I,NK)*TT0(I,NK)*(P00/PP0(I,NK))**ROCPQ
          QMIXG(I)=QMIXG(I)+DPP(I,NK)*Q00(I,NK)
          TKEMIXG(I)=TKEMIXG(I)+DPP(I,NK)*TKE0(I,NK)
          ZMIX=ZMIX+DPP(I,NK)*Z0G(I,NK)
          PMIXG(I)=PMIXG(I)+DPP(I,NK)*PP0(I,NK)
       enddo


       ! Average to obtain the mixed properties

       THMIXG(I)=THMIXG(I)/DPTHMXG(I)
       QMIXG(I)=QMIXG(I)/DPTHMXG(I)
       QMIXG(I)=AMAX1( QMIXG(I),1.0E-10 )
       TKEMIXG(I)=TKEMIXG(I)/DPTHMXG(I)
       TKEMIXG(I)=max(TKEMIXG(I),ETRMIN)
       ZMIX=ZMIX/DPTHMXG(I)
       PMIXG(I)=PMIXG(I)/DPTHMXG(I)
       ROCPQ=0.2854*(1.-0.28*QMIXG(I))
       TMIXG(I)=THMIXG(I)*(PMIXG(I)/P00)**ROCPQ
       THVMIXG(I)=THMIXG(I)*(1.+0.608*QMIXG(I))


       ! We have all the characteristics of the
       ! departure parcels (the mixed values).
       ! Now we need the temperature and pressure
       ! at the LCL.  First calculate dew point
       ! temperature TDPT by inverting the saturated
       ! water vapor equation (see Davis and Jones 1983).

       EMIX=QMIXG(I)*PMIXG(I)/(0.622+QMIXG(I))
       TLOG=ALOG(EMIX/ALIQ)
       TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
       TLCLG(I)=TDPT-(.212+1.571E-3*(TDPT-TRPL)-4.36E-4*(TMIXG(I)-TRPL))* &
            (TMIXG(I)-TDPT)
       TLCLG(I)=AMIN1(TLCLG(I),TMIXG(I))
       TVLCL= TLCLG(I)*(1.+0.608*QMIXG(I))
       CPORQ=1./ROCPQ
       PLCLG(I)=P00*(TLCLG(I)/THMIXG(I))**CPORQ



       ! Find the first level above the LCL
       cin = 0.
       do NK = LC,KL
          KLCL=NK
          if(PLCLG(I).ge.PP0(I,NK))goto 35
          if (nk > kpbl) then
             dzz = z0g(i,nk+1)-z0g(i,nk)
             udlbe = (2.*thvmixg(i)/(thv0(i,nk+1)+thv0(i,nk))-1.)*dzz
             if (udlbe < 0.) cin = cin - udlbe*GRAV
          endif
       enddo

       goto 325

35     continue
       KLCLM1=KLCL-1
       DLP=ALOG(PLCLG(I)/PP0(I,KLCLM1))/ALOG(PP0(I,KLCL)/PP0(I,KLCLM1))


       ! Estimate environmental temperature
       ! and mixing ratio at the LCL


       TENV=TT0(I,KLCLM1)+(TT0(I,KLCL)-TT0(I,KLCLM1))*DLP
       QENV=Q00(I,KLCLM1)+(Q00(I,KLCL)-Q00(I,KLCLM1))*DLP
       TVEN= TENV*(1.+0.608*QENV)
       TVBAR=0.5*(TV00(I,KLCLM1)+TVEN)

       ! Trigger velocity (represents subgrid
       ! variance) at the LCL
       if (associated(wklclp)) then
          wklcla(i) = wklcl0(i,KLCLM1) + (wklcl0(i,klcl)-wklcl0(i,KLCLM1))*dlp
       endif


       ! to avoid getting negative dzz later on:
       ! use simple linear interpolation to get
       ! the level zlcl

       !  ZLCL=Z0G(I,KLCLM1)+RGASD*TVBAR*ALOG(PP0(I,KLCLM1)/PLCLG(I))/GRAV

       ZLCL= Z0G(I,KLCLM1)+ (Z0G(I,KLCLM1+1) - Z0G(I,KLCLM1)) * &
            ( PLCLG(I)-PP0(I,KLCLM1) )/ ( PP0(I,KLCLM1+1) - PP0(I,KLCLM1) )

       ! we do not want zlcl to be at the level
       ! z0g(KLCLM1+1), otherwise we get dzz=0.
       ! Therefore we impose a difference
       ! of 1 metre.

       ZLCL= min (ZLCL, Z0G(I,KLCLM1+1) - 1.)


       ! Check to see if the parcel is buoyant

       select case (deeptrig(i))
       case ('ORIGINAL')
          WKLCLD=WKLCLA(I)
          wklclout(i) = wklcld
          WKL = WW0(I,KLCLM1)+(WW0(I,KLCL)-WW0(I,KLCLM1))*DLP-WKLCLD
          WABS = abs(WKL)+1.E-10
          WSIGNE = WKL/WABS
          DTLCL = 4.64*WSIGNE*WABS**0.33
       case ('TKE')
          dzz = z0g(i,klcl)-z0g(i,KLCLM1)
          dtenv = (thv0(i,klcl)-thv0(i,KLCLM1))/dzz
          ktlcl = kt0(i,KLCLM1)+(kt0(i,klcl)-kt0(i,KLCLM1))*dlp
          dtlcl = min(10.*ktlcl*dtenv/sqrt(tkemixg(i)),3.)
       case ('KAIN04')
          wkl = 100.*(ww0(i,KLCLM1)+(ww0(i,klcl)-ww0(i,KLCLM1))*dlp) - &
               2.*min(zlcl/2000.,1.)  !cm/s
          dtlcl = sign(1.*(abs(wkl))**(1./3.),wkl)
       case ('SCALED')
          dtlcl = tstar(i)
       case DEFAULT
          call physeterror('kfrpn', 'no trigger type '//trim(deeptrig(i)))
          return
       end select

       ! Parcel is convectively unstable
       ! Calculate the tendencies

       if (TLCLG(I)+DTLCL.gt.TENV) goto 45                                            !Initiation of new cloud
       if (deep_cloudobj .and. cowlcl(i) > WU_MIN) goto 45  !Pre-existing cloud triggering


       !  Parcel not buoyant.
       !  goto the next column.

       if(KPBL.ge.LLFC)goto 325
       goto 25


       ! 7. CHARACTERISTICS OF THE BUOYANT PARCELS AT THE LCL
       ! ====================================================

       !  THETEU = equivalent potential temperature
       !           of the updraft just below the LCL.

       !  Note that K here is the level just below
       !  the LCL, whereas KLCL is the level just
       !  above it.

45     THETEU(KLCLM1)=TMIXG(I)*(1.E5/PMIXG(I))**(0.2854*(1.-0.28*QMIXG(I)))* &
            exp((3374.6525/TLCLG(I)-2.5403)*QMIXG(I)*(1.+0.81*QMIXG(I)))
       ES=ALIQ*exp((TENV*BLIQ-CLIQ)/(TENV-DLIQ))
       ES=min( ES , 0.5*PLCLG(I) )


       ! Adjust the pressure at the LCL using the
       ! height of the LCL and an average temperature
       ! between the environmental and grid-scale
       ! temperature.

       TVAVG=0.5*(TV00(I,KLCL)+ TENV*(1.+0.608*QENV))
       PLCLG(I)=PP0(I,KLCL)*exp(GRAV/(RGASD*TVAVG)*(Z0G(I,KLCL)-ZLCL))


       ! Vertical motion of the parcel and equivalent
       ! pot. temperature of the environment at the
       ! LCL.

       QESE=0.622*ES/(PLCLG(I)-ES)
       QESE=max( min( QESE , 0.050 ) , 1.E-6 )
       select case (deeptrig(i))
       case ('ORIGINAL')
          GDT=GRAV*DTLCL*(ZLCL-Z0G(I,LC))/(TV00(I,LC)+TVEN)
          WLCL=1.+.5*WSIGNE*sqrt(abs(GDT)+1.E-10)
          if (deep_cloudobj) then

!!! TODO: compute an updraft velocity here that is a linear decrease
!!! below cloud base height (cozlcl) to proposed cloud base (zlcl) or
!!! is the cloud base updraft (cowlcl) above.  Use this instead of
!!! cowlcl here so that initiation with lower cloud bases isn't 
!!! artificially enhanced with the updraft velocity from remnant
!!! convection at a higher level.  (This is consistent with the linear
!!! decrease in UMF below cloud base.)

             if (cowlcl(i) > WLCL) then
                wlcl = cowlcl(i)  !Existing updraft is stronger
             else
                coage(i) = 0.     !New initiation is stronger
             endif
          endif
       case ('TKE')
          wkl = ww0(i,KLCLM1)+(ww0(i,klcl)-ww0(i,KLCLM1))*dlp
!!$         stdw = sqrt(tkemixg(i))
!!$         int_max = 100
!!$         dwz = (wkl+5.*stdw)/(int_max-1)
!!$         wlcl = 0.
!!$         do il=1,int_max
!!$            wz = il * dwz
!!$            wlcl = wlcl + wz*exp(-0.5*((wz-wkl)/stdw)**2) * dwz
!!$         enddo
!!$         wlcl = wlcl / (2.*sqrt(PI)*stdw)
!!$         if (wlcl < 0.1) wlcl = 0.  !at least 10 cm/s updraft ascent
          wlcl = wkl + 1.*sqrt(tkemixg(i))
       case ('KAIN04')
          gdt = (zlcl-z0g(i,lc))*dtlcl/tven
          wlcl = 1. + 1.1*sqrt(max(gdt,0.))
       case ('SCALED')
          wlcl = wstar(i)
       case DEFAULT
          call physeterror('kfrpn', 'no trigger type '//trim(deeptrig(i)))
          return
       end select
       THTES(KLCLM1)=TENV*(1.E5/PLCLG(I))**(0.2854*(1.-0.28*QESE))* &
            exp((3374.6525/TENV-2.5403)*QESE*(1.+0.81*QESE))
       WTW=WLCL*WLCL

       ! No upward motion ==> take next slab
       ! and test for convective instability.

       if(WLCL.le.0.) goto 25
       TVLCL=TLCLG(I)*(1.+0.608*QMIXG(I))
       RHOLCL=PLCLG(I)/(RGASD*TVLCL)

       LCL=KLCL
       LET = LCL


       ! 8. COMPUTE UPDRAFT PROPERTIES
       ! =============================

       ! Initial updraft mass flux (VMFU)

       ! AU0 = area of the updraft at cloud base
       ! UMF = updraft mass flux
       ! VMFLCL = vertical mass flux at the LCL
       ! RATIO2 = degree of glaciation

       WU(KLCLM1)=WLCL
       AU0=PI*RAD*RAD
       UMF(KLCLM1)=RHOLCL*AU0
       VMFLCL=UMF(KLCLM1)
       UPOLD=VMFLCL
       UPNEW=UPOLD
       RATIO2(KLCLM1)=0.


       ! UER = environmental entrainment rate
       ! ABE = the available buoyant energy
       ! TRPPT = total rate of precip. production


       UER(KLCLM1)=0.
       ABE=0.
       CAPE=0.
       TRPPT=0.

       ! TU, TVU, QU = updraft properties

       TU(KLCLM1)=TLCLG(I)
       TVU(KLCLM1)=TVLCL
       QU(KLCLM1)=QMIXG(I)

       ! EQFRC = fraction of environmental air in mixed subparcels

       EQFRC(KLCLM1)=1.

       !  RLIQ = liquid water in the updraft
       !  RICE = ice in the updraft

       RLIQ(KLCLM1)=0.
       RICE(KLCLM1)=0.

       !  QLQOUT = liquid water fallout from the updraft
       !  QICOUT = ice fallout from the updraft

       QLQOUT(KLCLM1)=0.
       QICOUT(KLCLM1)=0.

       !  DETLQ = detrained liquid water from the updraft
       !  DETIC = detrained ice from the updraft

       DETLQ(KLCLM1)=0.
       DETIC(KLCLM1)=0.

       !  PPTLIQ = liquid precipitation from conv. updraft
       !  PPTICE = solid precipitation from conv. updraft

       PPTLIQ(KLCLM1)=0.
       PPTICE(KLCLM1)=0.


       !  KFRZ = freezing level

       IFLAG = 0
       KFRZ=LC


       ! The amount of convective available potential
       ! energy (CAPE) is calculated with respect to
       ! an undiluted parcel ascent.
       ! THUDL = equivalent potential temperature of an
       ! undiluted parcel.
       ! TUDL = temperature of an undiluted parcel

       THTUDL=THETEU(KLCLM1)
       TUDL=TLCLG(i)

       ! TTEMP is used during the calculations of the
       ! linear glaciation process.
       ! This temperature is initially set to the
       ! temperature at which freezing is specified
       ! to occur within the glaciation interval.
       ! It is set equal to the UPDR temperature
       ! at the previous model level.

       TTEMP=TTFRZ

       ! LOOP FOR UPDRAFT CALCULATIONS.  We calculate
       ! here the updraft temperature, the mixing ratio,
       ! the vertical mass flux, the lateral detrainment
       ! of mass and moisture, and the precipitation
       ! rates at each model level.
       ! This loop is also pretty long...

       MAIN_UPDRAFT_LOOP60:  do NK=KLCLM1,KLM
          NK1=NK+1
          RATIO2(NK1)=RATIO2(NK)

          ! At a first stage, consider that THETEU,
          ! QU, RLIQ, and RICE, are conserved from
          ! one level to another (UNDILUTED ascent
          ! before the role of ENTRAINMENT is evaluated).

          FRC1=0.
          TU(NK1)=TT0(i,nk1)
          THETEU(NK1)=THETEU(NK)
          QU(NK1)=QU(NK)
          RLIQ(NK1)=RLIQ(NK)
          RICE(NK1)=RICE(NK)


          ! Call the mixing subroutine for the
          ! UNDILUTED ASCENT.

          call TPMIX(PP0(I,NK1),THETEU(NK1),TU(NK1),QU(NK1),RLIQ(NK1), &
               RICE(NK1),QNEWLQ,QNEWIC,RATIO2(NK1),RL,XLV0,XLV1,XLS0,XLS1, &
               ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

          TVU(NK1)=TU(NK1)*(1.+0.608*QU(NK1))

          ! Check to see if updraft temperature is
          ! within the freezing interval.  If it is,
          ! calculate the fractional conversion to
          ! glaciation and adjust QNEWLQ to reflect the
          ! gradual change in THETAU since the last
          ! model level.  The glaciation effects will be
          ! determined after the amount of condensate
          ! available after precipitation fallout is
          ! determined.  TTFRZ is the temperature at
          ! which glaciation begins, TBFRZ is the
          ! temperature at which it ends.

          ! FRC1 = fractional conversion to glaciation
          !        (considering the entire glaciation process,
          !         from TTFRZ to TBFRZ)
          ! R1   = fractional conversion to glaciation
          !        (but this time considering only the
          !         remaining glaciation from TU(NK) to
          !         TBFRZ).

          if(TU(NK1).le.TTFRZ.and.IFLAG.lt.1)then
             if(TU(NK1).gt.TBFRZ)then

                !   within the glaciation buffer zone...

                if(TTEMP.gt.TTFRZ)TTEMP=TTFRZ
                FRC1=(TTEMP-TU(NK1))/(TTFRZ-TBFRZ)
                R1=(TTEMP-TU(NK1))/(TTEMP-TBFRZ)
             else

                !   glaciation processes are over...

                FRC1=(TTEMP-TBFRZ)/(TTFRZ-TBFRZ)
                R1=1.
                IFLAG=1
             endif

             !   In this glaciation process, the new liquid
             !   generated in the updraft directly goes into
             !   QNEWFRZ.

             QNWFRZ=QNEWLQ

             !   Increase new ice and decrease new liquid
             !   following the R1 fractional glaciation

             QNEWIC=QNEWIC+QNEWLQ*R1*0.5
             QNEWLQ=QNEWLQ-QNEWLQ*R1*0.5
             EFFQ=(TTFRZ-TBFRZ)/(TTEMP-TBFRZ)
             TTEMP=TU(NK1)
          endif

          !  Updraft vertical velocity and
          !  precipitation fallout.

          if(NK.eq.KLCLM1)then

             !  For the level just below the LCL

             BE=(TVLCL+TVU(NK1))/(TVEN+TV00(I,NK1))-1.
             BOTERM=2.*(Z0G(I,NK1)-ZLCL)*GRAV*BE/1.5
             ENTERM=0.
             DZZ=Z0G(I,NK1)-ZLCL

          else

             !  And above...

             BE=(TVU(NK)+TVU(NK1))/(TV00(I,NK)+TV00(I,NK1))-1.
             BOTERM=2.*DZP(I,NK)*GRAV*BE/1.5
             ENTERM=2.*UER(NK)*WTW/UPOLD
             DZZ=DZP(I,NK)
          endif

          !  Condensation in the updraft.

          call CONDLOAD_SAFE(RLIQ(NK1),RICE(NK1),WTW,DZZ,BOTERM,ENTERM, &
               ACRATE(I),QNEWLQ,QNEWIC,QLQOUT(NK1),QICOUT(NK1), &
               GRAV)

          WABS=sqrt(max(abs(WTW),1.E-10))

          WU(NK1)=WTW/WABS

          !   If the vertical velocity is less
          !   than zero, exit the updraft loop
          !   and if the cloud is deep enough,
          !   finalize the updraft calculations.

          if(WU(NK1).le.0.) exit MAIN_UPDRAFT_LOOP60


          !   Update ABE for undilute ascent.

          !   Note that we still did not consider
          !   any entrainment effect -
          !   except for the calculation of w.

          THTES(NK1)=TT0(I,NK1)*(1.E5/PP0(I,NK1))**(0.2854*(1.-0.28*QST1(I,NK1)))* &
               exp((3374.6525/TT0(I,NK1)-2.5403)*QST1(I,NK1)*(1.+0.81*QST1(I,NK1)))
          UDLBE=((2.*THTUDL)/(THTES(NK)+THTES(NK1))-1.)*DZZ
          if(UDLBE.gt.0.) ABE = ABE + UDLBE*GRAV

          !   Update CAPE and CIN for undilute ascent
          if (be > 0.) then
             cape = cape + GRAV*be*dzz
          else
             cin = cin - GRAV*be*dzz
          endif

          !   Other effects of cloud
          !   glaciation if within the specified
          !   temperature interval
          !   (on temperature, water variables, etc...).

          !   Phase change effects are calculated in
          !   DTFRZNEW subroutine.

          if(FRC1.gt.1.E-6)then

             call DTFRZNEW2(TU(NK1),PP0(I,NK1),THETEU(NK1),QU(NK1),RLIQ(NK1), &
                  RICE(NK1),RATIO2(NK1),QNWFRZ,RL,FRC1,EFFQ, &
                  IFLAG,XLV0,XLV1,XLS0,XLS1, &
                  ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

          endif



          !  Call subroutine to calculate ENVIRONMENTAL
          !  potential temperature within glaciation
          !  interval.  THETAE must be calculated with
          !  respect to the same degree of glaciation for
          !  all entraining air.

          call ENVIRTHT(PP0(I,NK1),TT0(I,NK1),Q00(I,NK1), &
               THETEE(NK1),RATIO2(NK1),RL, &
               ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

          !  Trace the level of minimum eq. potential
          !  temperature for the downdraft.

          if(NK1.eq.KLCL)THTMIN=THTES(NK1)
          THTMIN=AMIN1(THTES(NK1),THTMIN)
          if(THTMIN.eq.THTES(NK1))KMIN=NK1



          !  BEGINNING OF ENTRAINMENT/DETRAINMENT CALCULATIONS


          !  REI is the rate of environmental inflow

          REI=VMFLCL*DPP(i,NK1)*0.03/RAD
          TVQU(NK1)=TU(NK1)*(1.+0.608*QU(NK1)-RLIQ(NK1)-RICE(NK1))


          !  UER   = updraft entrainment rate
          !  UDR   = updraft detrainment rate
          !  EQFRC = critical fraction of environmental
          !          air in the mixed subparcels at which
          !          they are neutrally stable.

          !  If cloud parcels are virtually colder than
          !  the environment, no entrainment is allowed
          !  at this level.

          if(TVQU(NK1).le.TV00(I,NK1))then
             UER(NK1)=0.0
             UDR(NK1)=REI
             EE2=0.
             UD2=1.
             EQFRC(NK1)=0.
             goto 55
          endif

          LET=NK1
          TTMP=TVQU(NK1)

          !  Determine the critical mixed fraction of
          !  updraft and environmental air.

          !  First suppose a subparcel with 5% updraft air
          !  and 95% environmental air.

          F1=0.95
          F2=1.-F1

          !  Properties of the subparcel

          !  THTMP = eq. pot. temperature
          !  RTMP  = specific humidity
          !  TMPLIQ = liquid water content
          !  TMPICE = ice content

          THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)
          RTMP=F1*Q00(I,NK1)+F2*QU(NK1)
          TMPLIQ=F2*RLIQ(NK1)
          TMPICE=F2*RICE(NK1)

          !  Find the temperature and new liquid for
          !  this temporary subparcel.

          call TPMIX(PP0(I,NK1),THTTMP,TTMP,RTMP,TMPLIQ,TMPICE,QNEWLQ, &
               QNEWIC,RATIO2(NK1),RL,XLV0,XLV1,XLS0,XLS1, &
               ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

          !  moist virtual temperature of the subparcel

          TU95=TTMP*(1.+0.608*RTMP-TMPLIQ-TMPICE)

          !  If the parcel is buoyant with these mixing
          !  proportions, then the critical fraction of
          !  environmental air is close to 1.
          !  (In fact, we set it to 1. in this case).
          !  No detrainment in this case, only entrainment.

          if(TU95.gt.TV00(I,NK1))then
             EE2=1.
             UD2=0.
             EQFRC(NK1)=1.0
             goto 50
          endif

          !  Calculate the critical fraction of environmental
          !  air by creating a supbarcel with 10% updraft and
          !  90* environmental air.

          F1=0.10
          F2=1.-F1
          THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)
          RTMP=F1*Q00(I,NK1)+F2*QU(NK1)
          TMPLIQ=F2*RLIQ(NK1)
          TMPICE=F2*RICE(NK1)

          call TPMIX(PP0(I,NK1),THTTMP,TTMP,RTMP,TMPLIQ,TMPICE,QNEWLQ, &
               QNEWIC,RATIO2(NK1),RL,XLV0,XLV1,XLS0,XLS1, &
               ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

          TU10=TTMP*(1.+0.608*RTMP-TMPLIQ-TMPICE)

          ! It is possible, from this subparcel, to
          ! evaluate the critical fraction EQFRC.

          ! tvdiff can be zero when updraft and env temperature are very close
          TVDIFF = TU10-TVQU(NK1)
          if(TVDIFF.ge.0.0) then
             EQFRC(NK1)=0.0
          else
             EQFRC(NK1)=(TV00(I,NK1)-TVQU(NK1))*F1/(min(TVDIFF,-eps))
             EQFRC(NK1)=AMAX1(0.0,EQFRC(NK1))
             EQFRC(NK1)=AMIN1(1.0,EQFRC(NK1))
          endif

          ! In the simple cases in which the subparcel
          ! mixing possibilities are all buoyant or all
          ! non-buoyant, no need to integrate the
          ! Gaussian distribution.

          if(EQFRC(NK1).eq.1)then

             ! All buoyant case

             EE2=1.
             UD2=0.
             goto 50
          elseif(EQFRC(NK1).eq.0.)then

             ! All non-buoyant case

             EE2=0.
             UD2=1.
             goto 50
          else

             ! Subroutine PROF5 integrates over the
             ! Gaussian distribution to determine the
             ! fractional entrainment and detrainment
             ! rates.

             call PROF5(EQFRC(NK1),EE2,UD2)

          endif

          ! EE1 and UD1 are the entrainment and detrainment
          ! rates from the previous level.

50        if(NK.eq.KLCLM1)then
             EE1=1.
             UD1=0.
          endif

          !  Net entrainment and detrainment rates are
          !  given by the average fractional values in
          !  the layer.

          UER(NK1)=0.5*REI*(EE1+EE2)
          UDR(NK1)=0.5*REI*(UD1+UD2)

          !  Modification of the updraft detrainment rate
          !  to include the effect of an "ensemble" of
          !  updrafts.  This is to avoid the problem of
          !  intense detrainment at a single layer at the
          !  cloud top level (between the LET and the top
          !  in fact).

!!$        if ( Z0G(I,NK1).gt.ZLCL+(1.-DETREG)*(ZTOP(I)-ZLCL) ) then
!!$          DETFRAC  = DETTOT / max( FLOAT( KTOP-NK1 + 1 ), 1. )
!!$          DETFRAC  = min( DETFRAC, 0.25 )
!!$          UDR(NK1) = max( UDR(NK1), DETFRAC*UMF(NK) )
!!$          UDR(NK1) = min( UDR(NK1), UMF(NK) - 11. )
!!$        end if

          !  If the calculated updraft detrainment rate
          !  is greater than the total updraft mass flux,
          !  all cloud mass detrains.  Exit updraft
          !  calculations.

55        if(UMF(NK)-UDR(NK1).lt.10.)then
             !  If the calculated detrained mass flux is
             !  greater than the total updraft flux, impose
             !  total detrainment of updraft mass at the
             !  previous model level.

             if (UDLBE > 0.) ABE=ABE-UDLBE*GRAV
             if (be > 0.) then
                cape = cape - GRAV*be*dzz
             else
                cin = cin + GRAV*be*dzz
             endif
             LET=NK
             exit MAIN_UPDRAFT_LOOP60
          endif


          !  New updraft mass flux (remove the part that
          !  was detrained - and add the part that was
          !  entrained).

          EE1=EE2
          UD1=UD2
          UPOLD=UMF(NK)-UDR(NK1)
          UPNEW=UPOLD+UER(NK1)
          UMF(NK1)=UPNEW

          ! DETLQ and DETIC are the rates of detrainment
          ! of liquid and ice in the detraining updraft mass

          DETLQ(NK1)=RLIQ(NK1)*UDR(NK1)
          DETIC(NK1)=RICE(NK1)*UDR(NK1)
          QDT(NK1)=QU(NK1)

          ! New properties of the updraft after
          ! entrainment processes

          QU(NK1)=(UPOLD*QU(NK1)+UER(NK1)*Q00(I,NK1))/UPNEW
          THETEU(NK1)=(THETEU(NK1)*UPOLD+THETEE(NK1)*UER(NK1))/UPNEW
          RLIQ(NK1)=RLIQ(NK1)*UPOLD/UPNEW
          RICE(NK1)=RICE(NK1)*UPOLD/UPNEW


          ! END OF ENTRAINMENT/DETRAINMENT CALCULATIONS


          ! KFRZ is the highest model level at which
          ! liquid condensate is generated.  PPTLIQ is the
          ! rate of generation (fallout) of liquid precipitation
          ! at a given model level.  PPTICE is the same for ice.
          ! TRPPT is the total rate of production of precip.
          ! up to the current model level.

          if(abs(RATIO2(NK1)-1.).gt.1.E-6)KFRZ=NK1
          PPTLIQ(NK1)=QLQOUT(NK1)*UMF(NK)
          PPTICE(NK1)=QICOUT(NK1)*UMF(NK)
          TRPPT=TRPPT+PPTLIQ(NK1)+PPTICE(NK1)
          if(NK1.le.KPBL)UER(NK1)=UER(NK1)+VMFLCL*DPP(i,NK1)/DPTHMXG(I)

       enddo MAIN_UPDRAFT_LOOP60


       !  END OF THE PRETTY LONG LOOP for
       !  the updraft calculations.

       !  Check cloud depth.  If the cloud is deep
       !  enough, estimate the equilibrium temperature
       !  level (ETL) and adjust the mass flux profile
       !  at cloud top so that the mass flux decreases
       !  to zero as a linear function of pressure
       !  difference the LET and cloud top.

       !  LTOP is the model level just below the level
       !  at which vertical velocity first becomes
       !  negative.

       LTOP=NK

       ! For the TKE trigger, compute energy available to overcome CIN
!!$      if (trigger_type == 'tke') then
!!$         trigger_energy = 0.
!!$         do il=1,int_max
!!$            wz = il * dwz
!!$            trigger_energy = trigger_energy + exp(-0.5*((wz-wkl)/stdw)**2) * dwz
!!$         enddo
!!$         trigger_energy = tkemixg(i) * trigger_energy / (2.*sqrt(PI)*stdw)
!!$      else
       trigger_energy = cin + 1. !always pass this part of trigger
!!$      endif

       CLDHGT=Z0G(I,LTOP)-ZLCL

       !  If the cloud top height is less than the
       !  specified minimum height, go back to the
       !  next highest 50 mb layer to see if a bigger
       !  cloud can be obtained from the parcel.

       if(CLDHGT < CDEPTH .or. ABE < 1. .or. trigger_energy < cin)then
          do NK=KLCLM1,LTOP
             UMF(NK)=0.
             UDR(NK)=0.
             UER(NK)=0.
             DETLQ(NK)=0.
             DETIC(NK)=0.
             PPTLIQ(NK)=0.
             PPTICE(NK)=0.
          enddo
          goto 25
       endif

       KKFC(I)=1.  !PV:if we made it here, the column is active

       ZBASEOUT(I) = ZLCL
       ZTOPOUT(I)  = Z0G(I,LTOP)

       ! ABOVE THE LET

       ! If the LET and LTOP are the same, detrain all
       ! of the updraft mass at this level.

       if(LET.eq.LTOP)then
          UDR(LTOP)=UMF(LTOP)+UDR(LTOP)-UER(LTOP)
          DETLQ(LTOP)=RLIQ(LTOP)*UDR(LTOP)*UPNEW/UPOLD
          DETIC(LTOP)=RICE(LTOP)*UDR(LTOP)*UPNEW/UPOLD
          UER(LTOP)=0.
          UMF(LTOP)=0.
          goto 85
       endif


       ! Begin total detrainment at the level above the LET.

       DPTT=0.
       do NJ=LET+1,LTOP
          DPTT=DPTT+DPP(I,NJ)
       enddo
       DUMFDP=UMF(LET)/DPTT


       ! Adjust mass flux profiles, detrainment rates,
       ! and precipitation fall rates to reflect the
       ! linear decrease in mass flux between the LET
       ! and cloud top.

       do NK=LET+1,LTOP
          ! A small part is detrained at each level above the LET.

          UDR(NK)=DPP(I,NK)*DUMFDP
          UMF(NK)=UMF(NK-1)-UDR(NK)
          DETLQ(NK)=RLIQ(NK)*UDR(NK)
          DETIC(NK)=RICE(NK)*UDR(NK)

          if(NK.ge.LET+2)then
             TRPPT=TRPPT-PPTLIQ(NK)-PPTICE(NK)
             PPTLIQ(NK)=UMF(NK-1)*QLQOUT(NK)
             PPTICE(NK)=UMF(NK-1)*QICOUT(NK)
             TRPPT=TRPPT+PPTLIQ(NK)+PPTICE(NK)
          endif
       enddo


85     continue


       ! Below departure level (LC:base of mixed layer)

       do NK=1,LC-1
          TU(NK)=0.
          QU(NK)=0.
          TVU(NK)=0.
          UMF(NK)=0.
          WU(NK)=0.
          UER(NK)=0.
       enddo

       ! Between the departure and cloud base levels (only entrainment)

       UMF(LC)=VMFLCL*DPP(i,LC)/DPTHMXG(I)        ! departure level
       UER(LC)=UMF(LC)                            ! departure level
       do NK=LC+1,KPBL                            ! Between the departure and top of mixed layer: linear increase
          UER(NK)=VMFLCL*DPP(i,NK)/DPTHMXG(I)     
          UMF(NK)=UMF(NK-1)+UER(NK)
       enddo
       do NK=KPBL+1,KLCLM1                        ! Between top of mixed layer and cloud base: constant umf
          UMF(NK)=VMFLCL
          UER(NK)=0.
       enddo
       do NK=LC,KLCLM1
          TU(NK)=TMIXG(I)+(Z0G(I,NK)-ZMIX)*GDRY
          QU(NK)=QMIXG(I)
          TVU(NK)= TU(NK)*(1.+0.608*QU(NK))
          WU(NK)=WLCL
       enddo
       !   Between first level and one level below cloud base(KLCLM1=KLCL-1)
       do NK=1,KLCLM1
          UDR(NK)=0.
          QDT(NK)=0.
          UU(NK)=0.
          VU(NK)=0.
          RLIQ(NK)=0.
          RICE(NK)=0.
          QLQOUT(NK)=0.
          QICOUT(NK)=0.
          PPTLIQ(NK)=0.
          PPTICE(NK)=0.
          DETLQ(NK)=0.
          DETIC(NK)=0.
          RATIO2(NK)=0.
          EE=Q00(I,NK)*PP0(I,NK)/(0.622+Q00(I,NK))
          TLOG=ALOG(EE/ALIQ)
          TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
          TSAT=TDPT-(.212+1.571E-3*(TDPT-TRPL)-4.36E-4*(TT0(I,NK)-TRPL))* &
               (TT0(I,NK)-TDPT)
          THTA=TT0(I,NK)*(1.E5/PP0(I,NK))**(0.2854*(1.-0.28*Q00(I,NK)))
          THETEE(NK)=THTA* &
               exp((3374.6525/TSAT-2.5403)*Q00(I,NK)*(1.+0.81*Q00(I,NK)))
          THTES(NK)=THTA* &
               exp((3374.6525/TT0(I,NK)-2.5403)*QST1(I,NK)*(1.+0.81*QST1(I,NK)))
          EQFRC(NK)=1.0
       enddo


       ! Some more definitions.

       LTOP1=LTOP+1
       LTOPM1=LTOP-1

       ! Make sure that the updraft vertical motion
       ! is not too small so that the cloud fraction
       ! becomes not too large just above the lifting
       ! condensation level (first level in the cloud
       ! above which WU is different from WLCL)

       do NK=LC+1,LTOPM1
          if (WU(NK).lt.WLCL) then
             WU(NK) = max(WLCL,WU(NK))
          else if (WU(NK).gt.WLCL) then
             goto 900
          endif
       end do
900    continue


       ! Define variables above cloud top.

       do NK = LTOP1, KX
          UMF(NK) = 0.
          UDR(NK)=0.
          UU(NK)  = 0.
          VU(NK)  = 0.
          UER(NK)=0.
          QDT(NK)=0.
          RLIQ(NK)=0.
          RICE(NK)=0.
          QLQOUT(NK)=0.
          QICOUT(NK)=0.
          DETLQ(NK)=0.
          DETIC(NK)=0.
          PPTLIQ(NK)=0.
          PPTICE(NK)=0.
          if(NK.gt.LTOP1)then
             TU(NK)=0.
             QU(NK)=0.
             TVU(NK)=0.
             WU(NK)=0.
          endif
          THTA0(NK)=0.
          THTAU(NK)=0.
          EMS(NK)=DPP(I,NK)*DXSQ/GRAV
          EMSD(NK)=1./EMS(NK)
          UG(NK) = U00(I,NK)
          VG(NK) = V00(I,NK)
          TG(NK) = TT0(I,NK)
          QG(NK) = Q00(I,NK)
          QLG(NK)= 0.
          QIG(NK)= 0.
          QTG(NK)= QG(NK)+QLG(NK)+QIG(NK)
          DTFM(NK)=0.

          OMGA(NK)=0.
          OMGAEFF(NK)=0.
          NUP(NK)=NK
          THADV(NK)=0.
          QADV(NK)=0.
          OMG(NK)=0.
       enddo
       OMG(KX+1)=0.
       P165=PP0(I,KLCL)-1.65E4
       ! PPTMLT=0.


       ! Initialize some other variables to be used later
       ! in the vertical advection scheme (from the
       ! surface to the cloud top).

       do NK=1,LTOP
          EMS(NK)=DPP(I,NK)*DXSQ/GRAV
          EMSD(NK)=1./EMS(NK)
          DTFM(NK)=0.

          EXN(NK)=(P00/PP0(I,NK))**(0.2854*(1.-0.28*QDT(NK)))
          THTAU(NK)=TU(NK)*EXN(NK)
          EXN(NK)=(P00/PP0(I,NK))**(0.2854*(1.-0.28*Q00(I,NK)))
          THTA0(NK)=TT0(I,NK)*EXN(NK)

          if(PP0(I,NK).gt.P165)LVF=NK
          !  PPTMLT=PPTMLT+PPTICE(NK)

          QTDT(NK) = QDT(NK)+RLIQ(NK)+RICE(NK) !updraft total water
          OMG(NK)=0.
       enddo

       CLVF=0.

       do NK = KLCL,LVF+1
          CLVF=CLVF+PPTLIQ(NK)+PPTICE(NK)
       enddo

       USR=UMF(LVF+1)*QU(LVF+1)+CLVF
       USR=AMIN1(USR,TRPPT)


       ! 9. CONVECTIVE TIME SCALE AND PRECIPITATION EFFICIENCY
       ! =====================================================


       !   Compute convective time scale (TIMEC).
       !   The mean wind at the LCL and midtroposphere !not true
       !   is used.                                    !not true

       WSPD(KLCL)=sqrt(U00(I,KLCL)*U00(I,KLCL)+V00(I,KLCL)*V00(I,KLCL))
       WSPD(L5)=sqrt(U00(I,L5)*U00(I,L5)+V00(I,L5)*V00(I,L5))
       WSPD(LTOP)=sqrt(U00(I,LTOP)*U00(I,LTOP)+V00(I,LTOP)*V00(I,LTOP))



       ! timec is the convective time scale, i.e. time over which we want to reduce CAPE by 90%
       ! timer is the time after which we want to recalculate(refresh) column conv activity and tendencies
       ! timer could be a fct of wind as intended in comments above
       ! timer should be min(timer,TIMEC)
       ! to refresh every timestep; set deep_timerefresh='1p';  will result in NIC = 1

       ! Compute the convective adjustment timescale (timec)
       if (deep_timeconv_sec > 0.) then
          nic = nint(deep_timeconv_sec/delt)
          timec = real(max(nic,1))*delt
       elseif (deep_timeconv == 'KFCTAUCAPE') then
          !AZ-PV Cape dependent timec; decrease timec for high cape values
          if(kfctaucape(1) > 0) then
             t1 = kfctaucape(1)     ! max timec
             t2 = kfctaucape(2)     ! min timec
             cmean = kfctaucape(3)  ! cape value  at which timec will be mean of t1 and t2
             dcape = kfctaucape(4)  ! decrease in timec from t1 to t2 will occur over range cmean-dcape to cmean+dcape
             tp = 0.5*(t2+t1)
             tm = 0.5*(t2-t1)
             dcape = max(1.,dcape)
             timec = tp + tm*tanh((2.*(abe-cmean)/dcape))
             nic   = nint(timec/delt)
             nic   = max(nic,1)
             timec = float(nic) * delt
          else
             call physeterror('kfrpn', 'kfctaucape must be provided for deep_timeconv='//trim(deep_timeconv))
             return
          endif
       elseif (deep_timeconv == 'BECHTOLD09') then
          if(kfctaucape(1) > 0) then
             t1 = kfctaucape(1)     ! max timec
             t2 = kfctaucape(2)     ! min timec
             wumean = 0.
             do k=klcl,let
                wumean = wumean + (wu(k)*dpp(i,k))/sum(dpp(i,klcl:let))
             enddo
             timec = 4.0 * (z0g(i,let)-zlcl) / (wumean)
             timec = min(t1,max(timec,t2))
             nic   = nint(timec/delt)
             nic   = max(nic,1)
             timec = float(nic) * delt
          else
             call physeterror('kfrpn','kfctaucape must be provided for deep_timeconv='//trim(deep_timeconv))
             return
          endif

       else
          call physeterror('kfrpn', trim(deep_timeconv)//' convective timescale option not implemented')
          return
       endif
       ! Compute the refresh timescale (timer) and resulting refresh count (nic)
       if (deep_timerefresh_sec > 0.) then
          timer = deep_timerefresh_sec
       elseif (deep_timerefresh == 'TIMECONV') then
          timer = timec
       else
          call physeterror('kfrpn', trim(deep_timeconv)//' refresh timescale option not implemented')
          return
       endif
       timer=min(timer,timec)
       nic = max(nint(timer/delt),1)

       !   Compute wind shear and precipitation
       !   efficiency.

       if(WSPD(LTOP).gt.WSPD(KLCL))then
          SHSIGN=1.
       else
          SHSIGN=-1.
       endif

       VWS=(U00(I,LTOP)-U00(I,KLCL))*(U00(I,LTOP)-U00(I,KLCL))+ &
            (V00(I,LTOP)-V00(I,KLCL))*(V00(I,LTOP)-V00(I,KLCL))
       VWS = 1.E3*SHSIGN*sqrt(VWS)/(Z0G(I,LTOP)-Z0G(I,LCL))
       PEF = 1.591+VWS*(-.639+VWS*(9.53E-2-VWS*4.96E-3))
       PEF = AMAX1(PEF,PEFMIN)
       PEF = AMIN1(PEF,PEFMAX)


       !   Precipitation efficiency as a
       !   function of cloud base.

       CBH = (ZLCL-Z0G(I,1))*3.281E-3

       if (CBH.lt.3.) then
          RCBH = .02
       else
          RCBH = .96729352+CBH*(-.70034167+CBH*(.162179896+CBH*(- &
               1.2569798E-2+CBH*(4.2772E-4-CBH*5.44E-6))))
       endif

       if (CBH.gt.25) RCBH = 2.4
       PEFCBH = 1./(1.+RCBH)
       PEFCBH=AMIN1(PEFCBH,.9)


       !   Mean precipitation efficiency used to
       !   compute rainfall.

       PEFF = .5*(PEF+PEFCBH)
       PEFFOUT(I) = PEFF


       ! 10. DOWNDRAFT PROPERTIES
       ! ========================


       ! Let the representative downdraft originate at
       ! the level of minimum saturation equivalent
       ! potential temperature (SEQT) in the cloud layer.
       ! Extend downward toward the surface, down to the
       ! level of downdraft buoyancy (LDB), where SEQT is
       ! less than min(SEQT) in the cloud layer.  Let downdraft
       ! detrain over a layer of specified pressure-depth (i.e., DPDD).

       TDER = 0.


       ! KMIN = level of minimum environmental saturation
       !        equivalent potential temperature

       ! If the level KMIN is higher than the cloud top
       ! ==> track again the level between LCL and LTOP.

       if(KMIN.ge.LTOP)then
          THTMIN=THTES(KLCL)
          KMIN=KLCL
          do NK=KLCL+1,LTOP-1
             THTMIN=AMIN1(THTMIN,THTES(NK))
             if(THTMIN.eq.THTES(NK))KMIN=NK
          enddo
       endif


       ! The LFS is defined as KMIN

       LFS=KMIN



       if(RATIO2(LFS).gt.0.) &
            call ENVIRTHT(PP0(I,LFS),TT0(I,LFS),Q00(I,LFS), &
            THETEE(LFS),0.,RL, &
            ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)


       ! The EQFRC fraction again represents the fraction
       ! of environmental air in mixed subparcels
       ! (from updraft and environmental) at the LFS.

       EQFRC(LFS)=(THTES(LFS)-THETEU(LFS))/(THETEE(LFS)-THETEU(LFS)+1.E-20)
       EQFRC(LFS)=AMAX1(EQFRC(LFS),0.)
       EQFRC(LFS)=AMIN1(EQFRC(LFS),1.)

       !  Eq. pot. temperature of the downdraft at the
       !  LFS (saturated).

       THETED(LFS)=THTES(LFS)


       !  Correction of the temperature due to melting
       !  effects.

       if(ML.gt.0)then
          DTMLTD=0.5*(QU(KLCL)-QU(LTOP))*CHLF/CPD
       else
          DTMLTD=0.
       endif

       TZ(LFS)=TT0(I,LFS)-DTMLTD


       !  Recalculate the eq. pot. temperature using
       !  this new temperature.

       ES = ALIQ*exp((TZ(LFS)*BLIQ-CLIQ)/(TZ(LFS)-DLIQ))
       QS = 0.622*ES/(PP0(I,LFS)-ES)
       QD(LFS)=EQFRC(LFS)*Q00(I,LFS)+(1.-EQFRC(LFS))*QU(LFS)
       THTAD(LFS)=TZ(LFS)*(P00/PP0(I,LFS))**(0.2854*(1.-0.28*QD(LFS)))

       if(QD(LFS).ge.QS)then
          THETED(LFS)=THTAD(LFS)* &
               exp((3374.6525/TZ(LFS)-2.5403)*QS*(1.+0.81*QS))
       else
          call ENVIRTHT(PP0(I,LFS),TZ(LFS),QD(LFS), &
               THETED(LFS),0.,RL, &
               ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
       endif



       !  Find the level of downdraft buoyancy (LDB)
       !  ==> level where THETAES for the downdraft
       !  at the LFS becomes larger then the environmental
       !  THETAES.

       LDB=1
       do NK=1,LFS-1
          ND=LFS-NK
          if(THETED(LFS).gt.THTES(ND) .or. ND.eq.1)then
             LDB=ND
             goto 110
          endif
       enddo


       !  Determine the LDT level (level where the
       !  downdraft detrains).  This level cannot be
       !  higher than the downdraft source level (LFS).
110    DPDDMX=kfcdpdd * dpddmult(i)
       DPT=0.
       NK=LDB
       LDT=-1
       !determined LDT must be less than LFS
       do while (LDT < 0 .and. NK < LFS)
          DPT=DPT+DPP(I,NK)
          if (DPT > DPDDMX) then
             LDT=NK
             DPDD=DPDDMX
          elseif (NK==LFS-1) then
             LDT=NK
             DPDD=DPT
          endif
          NK=NK+1
       enddo
       FRC=(DPDD+DPP(I,LDT)-DPT)/DPP(I,LDT)

       ! Properties of the downdrafts at the LFS:
       ! TVD   = moist virtual temperature
       ! RDD   = air density
       ! DMF   = downdraft mass flux
       ! DER   = downdraft entrainment rate
       ! DDR   = downdraft detrainment rate

       TVD(LFS)=TT0(I,LFS)*(1.+0.608*QST1(I,LFS))
       RDD=PP0(I,LFS)/(RGASD*TVD(LFS))
       A1=(1.-PEFF)*AU0
       DMF(LFS)=-A1*RDD
       DER(LFS)=EQFRC(LFS)*DMF(LFS)
       DDR(LFS)=0.


       !  Loop for the dowdraft mass flux

       DMF_LOOP140: do  ND=LFS-1,LDB,-1
          ND1=ND+1
          if(ND.le.LDT)then           ! below LDT; only detrainement
             DER(ND)=0.
             DDR(ND)=-DMF(LDT+1)*DPP(I,ND)*FRC/DPDD
             DMF(ND)=DMF(ND1)+DDR(ND)
             FRC=1.
             THETED(ND)=THETED(ND1)
             QD(ND)=QD(ND1)
             TZ(ND)=TPDD1(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
                  XLV0,XLV1, &
                  ALIQ,BLIQ,CLIQ,DLIQ)
             EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
             THTAD(ND)=TZ(ND)*EXN(ND)
          else                      ! above LDT; only entrainement

             DER(ND)=DMF(LFS)*0.03*DPP(I,ND)/RAD
             DDR(ND)=0.
             DMF(ND)=DMF(ND1)+DER(ND)

             !  Recalculate THETEE

             if(RATIO2(ND).gt.0.) &
                  call ENVIRTHT(PP0(I,ND),TT0(I,ND),Q00(I,ND), &
                  THETEE(ND),0.,RL, &
                  ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

             !  Downdraft properties deduced from properties
             !  at the previous level and entrained air.

             THETED(ND)=(THETED(ND1)*DMF(ND1)+THETEE(ND)*DER(ND))/DMF(ND)
             QD(ND)=(QD(ND1)*DMF(ND1)+Q00(I,ND)*DER(ND))/DMF(ND)
             TZ(ND)=TPDD1(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
                  XLV0,XLV1, &
                  ALIQ,BLIQ,CLIQ,DLIQ)
             EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
             THTAD(ND)=TZ(ND)*EXN(ND)
          endif
       enddo DMF_LOOP140

       TDER=0.


       ! Effect of downdraft reaching the LDB
       ! (cooling and moistening).  Relative humidity
       ! of 90% at this level.

       DETRN_LAYER: do ND=LDB,LDT
          TZ(ND)=TPDD1(PP0(I,ND),THETED(LDT),TT0(I,ND),QS,QD(ND),1.0, &
               XLV0,XLV1, &
               ALIQ,BLIQ,CLIQ,DLIQ)
          ES = ALIQ*exp((TZ(ND)*BLIQ-CLIQ)/(TZ(ND)-DLIQ))
          QS = 0.622*ES/(PP0(I,ND)-ES)
          DQSDT=(CLIQ-BLIQ*DLIQ)/((TZ(ND)-DLIQ)*(TZ(ND)-DLIQ))


          ! Adjust temperature and humidity for a
          ! relative humidity of 90%.

          RL=XLV0-XLV1*TZ(ND)
          DTMP=RL*QS*(1.-RHBC)/(cpd+RL*RHBC*QS*DQSDT)
          T1RH=TZ(ND)+DTMP
          ES=RHBC*ALIQ*exp((BLIQ*T1RH-CLIQ)/(T1RH-DLIQ))
          ES=min( ES , 0.5*PP0(I,ND) )
          QSRH=0.622*ES/(PP0(I,ND)-ES)
          QSRH=max( min( QSRH , 0.050 ) , 1.E-6 )


          ! Check to see if the mixing ratio at specified
          ! relative humidity is less than the actual
          ! mixing ratio.  If so, adjust to give zero
          ! evaporation.

          if(QSRH.lt.QD(ND))then
             QSRH=QD(ND)
             T1RH=TZ(ND)
          endif


          ! If this is not the case, use the adjusted
          ! values for temperature and humidity

          TZ(ND)=T1RH
          QS=QSRH
          TDER=TDER + (QS-QD(ND))*DDR(ND)
          QD(ND)=QS
          THTAD(ND)=TZ(ND)*(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
       enddo DETRN_LAYER

141    continue

       if (.not.DOWNDRAFT) TDER = 0.
       DOWNDRAFT_OFF: if (TDER < 1.) then

          ! No evaporation then no downdraft.

          ! PPTFLX = precipitation flux
          ! CPR    = condensation production rate
          ! TDER   = total downdraft evaporation rate (???)

          PPTFLX=TRPPT
          CPR=TRPPT
          TDER=0.
          CNDTNF=0.
          UPDINC=1.
          LDB=LFS
          do NDK=1,KX
             UD(NDK)=0.            
             VD(NDK)=0.            
             DMF(NDK)=0.
             DER(NDK)=0.
             DDR(NDK)=0.
             THTAD(NDK)=0.
             TZ(NDK)=0.
             QD(NDK)=0.
          enddo
          DMFMIN=0.
          goto 165
       endif DOWNDRAFT_OFF




       ! Adjust downdraft mass flux so that evaporation
       ! rate in downdraft is consistent with
       ! precipitation efficiency relationship

       ! PPTFLX = precipitation flux =
       !          updraft supply rate times the
       !          precipitation efficiency
       ! RCED   = Rate of condensate evaporation
       !          in downdraft =
       !          total precipitation rate -
       !          precipitation flux
       !          (this is the precip that evaporates
       !          in downdraft).
       ! PPR    = total fallout of liquid water and
       !          and ice in the updraft.
       ! DMFLFS = adjusted downward mass flux at
       !          the LFS
       ! CNDTNF = condensate in updraft air at the LFS.
       ! DPPTDF = total fallout precip rate in updraft.

       DEVDMF=TDER/DMF(LFS)
       PPR=0.

       PPTFLX=PEFF*USR
       RCED=TRPPT-PPTFLX

       do NM=KLCL,LFS
          PPR=PPR+PPTLIQ(NM)+PPTICE(NM)
       enddo

       if(LFS.ge.KLCL)then

          DPPTDF=(1.-PEFF)*PPR*(1.-EQFRC(LFS))/UMF(LFS)

       else
          DPPTDF=0.
       endif
       CNDTNF=(RLIQ(LFS)+RICE(LFS))*(1.-EQFRC(LFS))

       ! SEE REVISION 037
       DENOM=DEVDMF+DPPTDF+CNDTNF
       if (abs(DENOM) .ge. EPS) then
          DMFLFS=RCED/DENOM
       else
          DMFLFS=1.
       endif

       !  DMFLFS=RCED/(DEVDMF+DPPTDF+CNDTNF)

       if(DMFLFS.gt.0.)then
          TDER=0.
          goto 141
       endif


       ! Adjust the downdraft fluxes.

       ! DDINC = downdraft increase
       ! UPDINC = updraft increase

       DDINC=DMFLFS/DMF(LFS)
       if(LFS.ge.KLCL)then
          UPDINC=(UMF(LFS)-(1.-EQFRC(LFS))*DMFLFS)/UMF(LFS)
       else
          UPDINC=1.
       endif
       do NK=LDB,LFS
          DMF(NK)=DMF(NK)*DDINC
          DER(NK)=DER(NK)*DDINC
          DDR(NK)=DDR(NK)*DDINC
       enddo
       CPR=TRPPT+PPR*(UPDINC-1.)
       TDER=TDER*DDINC


       !  Adjust upward mass flux, mass detrainment
       !  rate, and liquid water detrainment rates to
       !  be consistent with the transfer of the
       !  estimate from the updraft to the downdraft
       !  at the LFS.

       do NK=LC,LFS
          UMF(NK)=UMF(NK)*UPDINC
          UDR(NK)=UDR(NK)*UPDINC
          UER(NK)=UER(NK)*UPDINC
          PPTLIQ(NK)=PPTLIQ(NK)*UPDINC
          PPTICE(NK)=PPTICE(NK)*UPDINC
          DETLQ(NK)=DETLQ(NK)*UPDINC
          DETIC(NK)=DETIC(NK)*UPDINC
       enddo



       !   Set values below the downdraft buoyancy level...

       if(LDB.gt.1)then
          do NK=1,LDB-1
             DMF(NK)=0.
             DER(NK)=0.
             DDR(NK)=0.
             TZ(NK)=0.
             QD(NK)=0.
             THTAD(NK)=0.
             UD(NK)=0.
             VD(NK)=0.
          enddo
       endif

       !   and above the LFS...

       do NK=LFS+1,KX
          DMF(NK)=0.
          DER(NK)=0.
          DDR(NK)=0.
          TZ(NK)=0.
          QD(NK)=0.
          THTAD(NK)=0.
          UD(NK)=0.
          VD(NK)=0.
       enddo
       do NK=LDT+1,LFS-1
          TZ(NK)=0.
          QD(NK)=0.
       enddo



       ! Set limits on the updraft and downdraft mass
       ! fluxes so that the influx into convective drafts
       ! from a given layer is no more than is available
       ! in that layer initially.  Also, do not allow
       ! updraft detrainment to exceed twice the mass
       ! in a layer initially.  Or downdraft detrainment
       ! exceed one time the initial mass in a layer.

       ! LMAX is either the LCL or the LFS.

165    AINCMX=1000.
       LMAX=MAX0(KLCL,LFS)

       ! EMS = total mass of one model grid box
       dtime = timec
       if (deep_timeent_sec > 0) then
          dtime = deep_timeent_sec
       else if (deep_timeent == 'TIMECONV') then
          dtime = timec
       else if (deep_timeent == 'TIMEREFRESH') then
          dtime = timer
       endif
       if (deep_timeent /= 'NIL') then
          do nk=lc,lmax
             if ((uer(nk)-der(nk)) > 0.) aincmx = min(aincmx,ems(nk)/((uer(nk)-der(nk))*dtime))
          enddo
       endif


       !  If the entrainment rates for the updraft
       !  are reasonable physically speaking, then
       !  we should have AINCMX > 1, thus AINC=1 most
       !  of the time here.  If, on the other hand,
       !  the entrainment rates are so large so that
       !  AINCMX < 1, then AINC is equal to this fraction.
       !  (this correction is done to avoid removing
       !  more air than available).

       AINC=1.
       if(AINCMX.lt.AINC)AINC=AINCMX


       !  Save the relevant variables for a unit
       !  updraft and downdraft.  They are adjusted by
       !  the factor AINC to satisfy the stabilization
       !  closure

       NCOUNT=0
       PPTMLT=0.
       TDER2=TDER
       PPTFL2=PPTFLX

       do NK=1,LTOP
          DETLQ2(NK)=DETLQ(NK)
          DETIC2(NK)=DETIC(NK)
          UDR2(NK)=UDR(NK)
          UER2(NK)=UER(NK)
          DDR2(NK)=DDR(NK)
          DER2(NK)=DER(NK)
          UMF2(NK)=UMF(NK)
          DMF2(NK)=DMF(NK)

          ! Between the melting and cloud top levels.

          ! DPLIN=0.5*DPP(I,NK)/(PP0(I,NK)-PP0(I,NK+1))
          ! DLP=ALOG((PP0(I,NK)-0.5*DPP(I,NK))/PP0(I,NK))/ALOG(PP0(I,NK+1)/PP0(I,NK))
          ! THMID(NK+1)=THTA0(NK)+(THTA0(NK+1)-THTA0(NK))*DPLIN
          ! QMID(NK+1)=Q00(I,NK)+(Q00(I,NK+1)-Q00(I,NK))*DLP
       enddo

       ! do NK=ML+1,LTOP-1
       !    PPTMLT=PPTMLT+PPTICE(NK+1)
       ! enddo
       ! PPTML2=PPTMLT
       ! THMID(1)=0.
       ! QMID(1)=0.
       ! DMFMIN=UER(LFS)-EMS(LFS)/TIMEC

       FABE=1.
       STAB=0.95
       if(AINC/AINCMX.gt.0.999)goto 255


       ! BEGINNING OF THE STABILIZATION ITERATION

175    NCOUNT=NCOUNT+1


       ! Evaluate the vertical velocity OMGA (Pa/s) from
       ! the entrainment and detrainment rates.  (vertical
       ! velocity for the environment).

       ! DOMGDP = D(omega) / Dp
       ! DTT    = timestep for the vertical advection
       !          process (later).

       DTT=TIMEC
       do NK = 1, LTOP
          DOMGDP(NK)=-(UER(NK)-DER(NK)-UDR(NK)-DDR(NK))*EMSD(NK)
          if (NK > 1) then
             OMG(NK)=OMG(NK-1)-DPP(I,NK-1)*DOMGDP(NK-1)
             DTT1 = 0.75*DPP(I,NK-1)/(abs(OMG(NK))+1.E-10)
             DTT=AMIN1(DTT,DTT1)
          endif
       enddo
       NSTEP=nint(TIMEC/DTT+1)
       DTIME=TIMEC/FLOAT(NSTEP)
       do NK=1,LTOP
          THPA(NK)=THTA0(NK)
          QPA(NK)=Q00(I,NK)
          QLPA(NK)=QL0(I,NK)
          QIPA(NK)=QI0(I,NK)
          QTPA(NK)=QPA(NK)+QLPA(NK)+QIPA(NK)
       enddo
       ! Revised estimate of compensating vertical motions in the environment (the
       ! resulting profile of OMGA is unrealistic and yields bizarre tendencies)
       OMGA(1:LTOP) = OMG(1:LTOP)

       ! Compute environmental mass flux by mass conservation
       do nk=1,ltop
          emf(nk) = -(umf(nk)+dmf(nk))
       enddo

       !  Determine type of numerical scheme to use for mass flux equations, used
       !  to compute new values for the prognostic variables due to vertical advection
       !  of the environmental properties, as well as detrainment effects from the
       !  updraft and downdraft

       ! Semi-implicit solution for mass flux equations
       do NTC=1,NSTEP


          ! Assign THETA and Q values at the top and bottom of each layer
          ! based on the sign of OMEGA
          do  NK = 1,LTOPM1
             if(OMGA(NK).le.0.)then
                if(NK.eq.1)then
                   NUP(NK)=NK+1
                   THADV(NK)=THPA(NK)
                   QADV(NK)=QPA(NK)
                else
                   NUP(NK)=NK-1
                   THADV(NK)=THPA(NK-1)
                   QADV(NK)=QPA(NK-1)
                endif
             else
                NUP(NK)=NK+1
                THADV(NK)=THPA(NK+1)
                QADV(NK)=QPA(NK+1)
             endif
          enddo
          if(OMGA(LTOP).le.0.)then
             NUP(LTOP)=LTOP-1
             THADV(LTOP)=THPA(LTOP-1)
             QADV(LTOP)=QPA(LTOP-1)
          endif

          ! Compute new values for the potential temperature
          ! including the tendencies due vertical advection
          ! of environmental properties, as well as
          ! detrainment effect from the updraft and downdraft.
          do NK = 2, LTOPM1
             NUPNK = NUP(NK)

             THPA(NK) = ( &
                  (PP0(I,NUPNK)-PP0(I,NK)) * &
                  ( THPA(NK) + UDR(NK)*EMSD(NK)*THTAU(NK)*DTIME &
                  + DDR(NK)*EMSD(NK)*THTAD(NK)*DTIME  ) &
                  - OMGA(NK)*THADV(NK)*DTIME &
                  ) &
                  / &
                  ( &
                  (PP0(I,NUPNK)-PP0(I,NK)) * &
                  (       1. + UDR(NK)*EMSD(NK)*DTIME &
                  + DDR(NK)*EMSD(NK)*DTIME            ) &
                  - OMGA(NK)*DTIME &
                  )

             QPA(NK) = ( &
                  (PP0(I,NUPNK)-PP0(I,NK)) * &
                  (  QPA(NK) + UDR(NK)*EMSD(NK)*QDT(NK)*DTIME &
                  + DDR(NK)*EMSD(NK)*QD(NK) *DTIME  ) &
                  - OMGA(NK)*QADV(NK)*DTIME &
                  ) &
                  / &
                  ( &
                  (PP0(I,NUPNK)-PP0(I,NK)) * &
                  (       1. + UDR(NK)*EMSD(NK)*DTIME &
                  + DDR(NK)*EMSD(NK)*DTIME            ) &
                  - OMGA(NK)*DTIME &
                  )
          enddo

          ! At the surface, no contribution from the updraft.
          NK = 1
          NUPNK = NUP(NK)

          THPA(NK) = ( &
               (PP0(I,NUPNK)-PP0(I,NK)) * &
               ( THPA(NK) + DDR(NK)*EMSD(NK)*THTAD(NK)*DTIME  ) &
               - OMGA(NK)*THADV(NK)*DTIME &
               ) &
               / &
               ( &
               (PP0(I,NUPNK)-PP0(I,NK)) * &
               (       1. + DDR(NK)*EMSD(NK)*DTIME            ) &
               - OMGA(NK)*DTIME &
               )
          QPA(NK) = ( &
               (PP0(I,NUPNK)-PP0(I,NK)) * &
               (  QPA(NK) + DDR(NK)*EMSD(NK)*QD(NK) *DTIME  ) &
               - OMGA(NK)*QADV(NK)*DTIME &
               ) &
               / &
               ( &
               (PP0(I,NUPNK)-PP0(I,NK)) * &
               (       1. + DDR(NK)*EMSD(NK)*DTIME            ) &
               - OMGA(NK)*DTIME &
               )

          ! At the cloud top, only detrainment from the updraft.
          NK = LTOP

          THPA(NK) = ( THPA(NK) + UDR(NK)*THTAU(NK)*EMSD(NK)*DTIME ) &
               / &
               ( 1. + UDR(NK)*EMSD(NK)*DTIME )

          QPA(NK) = (  QPA(NK) + UDR(NK)*QDT(NK)  *EMSD(NK)*DTIME ) &
               / &
               ( 1. + UDR(NK)*EMSD(NK)*DTIME )

       enddo

       ! Specify the "G" grid values after convective adjustment, with
       ! new values after each stabilization iteration.
       do NK=1,LTOP
          THTAG(NK)=THPA(NK)
          QG(NK)=QPA(NK)
          QG(NK)=AMAX1(QG(NK),1.E-9)
       enddo



       !   Specify the "G" grid values after convective adjustment.
       !   New values after each stabilization iteration.

       do NK=1,LTOP
          THTAG(NK)=THPA(NK)
          if (TOTAL_WATER) then
             QLG(NK)=max(QLPA(NK),0.)
             QIG(NK)=max(QIPA(NK),0.)
             QTG(NK)=max(QTPA(NK),0.)
             QG(NK) = max(QTG(NK) - (QLG(NK)+QIG(NK)),0.)
          else
             QLG(NK)=QL0(I,NK)
             QIG(NK)=QI0(I,NK)
             QTG(NK)=Q00(I,NK)+QLG(NK)+QIG(NK)
             QG(NK)=max(QPA(NK),0.)
          endif
       enddo

       ! TOPOMG = vertical motion at the cloud top.

       ! TOPOMG = (UDR(LTOP)-UER(LTOP))*DPP(I,LTOP)*EMSD(LTOP)


       ! Convert THETA to T, freeze all supercooled detrained liquid water
       ! because no supercooled water is supposed to be allowed in the
       ! explicit condensation scheme.

       do NK=1,LTOP
          EXN(NK)=(P00/PP0(I,NK))**(0.2854*(1.-0.28*QG(NK)))
          TG(NK)=THTAG(NK)/EXN(NK)

          ! Temperature effect of melting above the
          ! melting level (ML).

          ! RON is this ok?  ... originally this stuff stayed as liquid, which is weird.  Why don't we let the grid-scale scheme do this?
          if(NK.gt.ML)then
             if (TOTAL_WATER) then
                DTFM(NK) = QLG(NK)*CHLF*EMSD(NK)/CPD
                QIG(NK) = QIG(NK) + QLG(NK)
                QLG(NK) = 0.
             else
                DTFM(NK)=DETLQ(NK)*CHLF*EMSD(NK)/CPD
             endif
          endif
          TG(NK)=TG(NK)+DTFM(NK)
          TVG(NK)=TG(NK)*(1.+0.608*QG(NK))
       enddo


       ! Allow frozen precip to melt over a layer of 200 mb below the melting
       ! level if precipitation is not back to the explicit condensation scheme.
       ! In the current version, this melting effect is not considered
       ! (i.e., DTMLTE = 0). Note here that it should be considered.

       !PV following code does not look complete; does not exist in WRF kf versions; never been used
       !  DTMLTE=0.
       !  TDP=0.
       !  DO K = 1,ML
       !    NK = ML-K+1
       !    TDP=TDP+DPP(I,NK)

       !    Check for a 200 mb layer below the ML.

       !    IF(TDP.LT.2.E4)THEN
       !      DTFM(NK)=DTMLTE
       !      TG(NK)=TG(NK)+DTMLTE*TIMEC
       !    ELSE
       !      DTFM(NK)=DTMLTE*(2.E4+DPP(I,NK)-TDP)/DPP(I,NK)
       !      TG(NK)=TG(NK)+DTMLTE*TIMEC*(2.E4+DPP(I,NK)-TDP)/DPP(I,NK)
       !      GOTO 232
       !    ENDIF
       !   ENDDO
       !232     CONTINUE


       !  11.  COMPUTE NEW CLOUD AND CHANGE IN AVAILABLE BUOYANT ENERGY
       !  =============================================================

       !  The following computations are similar to those for the updraft,
       !  except they are done using the convectively adjusted TG and QG.

       !  Redo the calculations associated with thetrigger function.

       ! (not much comments in the following -
       !  the reader is refered to the sections above).


       THMIXG(I)=0.
       QMIXG(I)=0.
       PMIXG(I)=0.
       do NK = LC,KPBL
          ROCPQ=0.2854*(1.-0.28*QG(NK))
          THMIXG(I)=THMIXG(I)+DPP(I,NK)*TG(NK)*(P00/PP0(I,NK))**ROCPQ
          QMIXG(I)=QMIXG(I)+DPP(I,NK)*QG(NK)
          PMIXG(I)=PMIXG(I)+DPP(I,NK)*PP0(I,NK)
       enddo
       THMIXG(I)=THMIXG(I)/DPTHMXG(I)
       QMIXG(I)=QMIXG(I)/DPTHMXG(I)

       QMIXG(I)=AMAX1( QMIXG(I),1.0E-10 )

       PMIXG(I)=PMIXG(I)/DPTHMXG(I)
       ROCPQ=0.2854*(1.-0.28*QMIXG(I))
       TMIXG(I)=THMIXG(I)*(PMIXG(I)/P00)**ROCPQ
       ES=ALIQ*exp((TMIXG(I)*BLIQ-CLIQ)/(TMIXG(I)-DLIQ))
       ES=min( ES , 0.5*PMIXG(I) )
       QS=0.622*ES/(PMIXG(I)-ES)
       QS=max( min( QS , 0.050 ) , 1.E-6 )
       if(QMIXG(I).gt.QS)then
          RL=XLV0-XLV1*TMIXG(I)
          CPM=CPD*(1.+0.887*QMIXG(I))
          DQSDT = QS*(CLIQ-BLIQ*DLIQ)/((TMIXG(I)-DLIQ)*(TMIXG(I)-DLIQ))
          DQ = (QMIXG(I)-QS)/(1.+RL*DQSDT/CPM)
          TMIXG(I) = TMIXG(I)+RL/CPD*DQ
          QMIXG(I) = QMIXG(I)-DQ
          ROCPQ = 0.2854*(1.-0.28*QMIXG(I))
          THMIXG(I) = TMIXG(I)*(P00/PMIXG(I))**ROCPQ
          TLCLG(I) = TMIXG(I)
          PLCLG(I) = PMIXG(I)
       else
          QMIXG(I)=AMAX1(QMIXG(I),0.)
          EMIX=QMIXG(I)*PMIXG(I)/(0.622+QMIXG(I))
          TLOG=ALOG(EMIX/ALIQ)
          TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
          TLCLG(I)=TDPT-(.212+1.571E-3*(TDPT-TRPL)-4.36E-4*(TMIXG(I)-TRPL))* &
               (TMIXG(I)-TDPT)
          TLCLG(I)=AMIN1(TLCLG(I),TMIXG(I))
          CPORQ=1./ROCPQ
          PLCLG(I)=P00*(TLCLG(I)/THMIXG(I))**CPORQ
       endif
       TVLCL=TLCLG(I)*(1.+0.608*QMIXG(I))
       do NK = LC,KL
          KLCL=NK
          if(PLCLG(I).ge.PP0(I,NK))goto 240
       enddo
240    continue

       KLCLM1=KLCL-1
       DLP=ALOG(PLCLG(I)/PP0(I,KLCLM1))/ALOG(PP0(I,KLCL)/PP0(I,KLCLM1))


       !   Estimate the environmental temperature
       !   and mixing ratio at the LCL

       TENV=TG(KLCLM1)+(TG(KLCL)-TG(KLCLM1))*DLP
       QENV=QG(KLCLM1)+(QG(KLCL)-QG(KLCLM1))*DLP
       TVEN=TENV*(1.+0.608*QENV)
       TVBAR=0.5*(TVG(KLCLM1)+TVEN)


       !  Characteristics of the updraft at the LCL.  (height,
       !  virtual temperature, eq. pot. temperature, pressure,
       !  temperature and detrainment).

       ZLCL=Z0G(I,KLCLM1)+RGASD*TVBAR*ALOG(PP0(I,KLCLM1)/PLCLG(I))/GRAV
       TVAVG=0.5*(TVEN+TG(KLCL)*(1.+0.608*QG(KLCL)))
       PLCLG(I)=PP0(I,KLCL)*exp(GRAV/(RGASD*TVAVG)*(Z0G(I,KLCL)-ZLCL))
       THETEU(KLCLM1)=TMIXG(I)*(1.E5/PMIXG(I))**(0.2854*(1.-0.28*QMIXG(I)))* &
            exp((3374.6525/TLCLG(I)-2.5403)*QMIXG(I)* &
            (1.+0.81*QMIXG(I)))
       ES=ALIQ*exp((TENV*BLIQ-CLIQ)/(TENV-DLIQ))
       ES=min( ES , 0.5*PLCLG(I) )
       QESE=0.622*ES/(PLCLG(I)-ES)
       QESE=max( min( QESE , 0.050 ) , 1.E-6 )
       THTESG(KLCLM1)=TENV*(1.E5/PLCLG(I))**(0.2854*(1.-0.28*QESE))* &
            exp((3374.6525/TENV-2.5403)*QESE*(1.+0.81*QESE))


       ! Compute adjusted ABE (ABEG)

       ABEG=0.
       THTUDL=THETEU(KLCLM1)
       THATA = TMIXG(I)*(1.E5/PMIXG(I))**0.286
       THTFC = THATA*exp((XLV0-XLV1*TLCLG(I))*QMIXG(I)/(CPD*TLCLG(I)))


       ! Loop for ABEG (adjusted).

       do NK=KLCLM1,LTOPM1
          NK1=NK+1
          ES=ALIQ*exp((TG(NK1)*BLIQ-CLIQ)/(TG(NK1)-DLIQ))
          ES=min( ES , 0.5*PP0(I,NK1) )
          QESE=0.622*ES/(PP0(I,NK1)-ES)
          QESE=max( min( QESE , 0.050 ) , 1.E-6 )
          THTESG(NK1)=TG(NK1)*(1.E5/PP0(I,NK1))**(0.2854*(1.-0.28*QESE))* &
               exp((3374.6525/TG(NK1)-2.5403)*QESE*(1.+0.81*QESE))
          if(NK.eq.KLCLM1)then
             DZZ=Z0G(I,KLCL)-ZLCL
          else
             DZZ=DZP(I,NK)
          endif
          BE=((2.*THTUDL)/(THTESG(NK1)+THTESG(NK))-1.)*DZZ
          if(BE.gt.0.)ABEG=ABEG+BE*GRAV
       enddo


       ! Assume greater than 90% of CAPE is
       ! removed by convection during the period TIMEC

       ! DABE = difference between original ABE and adjusted ABE.
       !        Cannot be smaller than 0.1 ABE).
       ! FABE = fraction of ABE left after the convective adjustment.
       ! STAB = 0.95 = stabilization factor.

       ABEOUT(I) = ABE
       CAPEOUT(I) = CAPE
       TAUC(I) =  TIMEC
       CINOUT(I) = CIN
       DABE=AMAX1(ABE-ABEG,0.1*ABE)
       FABE=ABEG/(ABE+1.E-8)

       if(AINC/AINCMX.gt.0.999 .and. FABE.gt.1.05-STAB)then
          goto 265
       endif


       ! IF THE COLUMN IS STABILIZED, THEN OUT OF THE ITERATION LOOP.

       if(FABE.le.1.05-STAB.and.FABE.ge.0.95-STAB)goto 265
       if(NCOUNT.gt.10)then
          goto 265
       endif


       ! If more than 10% of the original CAPE remains,
       ! increase the convective mass flux by the factor AINC

       AINC=AINC*STAB*ABE/(DABE+1.E-8)
255    AINC=AMIN1(AINCMX,AINC)
       AINC=AMAX1(AINC,0.0)

       ! Adjustment of all the updraft/downdraft characteristics.

       !  PPTMLT=PPTML2*AINC
       TDER=TDER2*AINC
       PPTFLX=PPTFL2*AINC

       do NK=1,LTOP
          UMF(NK)=UMF2(NK)*AINC
          DMF(NK)=DMF2(NK)*AINC
          DETLQ(NK)=DETLQ2(NK)*AINC
          DETIC(NK)=DETIC2(NK)*AINC
          UDR(NK)=UDR2(NK)*AINC
          UER(NK)=UER2(NK)*AINC
          DER(NK)=DER2(NK)*AINC
          DDR(NK)=DDR2(NK)*AINC
          DTFM(NK)=0.
       enddo

       !   If the downdraft overdraws the initial mass
       !   available at the LFS, allow PEFF to change
       !   so that updraft mass flux is not the
       !   limiting factor in the total convective
       !   mass flux.  This is usually necessary only
       !   when the downdraft is very shallow.

       DMFMIN=UER(LFS)-EMS(LFS)/TIMEC

       ! DER is large enough (don't forget that DER < 0).
       if(DER(LFS).lt.DMFMIN .and. DMFMIN.lt.0.)then

          !  RF is the proportionality constant for the
          !  decrease of downdraft fluxes.  (RF < 1).

          RF=DMFMIN/DER(LFS)

          do NK=LDB,LFS
             DER2(NK)=DER2(NK)*RF
             DDR2(NK)=DDR2(NK)*RF
             DMF2(NK)=DMF2(NK)*RF
             DER(NK)=DER2(NK)*AINC
             DDR(NK)=DDR2(NK)*AINC
             DMF(NK)=DMF2(NK)*AINC
          enddo

          TDER2=TDER2*RF
          TDER=TDER2*AINC


          ! Modify the UPDINC factor using the new downdraft fluxes.

          if(LFS.ge.KLCL)then
             UPDIN2=1.-(1.-EQFRC(LFS))*DMF(LFS)*UPDINC/UMF(LFS)
          else
             UPDIN2=1.
          endif


          ! Recalculate the condensate production rate
          ! (CPR) and precipitation efficiency.

          CPR=TRPPT+PPR*(UPDIN2-1.)
          PEFF=1.-(TDER+(1.-EQFRC(LFS))*DMF(LFS)*(RLIQ(LFS)+RICE(LFS)))/ &
               (CPR*AINC)

          ! Adjust precipitation fluxes.

          PPTFL2=PEFF*CPR
          PPTFLX=PPTFL2*AINC
          F1=UPDIN2/(AINC*UPDINC)


          ! Updraft characteristics have to change accordingly.
          !VDIR NODEP(UMF)
          do NK=LC,LFS
             UMF2(NK)=UMF(NK)*F1
             UMF(NK)=UMF2(NK)*AINC
             UDR2(NK)=UDR(NK)*F1
             UDR(NK)=UDR2(NK)*AINC
             UER2(NK)=UER(NK)*F1
             UER(NK)=UER2(NK)*AINC
             DETLQ2(NK)=DETLQ(NK)*F1
             DETLQ(NK)=DETLQ2(NK)*AINC
             DETIC2(NK)=DETIC(NK)*F1
             DETIC(NK)=DETIC2(NK)*AINC
             !  PPTML2=PPTML2-PPTICE(NK)*(1.-UPDIN2/UPDINC)
             PPTLIQ(NK)=PPTLIQ(NK)*F1*AINC
             PPTICE(NK)=PPTICE(NK)*F1*AINC
          enddo

          !  PPTMLT=PPTML2*AINC
          UPDINC=UPDIN2
       endif


       ! REDO THE CALCULATIONS FOR CONVECTIVE ADJUSTMENT.

       goto 175

265    continue



       ! Diagnostic outputs WUMAXOUT, RLIQOUT, RICEOUT, AREAUP

       WUMAXOUT(I) = 0.
       do K=1,KX
          NK = KX-K+1
          RLIQOUT(I,NK)  = RLIQ(K)
          RICEOUT(I,NK)  = RICE(K)
          !PV - change units from Kg s-1 to Kg s-1 m-2 for compatibility with BKF and most publications
          UMFOUT(I,NK) = UMF(K)/DXSQ
          DMFOUT(I,NK) = DMF(K)/DXSQ

          ! Vertical integral of rliqout and riceout (in-cloud values)

          RLIQ_INT(I) = RLIQ_INT(I) + RLIQOUT(I,NK)/GRAV * DPP(I,K)
          RICE_INT(I) = RICE_INT(I) + RICEOUT(I,NK)/GRAV * DPP(I,K)

          WUMAXOUT(I) = max( WU(K), WUMAXOUT(I) )

          if (WU(K).gt.0.1) &
               AREAUP(I,NK) = ( UMF(K)*RGASD*TVU(K) ) / &
               ( PP0(I,K) * WU(K)    )


          ! Convert updraft areas into cloud fractions.
          ! Clip cloud fraction to eliminate negative values.

          CLOUDS (I,NK) = min(1.,max(0.,AREAUP(I,NK)/DXDY(I)))
          if (K < KLCL) CLOUDS (I,NK) = 0.

          ! Normalize condensate mixing ratio

          RLIQOUT(I,NK) = RLIQOUT(I,NK)*CLOUDS(I,NK)
          RICEOUT(I,NK) = RICEOUT(I,NK)*CLOUDS(I,NK)

       end do


       ! Activating Convective Momentum Transport or not
       ! First initialize the momentum variables
       ! for the vertical advection calculations
       ! note: same as loop 495 but for winds

       CALC_CMT: if (cmt_type_i /= CMT_NONE) then

          ! Should recompute updraft and downdraft mass fluxes to ensure consistency with entr/detr rates but it crashes!
          ! incoherence entre umf calcule plus haut et le calcul ci-bas pour 2 raisons:
          ! i)Voir SG: prob quand klcl < kpbl; i.e. base du nuage a l'int de la couche melangee
          ! ii) emprunt de umf pour dmf a LFS semble aussi creer une incoherence

          UU(LC)=U00(I,LC)
          VU(LC)=V00(I,LC)
          ! next lines are true given that umf=udr=0 at k=lc-1 (see eqn below for uu)
          UU(LC+1)=U00(I,LC)
          VU(LC+1)=V00(I,LC)

          !  cmt_ecmwf_lambda = 2.*(1. + tanh((nk-let+10)/2.)) !test non-constant lambda

          IF_CMT_ECMWF: if (cmt_type_i == CMT_ECMWF_PH2) then
             !# cmt implemented in phase 2 (ecmwf lambda with bug i.e. udr-uer not updated)
             do NK = LC+1,LTOPM1
                if (nk .gt. let) then  ! apply lambda approach to anvil only (OPS-2019)
                   nudr=(cmt_ecmwf_lambda+1.)*UDR(NK)            
                   nuer=uer(nk)+cmt_ecmwf_lambda*UDR(NK)         
                else
                   nudr=UDR(NK)
                   nuer=uer(nk)
                endif
                UPOLD=UMF(NK-1)-nudr
                UPNEW=UPOLD+nuer
                UU(NK+1)=(UU(NK)*UPOLD+U00(I,NK)*nuer)/UPNEW
                VU(NK+1)=(VU(NK)*UPOLD+V00(I,NK)*nuer)/UPNEW
             enddo
          elseif (cmt_type_i == CMT_ECMWF) then
             !# ecmwf lambda=2 approach without bug and over whole cloud
             do NK = LC+1,LTOPM1
                nudr=(cmt_ecmwf_lambda+1.)*UDR(NK)            
                nuer=uer(nk)+cmt_ecmwf_lambda*UDR(NK)         
                udr(nk)=nudr              
                uer(nk)=nuer              
                UPOLD=UMF(NK-1)-nudr
                UPNEW=UPOLD+nuer
                UU(NK+1)=(UU(NK)*UPOLD+U00(I,NK)*nuer)/UPNEW
                VU(NK+1)=(VU(NK)*UPOLD+V00(I,NK)*nuer)/UPNEW
             enddo
          else 
             !# for all other cmt_type(lambda=0)
             do NK = LC+1,LTOPM1
                UPOLD=UMF(NK-1)-udr(nk)
                UPNEW=UPOLD+uer(nk)
                UU(NK+1)=(UU(NK)*UPOLD+U00(I,NK)*uer(nk))/UPNEW
                VU(NK+1)=(VU(NK)*UPOLD+V00(I,NK)*uer(nk))/UPNEW
             enddo
          endif IF_CMT_ECMWF

          if (DOWNDRAFT) then
             UD(LFS)=EQFRC(LFS)*U00(I,LFS)+(1.-EQFRC(LFS))*UU(LFS+1)
             VD(LFS)=EQFRC(LFS)*V00(I,LFS)+(1.-EQFRC(LFS))*VU(LFS+1)
             do  NK=1,LFS-LDB
                NJ=LFS-NK
                if(NJ.gt.LDT)then
                   UD(NJ)=(UD(NJ+1)*DMF(NJ+1)+U00(I,NJ)*DER(NJ))/DMF(NJ)
                   VD(NJ)=(VD(NJ+1)*DMF(NJ+1)+V00(I,NJ)*DER(NJ))/DMF(NJ)
                else
                   UD(NJ)=UD(LDT+1)
                   VD(NJ)=VD(LDT+1)
                endif
             enddo
          endif


          do NK=1,LTOP+1    !monte de 1 niv pour avoir uadv(ltop)
             UPA(NK)=U00(I,NK)
             VPA(NK)=V00(I,NK)
          enddo

          OMGAEFF(:)=OMGA(:)
          DOMGDP(1)=- ( ONEMINCU*(UER(1)- UDR(1)) + ONEMINCD*(-DER(1)-DDR(1)) )*EMSD(1)
          do NK = 2, LTOP
             DOMGDP(NK)=- ( ONEMINCU*(UER(NK)- UDR(NK)) + ONEMINCD*(-DER(NK)-DDR(NK)) )*EMSD(NK)
             OMGAEFF(NK)=OMGAEFF(NK-1)-DPP(I,NK-1)*DOMGDP(NK-1)
             DTT1 = 0.75*DPP(I,NK-1)/(abs(OMGAEFF(NK))+1.E-10)
             DTT=AMIN1(DTT,DTT1)
          enddo
          NSTEP=nint(TIMEC/DTT+1)
          DTIME=TIMEC/FLOAT(NSTEP)

          ! Determine the type of numerical scheme to use for the momentum mass flux
          ! equations, used to compute new values for the prognostic wind variables
          ! due to vertical advection of the environmental properties, as well as
          ! detrainment effects from the updraft and downdraft.

          ! Semi-implicit solution for the momentum mass flux equations

          ! calculations that don't need to be repeated every ntc steps 
          do nk = 1,ltopm1
             if (omgaeff(nk) > 0. .or. nk == 1) then
                delpup(nk)=pp0(i,nk+1)-pp0(i,nk)
             else
                delpup(nk)=pp0(i,nk-1)-pp0(i,nk)
             endif
             ct3(nk)= onemincu*udr(nk)*emsd(nk)*dtime
             ct4(nk)= onemincd*ddr(nk)*emsd(nk)*dtime
             cmtdenom(nk)= ( delpup(nk) * &
                  (1. + ct3(nk) + ct4(nk) ) &
                  - omgaeff(nk)*dtime     )
             ct1(nk)=   delpup(nk) / cmtdenom(nk)
             ct2(nk)= - omgaeff(nk)*dtime/cmtdenom(nk)
          enddo
          nk = ltop
          ct1(nk)=1./( 1. + onemincu*udr(nk)*emsd(nk)*dtime )
          ct3(nk)=          onemincu*udr(nk)*emsd(nk)*dtime
          ct2(nk)=0.
          ct4(nk)=0.

          SEMI_MOM_STEP: do ntc=1,nstep

             ! Assign upwind wind values based on the sign of OMEGA
             do nk = 1,ltop+1    !monte de 1 niv pour avoir uadv(ltop)
                if(omgaeff(nk).le.0.)then
                   if(nk.eq.1)then
                      uadv(nk)=upa(nk)
                      vadv(nk)=vpa(nk)
                   else
                      uadv(nk)=upa(nk-1)
                      vadv(nk)=vpa(nk-1)
                   endif
                else
                   uadv(nk)=upa(nk+1)
                   vadv(nk)=vpa(nk+1)
                endif
             enddo

             ! Compute new values for the winds
             ! including the tendencies due vertical advection
             ! of environmental properties, as well as
             ! detrainment effect from the updraft and downdraft.

             do nk = 1, ltopm1

                UPA(NK) = & 
                     ( delpup(nk) * &
                     ( UPA(NK)  &
                     + CT3(NK)*UU(NK)  &
                     + CT4(NK)*UD(NK)) &
                     - OMGAEFF(NK)*UADV(NK)*DTIME &
                     ) /cmtdenom(nk) 

                VPA(NK) = &
                     ( delpup(nk) * &
                     ( VPA(NK)  &
                     + CT3(NK)*VU(NK)  &
                     + CT4(NK)*VD(NK)) &
                     - OMGAEFF(NK)*VADV(NK)*DTIME &
                     ) / cmtdenom(nk)
             enddo

             ! At the cloud top, only detrainment from the updraft.
             ! je ne peux pas inclure ca dans la boucle principale sur nk parce que omga(ltop) .ne.0.0
             NK = LTOP
             UPA(NK) = ( UPA(NK) + CT3(NK)*UU(NK) ) / ( 1. + CT3(NK) )
             VPA(NK) = ( VPA(NK) + CT3(NK)*VU(NK) ) / ( 1. + CT3(NK) )

          enddo SEMI_MOM_STEP

          do NK=1,LTOP
             UG(NK)=UPA(NK)
             VG(NK)=VPA(NK)
          enddo

          do K= 1, KX
             NK = KX-K+1
             DUDT(I,NK) = ( UG(K)-U00(I,K) ) / TIMEC
             DVDT(I,NK) = ( VG(K)-V00(I,K) ) / TIMEC
          enddo

          ! Use AZ approach to calculate components(does not work with implicit solver)
          IF_CMT_DIAG: if (cmt_comp_diag) then
             call azcmtcomp(nstep, dtime, ltop, u00(i,:), uu, ud, OMGAEFF, &
                  ct1, ct2, ct3, ct4, idudt1(i,:), idudt2(i,:), idudt3(i,:))
             call azcmtcomp(nstep, dtime, ltop, v00(i,:), vu, vd, OMGAEFF, &
                  ct1, ct2, ct3, ct4, idvdt1(i,:), idvdt2(i,:), idvdt3(i,:))

             ! inverser la coord verticale
             do K = 1, KX
                NK = KX-K+1
                DUDT1(I,NK) = IDUDT1(I,K)
                DVDT1(I,NK) = IDVDT1(I,K)
                DUDT2(I,NK) = IDUDT2(I,K)
                DVDT2(I,NK) = IDVDT2(I,K)
                DUDT3(I,NK) = IDUDT3(I,K)
                DVDT3(I,NK) = IDVDT3(I,K)
             enddo

             do K = 1, KX
                SUMDUDT(I,K) = DUDT1(I,K) + DUDT2(I,K) + DUDT3(I,K)
                SUMDVDT(I,K) = DVDT1(I,K) + DVDT2(I,K) + DVDT3(I,K)
             enddo
          endif IF_CMT_DIAG

          ! NOTE: conventions pour les mass fluxes
          ! - unites de (umf,dmf,udr,uer,ddr,der) = kg/s
          ! - umf est positif ; dmf est negatif (exceptions aux cond frontieres)
          ! - udr,uer,ddr sont positifs
          ! - der est negatif !?!

          ! Impose CMT conservation
          ! Calculate CMT budget from cloud top to sfc (attention a l'ordre des indices)
          INTDUDT = 0.
          INTDVDT = 0.
          INTDP   = 0.
          do K=KX-LTOP+1,KX
             NK = KX-K+1
             INTDUDT  = INTDUDT  + DUDT(I,K)  / GRAV * DPP(I,NK)
             INTDVDT  = INTDVDT  + DVDT(I,K)  / GRAV * DPP(I,NK)
             INTDP = INTDP + DPP(I,NK)
          end do
          ! Apply correction from cloud top to sfc
          do K=KX-LTOP+1,KX
             DUDT(I,K) = DUDT(I,K) - INTDUDT*GRAV/INTDP
             DVDT(I,K) = DVDT(I,K) - INTDVDT*GRAV/INTDP
          end do

       else !# CALC_CMT

          ! If CMT is inactive, initialize updated wind profile as original
          do K=1,KX
             UG(K) = U00(I,K)
             VG(K) = V00(I,K)
          enddo

       endif CALC_CMT

       ! Update the convective counter FLAGCONV

       FLAGCONV(I) = NIC - 1


       ! Note that DTFM is the temperature change due
       ! to freezing of liquid condensate detrained
       ! above the freezing level and melting of
       ! precipitation in the environment below
       ! the melting level.

       do K = 1, KX
          NK = KX-K+1
          DTDT(I,NK) = (TG(K)-TT0(I,K))/TIMEC
          DQDT(I,NK) = (QG(K)-Q00(I,K))/TIMEC
          DQCDT(I,NK) = (QLG(K)-QL0(I,K))/TIMEC
          DQIDT(I,NK) = (QIG(K)-QI0(I,K))/TIMEC
          if(IFEXFB.eq.0)then
             DQRDT(I,NK)=0.
          else
             if(K.lt.KX)then
                DQRDT(I,NK)=(PPTLIQ(K)+PPTICE(K))*PEFF* &
                     AINC*EMSD(K)
             else
                DQRDT(I,NK)=0.
             endif
          endif
          if (.not.TOTAL_WATER) then
             DQCDT(I,NK)=DETLQ(K)*EMSD(K)
             DQIDT(I,NK)=DETIC(K)*EMSD(K)
             DQCDT(I,NK)=max( DQCDT(I,NK) , 0. )
             DQIDT(I,NK)=max( DQIDT(I,NK) , 0. )
          endif
          if (associated(wklclp)) WKLCLP(I,NK) = WKLCL0(I,K)
       enddo

       ! Produce outputs of precipitation fluxes for assimilation purposes

       SUMFLX = 0.
       RNFLX(I,1)  = 0.
       SNOFLX(I,1) = 0.
       do K=KX-1,1,-1
          NK = KX-K+1
          LIQFRAC = ( TT0(I,K)-TBFRZ ) / ( TTFRZ - TBFRZ )
          ICEFRAC = 1 - LIQFRAC
          SUMFLX = SUMFLX + (PPTLIQ(K)+PPTICE(K)) / DXSQ
          RNFLX(I,NK)  = max( min( 1. , LIQFRAC ), 0. ) * SUMFLX
          SNOFLX(I,NK) = max( min( 1. , ICEFRAC ), 0. ) * SUMFLX
       end do

       !   Impose conservation upon request
       ZCRR(I) = PPTFLX/DXSQ
       PRECIP_ADJUSTMENT: if (deep_conserve == 'PRECIP') then
          ! Adjust precipitation to ensure that total water budget is closed

          DQDTDK = 0.
          DQCDTDK = 0.
          do K=1,KX
             NK = KX-K+1
             DQDTDK  = DQDTDK  + DQDT(I,K)  / GRAV * DPP(I,NK)
             DQCDTDK = DQCDTDK + (DQCDT(I,K)+DQIDT(I,K)) / GRAV * DPP(I,NK)
          end do
          ZCRR(I) = max(0.0, -(DQDTDK + DQCDTDK))

          ! Normalize precipitation fluxes with new
          ! (i.e., CONSERVATIVE) precipitation at the
          ! surface

          do K=1,KX
             NORM        = ( ZCRR(I) &
                  / max(RNFLX(I,KX)+SNOFLX(I,KX),1.E-10) )
             RNFLX(I,K)  = NORM * RNFLX(I,K)
             SNOFLX(I,K) = NORM * SNOFLX(I,K)
          end do

       endif PRECIP_ADJUSTMENT

       ! Cloud object handling
       CO_SETUP: if (deep_cloudobj) then

          ! Compute cloud layer mean winds for cloud object advection.
          intu = 0.
          intv = 0.
          intdp = 0.
          do k=klcl,ltop
             intu = intu + ug(k)*dpp(i,k)
             intv = intv + vg(k)*dpp(i,k)
             intdp = intdp + dpp(i,k)
          enddo
          coadvu(i) = intu / intdp
          coadvv(i) = intv / intdp

          ! Update cloud object properties
          wlcl0 = wlcl*exp(coage(i)/deep_codecay)
          coage(i) = coage(i) + delt
          if (deep_codecay > 0.) then
             cowlcl(i) = wlcl0*exp(-coage(i)/deep_codecay)
          else
             cowlcl(i) = 0.
          endif
          cozlcl(i) = zlcl

       endif CO_SETUP

       ! Column has active convection
       ACTIV(I) = .true.

325    continue

       ! Reset cloud object properties if the column is no longer active
       if (.not.ACTIV(I)) then
          coage(i) = 0.
          cowlcl(i) = 0.
          cozlcl(i) = 0.
       endif

    enddo MAIN_I_LOOP

    !  POST-PROCESSING
    !  - - - - - - - -
    !  Convert updraft areas into cloud fractions.
    !  Must be done every timestep because CCFCP
    !  is an automatic array in caller.

    do K=1,KX
       do I=1,IX
          CLOUDS (I,K) = min(1.,max(0.,AREAUP(I,K)/DXDY(I)))
       end do
    end do

    !  Convert tendencies from mixing ratio back to specific humidity
    if (MIXING_RATIO) then
       do k=1,kx
          nk = kx-k+1
          dqdt(:,k) = dqdt(:,k) / (1.+Q00(:,nk))**2
          dqcdt(:,k) = dqcdt(:,k) / (1.+Q00(:,nk))**2
          dqidt(:,k) = dqidt(:,k) / (1.+Q00(:,nk))**2      
       enddo
    endif

!!$       do i=1,ix
!!$          clouds(i,1:kx) = areaup(i,1:kx)/dxdy(i)
!!$       end do


    !  COUNTER FOR CONVECTION
    !  ----------------------

    !  COUNTER IS SET TO ZERO AT FIRST TIMESTEP FOR CONSISTENCY
    !  SINCE PRECIPITATION RATE IS SET TO ZERO IN caller AND
    !  TEMPERATURE/HUMIDITY TENDENCIES ARE NOT APPLIED IN DYNAMICS.
    !  HOWEVER, CLOUD FRACTION AND LIQUID/SOLID WATER CONTENT
    !  ARE KEPT FOR RADIATION CALCULATIONS AT KOUNT.EQ.1.

    if (KOUNT.eq.0) then
       do I=1,IX
          FLAGCONV(I) = 0.
       end do
    endif

    do I=1,IX
       if (FLAGCONV(I).gt.0.) KKFC(I)=1.
    end do

!!$       if (kount.eq.0) then
!!$          flagconv(1:ix) = 0.
!!$       else
!!$          where(flagconv(1:ix) > 0.) kkfc(1:ix) = 1.
!!$       endif

    return
  end subroutine kfrpn4

end module kfrpn
