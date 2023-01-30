!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
! version 3; Last Modified: May 7, 2008.
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

module kfmid
   implicit none
   private
   public :: kfmid1

contains

   subroutine kfmid1(tp1, qp1, ub, vb, wb, gzm, sigt, ps, dxdy, latr, &
        deep_cape, deep_active, cond_prflux, &
        dtdt, dqdt, dudt, dvdt, dqcdt, dqidt, &
        active, zcrr, clouds, rliqout, riceout, rnflx, snoflx, mcd, mpeff, mainc, &
        mrk2, delt, ix, kx)
      use, intrinsic :: iso_fortran_env, only: INT64
      use cnv_options, only: mid_peff, mid_maxcape, mid_dpdd, mid_conserve, &
           mid_minbase, mid_depth, mid_minemf, mid_peff_const, mid_emffrac, &
           mid_emfmod
      use tdpack, only: fohr, RGASD, GRAV, TRPL, PI, CHLF, CPD
      use integrals, only: int_profile, INT_OK
      use debug_mod, only: init2nan
      use tpdd, only: tpdd1
      use ens_perturb, only: ens_spp_get
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      ! Argument declaration
      integer, intent(in) :: ix                                    !Row length
      integer, intent(in) :: kx                                    !Number of prognostic levels
      real, intent(in) :: delt                                     !Time step (s)
      real, dimension(:,:), intent(in) :: tp1                      !Dry air temperature (K)
      real, dimension(:,:), intent(in) :: qp1                      !Specific humidity (kg/kg)
      real, dimension(:,:), intent(in) :: ub                       !U-component wind speed (m/s)
      real, dimension(:,:), intent(in) :: vb                       !V-component wind speed (m/s)
      real, dimension(:,:), intent(in) :: wb                       !Vertical motion (m/s)
      real, dimension(:,:), intent(in) :: gzm                      !Geopotential heights (m)
      real, dimension(:,:), intent(in) :: sigt                     !Eta level value (0-1)
      real, dimension(:), intent(in) :: ps                         !Surface pressure (Pa)
      real, dimension(:), intent(in) :: dxdy                       !Grid cell area (m2)
      real, dimension(:), intent(in) :: latr                       !Latitude (radians)
      real, dimension(:), intent(in) :: deep_cape                  !CAPE from deep convective scheme (J/kg/m2)
      real, dimension(:), intent(in) :: deep_active                !Deep convective activity mask (0 or 1)
      real, dimension(:,:), intent(in) :: cond_prflux              !Precipitation flux from gridscale scheme
      real, dimension(:,:), pointer :: mrk2                        !Markov chains for stochastic perturbations
      real, dimension(:,:), intent(out) :: dtdt                    !Temperature tendency (K/s)
      real, dimension(:,:), intent(out) :: dqdt                    !Humidity tendency (kg/kg/s)
      real, dimension(:,:), intent(out) :: dudt                    !U-component wind tendency (m/s2)
      real, dimension(:,:), intent(out) :: dvdt                    !V-component wind tendency (m/s2)
      real, dimension(:,:), intent(out) :: dqcdt                   !Grid-scale liquid condensate tendency (kg/kg/s)
      real, dimension(:,:), intent(out) :: dqidt                   !Grid-scale ice condensate tendency (kg/kg/s)
      real, dimension(:), intent(out) :: active                    !Scheme activity mask (0 or 1)
      real, dimension(:), intent(out) :: zcrr                      !Rainfall rate (flux units: kg/m2/s)
      real, dimension(:,:), intent(out) :: clouds                  !Cloud fraction (0-1)
      real, dimension(:,:), intent(out) :: rliqout                 !Updraft liquid water mixing ratio (kg/kg)
      real, dimension(:,:), intent(out) :: riceout                 !Updraft ice water mixing ratio (kg/kg)
      real, dimension(:,:), intent(out) :: rnflx                   !Rain water flux (kg/m2/s)
      real, dimension(:,:), intent(out) :: snoflx                  !Snow flux (kg/m2/s)
      real, dimension(:), intent(out) :: mcd                       !Cloud depth (m)
      real, dimension(:), intent(out) :: mpeff                     !Precipitation efficiency (fractional)
      real, dimension(:), intent(out) :: mainc                     !Multiplicative closure adjustment

      !@Author Jack Kain and JM Fritsch (Oct 14,1990)

      !@Object
      !          compute the effects of mid-level convection using a modified version of
      !          the Kain-Fritsch convective paramerization scheme.

      !Local variable declarations
      integer :: i,k,iflag,klcl,kmin,kmix,kfrz,kpbl,lc,lcl,ldt,  &
           ldb,let,lmax,ltop,ltop1,ltopm1,lvf,lfs,nd,nd1,ndk,    &
           nj,nk,nk1,nm,nstep,ntc,nlayrs,ktop,klm,kl,mxlayr,     &
           timec,klclm1
      integer, dimension(1) :: cond_maxlev
      integer, dimension(kx) :: nup
      real :: a1,abe,ainc,aincm1,aincm2,aincmx,au0,be,boterm,      &
           cldhgt,clvf,cndtnf,cporq,cpr,dlp,dpt,dtt,dtt1,dzz,&
           devdmf,dmfmin,dmflfs,dptt,zlcl,dpdd,dpddmx,dpptdf,ddinc,     &
           dqsdt,dtmp,dtime,dumfdp,es,ee,ee1,ee2,emix,effq,       &
           enterm,f1,f2,frc,frc1,liqfrac,icefrac,p165,peff,&
           ppr,pptflx,qenv,qese,qnewic,qnewlq,&
           qnwfrz,qs,qsrh,r1,rl,rf,rei,rdd,rocpq,rced,   &
           rtmp,sumflx,tvlcl,tlog,tenv,tven,tdpt,&
           tvavg,tvbar,thtudl,thtmin,thttmp,ttmp,ttemp,     &
           tmpliq,tmpice,trppt,tsat,tder,t1rh,     &
           tu10,tu95,tudl,ud1,ud2,usr,udlbe,updinc,updin2,upold,&
           upnew,vmflcl,norm,wlcl,wtw,dpup,     &
           dpdown,dqdtdk,dqcdtdk,denom,tvdiff,intdudt,      &
           intdvdt,intdp,dxsq,rholcl,wabs,zmix,thta,delp,    &
           nuer,nudr,lambda,rhmixg,rad,ww
      real, dimension(1) :: subcloud_pw, subcloud_pws
      real, dimension(ix) :: dpthmxg,pmixg,tmixg,qmixg,thmixg,tlclg,    &
           plclg,psb,thvmixg,minemf,radmult,dpddmult
      real, dimension(kx) :: ddr,der,detic,detlq,dmf,domgdp,dtfm,ems,   &
           emsd,eqfrc,exn,omga,pptice,pptliq,qadv,qd,qdt,qg,qicout,     &
           qlqout,qpa,qu,ratio2,rice,rliq,tg,theted,thetee,theteu,thadv,&
           thpa,thtad,thtag,thtau,thta0,thtes,tu,tvd,tvqu,   &
           tvu,tz,udr,uer,umf,wu,uu,vu,ud,vd,ug,vg,upa,vpa,qlpa,        &
           qipa,qlg,qig,qtpa,qtdt,qtg,thpai,qtpai,qlpai,qipai,ta,tb,tc, &
           td_thta,td_qt,td_ql,td_qi,tri_thta,tri_qt,tri_ql,tri_qi,     &
           workk,emf,td_q,tri_q,qpai,td_u,td_v,tri_u,tri_v,upai,vpai
      real, dimension(kx+1)  :: omg
      real, dimension(1,kx) :: pp0c,q00c,qst1c
      real, dimension(ix,kx) :: tt0,tv00,q00,u00,v00,wz0,dzp,dpp,qst1,  &
           pp0,z0g,sigkfc,ql0,qi0,thv0,areaup
      logical :: need_buoyancy_sorting,need_above_let,need_mixing
!!$      character(len=8) :: cdmf
!!$      character(len=2) :: cderl,cderr,cddrl,cddrr
!!$      character(len=512) :: clfs

      ! External subprograms
      external tpmix
      external condload
      external envirtht
      external prof5

      ! Basic parameters
#include "clefcon.cdk"
      include "phyinput.inc"

      ! User-adjustable parameters
      real, parameter :: RADIUS = 1500.                         !Cloud radius for cloud model (m)
      real, parameter :: WLCL_MULT = 1.                         !Multiplier of wz_env to get cloud base vertical motion r(realted to width of vertical motion spectrum)
      real, parameter :: RHMIN = -1.                            !Relative humidity threshold for triggering (fraction)
      real, parameter :: DPMIX = 10e2                           !Depth of the mixing layer (Pa)
      logical, parameter :: USE_DOWNDRAFT = .true.              !Include downdraft calculations

      ! Internal parameters
      real, parameter :: EPS = 1.E-8
      real, parameter :: P00 = 1.E5                             !Reference pressure (Pa)
      REAL, PARAMETER :: CV = 717.
      real, parameter :: RHIC = 1.
      real, parameter :: RHBC = 0.90
      real, parameter :: TTFRZ = 268.16
      real, parameter :: TBFRZ = 248.16
      real, parameter :: C5 = 1.0723E-3
      real, parameter :: RATE = 0.005                           !Autoconversion rate
      real, parameter :: MAXZPAR = 10000.                       !Maximum departure height for parcel tests (m)
      real, parameter :: BETAENT = 1.05
      real, parameter :: PEFFMAX = 0.9                          !Maximum precipitation efficiency (fractional)
      real, parameter :: ROVG = RGASD/GRAV
      real, parameter :: GDRY = -GRAV/CPD
      real, parameter :: ONEOVG = 101.9368
      real, parameter :: BETAW = 6822.459384
      real, parameter :: GAMW  = 5.13948
      real, parameter :: BETAI = 6295.421
      real, parameter :: GAMI = 0.56313
      real, parameter :: ONEMINC = 1.0 !PV Hor pressure gradient effects on CMT; value of C=0.0 suggested in Romps(2012)
      real, parameter :: XLV0 = 3.147E+6                        !Constant for latent heat of vapourization
      real, parameter :: XLV1 = 2369.                           !Constant for latent heat of vapourization
      real, parameter :: XLS0 = 2.905E+6                        !Constant for latent heat of vapourization
      real, parameter :: XLS1 = 259.532                         !Constant for latent heat of vapourization
      real, parameter :: ALIQ = 613.3                           !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: BLIQ = 17.502                          !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: CLIQ = 4780.8                          !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: DLIQ = 32.19                           !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: AICE = 613.2                           !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: BICE = 22.452                          !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: CICE = 6133.0                          !Constant for saturation vapour press (Buck 1981; JAM)
      real, parameter :: DICE = 0.61                            !Constant for saturation vapour press (Buck 1981; JAM)

      call init2nan(dpthmxg,pmixg,tmixg,qmixg,thmixg,tlclg)
      call init2nan(plclg,psb,thvmixg, omg,theteu,thadv)
      call init2nan(ddr,der,detic,detlq,dmf,domgdp,dtfm,ems,qicout)
      call init2nan(emsd,eqfrc,exn,omga,pptice,pptliq,qadv,qd,qdt,qg)
      call init2nan(qlqout,qpa,qu,ratio2,rice,rliq,tg,theted,thetee)
      call init2nan(thpa,thtad,thtag,thtau,thta0,thtes,tu,tvd,tvqu)
      call init2nan(tvu,tz,udr,uer,umf,wu,uu,vu,ud,vd)
      call init2nan(ug,vg,upa,vpa,qlpa,ta,tb,tc,vpai)
      call init2nan(qipa,qlg,qig,qtpa,qtdt,qtg,thpai,qtpai,qlpai,qipai)
      call init2nan(td_thta,td_qt,td_ql,td_qi,tri_thta,tri_qt,tri_ql,tri_qi)
      call init2nan(workk,emf,td_q,tri_q,qpai,td_u,td_v,tri_u,tri_v,upai)
      call init2nan(pp0c,q00c,qst1c)
      call init2nan(tt0,tv00,q00,u00,v00,wz0,dzp,dpp,qst1)
      call init2nan(pp0,z0g,sigkfc,ql0,qi0,thv0,areaup)

      ! Basic initializations
      kl = kx
      klm = kl-1
      active(:) = 0.
      zcrr(:) = 0.
      dtdt(:,:) = 0.
      dqdt(:,:) = 0.
      dqcdt(:,:) = 0.
      dqidt(:,:) = 0.
      dudt(:,:) = 0.
      dvdt(:,:) = 0.
      rliqout(:,:) = 0.
      riceout(:,:) = 0.
      rnflx(:,:) = 0.
      snoflx(:,:) = 0.
      areaup(:,:) = 0.
      clouds(:,:) = 0.
      mcd = 0.
      mpeff = 0.
      mainc = 0.

      ! The time scale for this closure is set by the time step
      timec = delt

      ! Stochastic perturbations of sensitive parameters
      minemf(:) = ens_spp_get('mid_minemf', mrk2, mid_minemf)
      radmult(:) = ens_spp_get('crad_mult', mrk2, 1.) 
      dpddmult(:) = ens_spp_get('dpdd_mult', mrk2, 1.)
      
      ! 4. INPUT A VERTICAL SOUNDING (ENVIRONMENTAL)
      ! ============================================
      do I = 1, IX
         PSB(I)=PS(I)*1.E-3
      enddo
      DO K=1,KX
         DO I=1,IX
            SIGKFC(I,K) = SIGT(I,K)
         END DO
      END DO

      ! Construct the sounding
      do K = 1, KX
         NK = KX-K+1
         do I = 1, IX
            ! In the RPN physics, the K indices goes
            ! from 1 at the top of the model to KX
            ! at the surface.
            ! In the Kain-Fritsch subroutine, we use
            ! the opposite (of course).
            ! The switch is done in this loop.
            PP0(I,K) = 1.E3*(SIGKFC(I,NK)*PSB(I))
            TT0(I,K) = TP1(I,NK)
            U00(I,K) = UB(I,NK)
            V00(I,K) = VB(I,NK)
            Q00(I,K) = max(qp1(i,nk),1.0E-10 )
            QL0(I,K) = 0.
            QI0(I,K) = 0.
            Z0G(I,K) = GZM(I,NK)

            ! If Q0 is above saturation value, reduce it
            ! to saturation level
            ES=ALIQ*exp((BLIQ*TT0(I,K)-CLIQ)/(TT0(I,K)-DLIQ))
            ES=min( ES, 0.5*PP0(I,K) )
            QST1(I,K)=0.622*ES/(PP0(I,K)-ES)
            QST1(I,K) = max( min(QST1(I,K),0.050) , 1.E-6 )
            Q00(I,K)  = min(QST1(I,K),Q00(I,K))
            TV00(I,K) = TT0(I,K) * (1. + 0.608*Q00(I,K) )
            ROCPQ=0.2854*(1.-0.28*Q00(I,K))
            THV0(I,K) = TV00(I,K)*(1.E5/PP0(I,K))**ROCPQ
            WZ0(I,K)  = WB(I,NK)

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

      do K = 2, KX
         do I = 1, IX
            DZP(I,K-1) = Z0G(I,K)-Z0G(I,K-1)
         end do
      end do

      ! At the top: do not put the thickness to zero but instead
      ! put the level kx-1 in kx
      do I=1,IX
         DZP(I,KX) = DZP(I,KX-1)
      end do


      ! MAIN LOOP ACROSS I FOR KFMID
      ! ================
      MAIN_I_LOOP: do I = 1, IX

         ! Column initializations
         RLIQ(:) = 0.
         RICE(:) = 0.
         DXSQ = DXDY(I)

         ! If there is sufficient CAPE for active deep convection, the mid-level
         ! convective scheme does not activate
         if (mid_maxcape > 0 .and. deep_cape(i) > mid_maxcape .and. nint(deep_active(i)) == 1) cycle

         ! Highest TOP
         do K = 2, KX
            if(Z0G(I,K) - Z0G(I,1) .lt. MAXZPAR ) KTOP=K
         end do


         ! 6. TRIGGER FUNCTION
         ! ===================

         ! Parameterize the minimum environmental mass flux on request
         if (mid_emfmod == 'LATITUDE') then
            minemf(i) = min(1./max(abs(sin(latr(i))),epsilon(latr)),10.) * minemf(i)
         elseif (mid_emfmod == 'LATMOD1') then
            ww = max(0., min(1., 0.1*( 40. - 57.2958*abs(latr(i)) ) ))
            minemf(i) = min(1./max(abs(sin(latr(i))),epsilon(latr)),10.) * &
                 ( ww*(1.5*minemf(i)) + (1.-ww)*(1.*minemf(i)) )
         elseif (mid_emfmod == 'LATMOD2') then
            ww = max(0., min(1., 0.1*( 40. - 57.2958*abs(latr(i)) ) ))
            minemf(i) = min(1./max(abs(sin(latr(i))),epsilon(latr)),10.) * &
                 ( ww*(1.25*minemf(i)) + (1.-ww)*(0.8*minemf(i)) )
         endif

         ! Find the first level that meets the vertical motion criterion
         k = 1
         do while (dxdy(i)*wz0(i,k)*(pp0(i,k)/(RGASD*tv00(i,k))) < minemf(i) .and. k < kx)
            k = k+1
         enddo
         if (k < kx) then
            kmix = k
         else
            cycle MAIN_I_LOOP
         endif

         ! Begin looping vertically to find an active layer
         TRIGGER: do while (kmix < kx .and. nint(active(i)) == 0)
            lc = kmix

            ! Climb to at least the minimum parcel departure height
            do while (Z0G(I,LC)-Z0G(I,1) < mid_minbase .and. LC < KX)
               LC = LC+1
            enddo

            ! No levels (or parcels) in the lowest
            ! 300 mb were found to be unstable
            ! ==> go to the next grid point
            if (Z0G(I,LC)-Z0G(I,1) >= MAXZPAR) cycle MAIN_I_LOOP
            MXLAYR=0

            ! Find the number of layers needed to
            ! have a mixed layer of at least 60 mb
            NLAYRS=0             ! number of model layers in the mixed layer
            DPTHMXG(I)=0.        ! pressure depth of the mixed layer
            nk = lc
            do while (dpthmxg(i) <= DPMIX .and. nk <= kx)
               dpthmxg(i) = dpthmxg(i) + dpp(i,nk)
               nlayrs = nlayrs + 1
               nk = nk + 1
            enddo
            if (dpthmxg(i) <= DPMIX) cycle MAIN_I_LOOP
            KPBL=LC+NLAYRS-1
            KMIX=LC+1
            THMIXG(I)=0.         ! potential temperature of the mixed layer
            QMIXG(I)=0.
            ZMIX=0.
            PMIXG(I)=0.
            DPTHMXG(I)=0.

            ! Integrated properties
            do NK = LC,KPBL
               DPTHMXG(I)=DPTHMXG(I)+DPP(I,NK)
               ROCPQ=0.2854*(1.-0.28*Q00(I,NK))
               THMIXG(I)=THMIXG(I)+DPP(I,NK)*TT0(I,NK)*(P00/PP0(I,NK))**ROCPQ
               QMIXG(I)=QMIXG(I)+DPP(I,NK)*Q00(I,NK)
               ZMIX=ZMIX+DPP(I,NK)*Z0G(I,NK)
               PMIXG(I)=PMIXG(I)+DPP(I,NK)*PP0(I,NK)
            enddo

            ! Average to obtain the mixed properties
            THMIXG(I)=THMIXG(I)/DPTHMXG(I)
            QMIXG(I)=QMIXG(I)/DPTHMXG(I)
            QMIXG(I)=AMAX1( QMIXG(I),1.0E-10 )
            ZMIX=ZMIX/DPTHMXG(I)
            PMIXG(I)=PMIXG(I)/DPTHMXG(I)
            ROCPQ=0.2854*(1.-0.28*QMIXG(I))
            TMIXG(I)=THMIXG(I)*(PMIXG(I)/P00)**ROCPQ
            THVMIXG(I)=THMIXG(I)*(1.+0.608*QMIXG(I))
            RHMIXG=FOHR(QMIXG(I),TMIXG(I),PMIXG(I))

            ! If parcel does not meet the relative humidity
            ! threshold, move to next level and continue
            ! trigger search.
            if (rhmixg < RHMIN) cycle TRIGGER

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
            nk = lc
            do while (plclg(i) < pp0(i,nk) .and. nk <= kl)
               if (nk > kpbl) then
                  dzz = z0g(i,nk+1)-z0g(i,nk)
                  udlbe = (2.*thvmixg(i)/(thv0(i,nk+1)+thv0(i,nk))-1.)*dzz
               endif
               nk = nk + 1
            enddo
            klcl = min(nk,kl)
            klclm1 = klcl - 1
            dlp = log(plclg(i)/pp0(i,klclm1)) / log(pp0(i,klcl)/pp0(i,klclm1))

            ! Estimate environmental temperature and mixing ratio at the LCL
            TENV=TT0(I,KLCLM1)+(TT0(I,KLCL)-TT0(I,KLCLM1))*DLP
            QENV=Q00(I,KLCLM1)+(Q00(I,KLCL)-Q00(I,KLCLM1))*DLP
            TVEN= TENV*(1.+0.608*QENV)
            TVBAR=0.5*(TV00(I,KLCLM1)+TVEN)

            ! Height of the LCL (to  avoid getting negative dzz later on
            ! use simple linear interpolation to get zlcl)
            ZLCL= Z0G(I,KLCLM1)+ (Z0G(I,KLCL) - Z0G(I,KLCLM1)) * &
                 ( PLCLG(I)-PP0(I,KLCLM1) )/ ( PP0(I,KLCL) - PP0(I,KLCLM1) )

            ! We do not want zlcl to be at the level z0g(klcl), otherwise we get dzz=0.
            ! Therefore we impose a difference of 1 metre.
            ZLCL= min (ZLCL, Z0G(I,KLCL) - 1.)


            ! 7. CHARACTERISTICS OF PARCELS AT THE LCL
            ! ====================================================
            ! Note that KLCLM1 here is the level just below
            ! the LCL, whereas KLCL is the level just
            ! above it.

            ! THETEU = equivalent potential temperature
            ! of the updraft just below the LCL.
            THETEU(KLCLM1)=TMIXG(I)*(1.E5/PMIXG(I))**(0.2854*(1.-0.28*QMIXG(I)))* &
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
            ! pot. temperature of the environment at the LCL
            QESE=0.622*ES/(PLCLG(I)-ES)
            QESE=max( min( QESE , 0.050 ) , 1.E-6 )
            wlcl = WLCL_MULT * (wz0(i,klclm1)+(wz0(i,klcl)-wz0(i,klclm1))*dlp)
            THTES(KLCLM1)=TENV*(1.E5/PLCLG(I))**(0.2854*(1.-0.28*QESE))* &
                 exp((3374.6525/TENV-2.5403)*QESE*(1.+0.81*QESE))
            WTW=WLCL*WLCL

            ! No upward motion ==> take next slab
            ! and test for convective instability.
            if(WLCL.le.0.) cycle TRIGGER

            ! Diagnostic properties at the LCL
            TVLCL=TLCLG(I)*(1.+0.608*QMIXG(I))
            RHOLCL=PLCLG(I)/(RGASD*TVLCL)
            LCL=KLCL
            LET = LCL

            ! 8. COMPUTE UPDRAFT PROPERTIES
            ! =============================

            ! Initial updraft properties at the LCL
            RAD = RADIUS * RADMULT(I)
            WU(KLCLM1)=WLCL
            AU0=PI*RAD*RAD                      ! area of updraft at cloud base
            UMF(KLCLM1)=RHOLCL*AU0*WU(KLCLM1)   ! updraft mass flux
            VMFLCL=UMF(KLCLM1)                  ! vertical mass flux at the LCL
            UPOLD=VMFLCL
            UPNEW=UPOLD
            RATIO2(KLCLM1)=0.                   ! degree of glaciation

            ! Updraft characteristics
            UER(KLCLM1)=0.                      ! environmental entrainment rate
            ABE=0.                              ! available buoyant energy
            TRPPT=0.                            ! total rate of precip. production
            EQFRC(KLCLM1)=1.                    ! fraction of environmental air in mixed parcels

            ! Updraft properties
            TU(KLCLM1)=TLCLG(I)                 ! updraft temperature
            TVU(KLCLM1)=TVLCL                   ! updraft virtual temperature
            QU(KLCLM1)=QMIXG(I)                 ! updraft mixing ratio
            RLIQ(KLCLM1)=0.                     ! updraft liquid water
            RICE(KLCLM1)=0.                     ! updraft ice

            ! Updraft/environment interactions
            QLQOUT(KLCLM1)=0.                   ! liquid water fallout from the updraft
            QICOUT(KLCLM1)=0.                   ! ice fallout from the updraft
            DETLQ(KLCLM1)=0.                    ! detrained liquid water from the updraft
            DETIC(KLCLM1)=0.                    ! detrained ice from the updraft
            PPTLIQ(KLCLM1)=0.                   ! liquid precipitation from conv. updraft
            PPTICE(KLCLM1)=0.                   ! solid precipitation from conv. updraft

            ! Phase transitions
            IFLAG = 0                           ! freezing flag
            KFRZ=LC                             ! freezing level

            ! The amount of convective available potential
            ! energy (CAPE) is calculated with respect to
            ! an undiluted parcel ascent.
            THTUDL=THETEU(KLCLM1)               ! equiv. pot. temp. of an undiluted parcel
            TUDL=TLCLG(i)                       ! temperature of an undiluted parcel

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
            UPDRAFT: do NK=KLCLM1,KLM
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
               ! (considering the entire glaciation process,
               ! from TTFRZ to TBFRZ)
               ! R1   = fractional conversion to glaciation
               ! (but this time considering only the
               ! remaining glaciation from TU(NK) to
               ! TBFRZ).
               FREEZING: if(TU(NK1).le.TTFRZ.and.IFLAG.lt.1)then

                  ! Fractional freezing processes
                  if(TU(NK1).gt.TBFRZ)then
                     ! Within the glaciation buffer zone
                     if(TTEMP.gt.TTFRZ)TTEMP=TTFRZ
                     FRC1=(TTEMP-TU(NK1))/(TTFRZ-TBFRZ)
                     R1=(TTEMP-TU(NK1))/(TTEMP-TBFRZ)
                  else
                     ! Glaciation processes are over...
                     FRC1=(TTEMP-TBFRZ)/(TTFRZ-TBFRZ)
                     R1=1.
                     IFLAG=1
                  endif

                  ! In this glaciation process, the new liquid
                  ! generated in the updraft directly goes into
                  ! QNEWFRZ.
                  QNWFRZ=QNEWLQ

                  ! Increase new ice and decrease new liquid
                  ! following the R1 fractional glaciation
                  QNEWIC=QNEWIC+QNEWLQ*R1*0.5
                  QNEWLQ=QNEWLQ-QNEWLQ*R1*0.5
                  EFFQ=(TTFRZ-TBFRZ)/(TTEMP-TBFRZ)
                  TTEMP=TU(NK1)

               endif FREEZING

               ! Updraft vertical velocity and precipitation fallout.
               if(NK.eq.KLCLM1)then
                  ! For the level just below the LCL
                  BE=(TVLCL+TVU(NK1))/(TVEN+TV00(I,NK1))-1.
                  BOTERM=2.*(Z0G(I,NK1)-ZLCL)*GRAV*BE/1.5
                  ENTERM=0.
                  DZZ=Z0G(I,NK1)-ZLCL
               else
                  ! Above the LCLC
                  BE=(TVU(NK)+TVU(NK1))/(TV00(I,NK)+TV00(I,NK1))-1.
                  BOTERM=2.*DZP(I,NK)*GRAV*BE/1.5
                  ENTERM=2.*UER(NK)*WTW/UPOLD
                  DZZ=DZP(I,NK)
               endif

               ! Condensation in the updraft.
               call CONDLOAD(RLIQ(NK1),RICE(NK1),WTW,DZZ,BOTERM,ENTERM, &
                    RATE,QNEWLQ,QNEWIC,QLQOUT(NK1),QICOUT(NK1), &
                    GRAV)
               WABS=sqrt(max(abs(WTW),1.E-10))
               WU(NK1)=WTW/WABS

               ! If the vertical velocity is less than zero, exit the updraft loop
               ! and if the cloud is deep enough, finalize the updraft calculations.
               if(WU(NK1).le.0.) exit UPDRAFT

               ! Update ABE for undilute ascent.  (Note that we still did not
               ! consider any entrainment effect except for the calculation of w.)
               THTES(NK1)=TT0(I,NK1)*(1.E5/PP0(I,NK1))**(0.2854*(1.-0.28*QST1(I,NK1)))* &
                    exp((3374.6525/TT0(I,NK1)-2.5403)*QST1(I,NK1)*(1.+0.81*QST1(I,NK1)))
               UDLBE=((2.*THTUDL)/(THTES(NK)+THTES(NK1))-1.)*DZZ
               if(UDLBE.gt.0.) ABE = ABE + UDLBE*GRAV

               ! Other effects of cloud glaciation if within the specified
               ! temperature interval (on temperature, water variables, etc...).
               if(FRC1.gt.1.E-6)then
                  ! Phase transition effects
                  call DTFRZNEW2(TU(NK1),PP0(I,NK1),THETEU(NK1),QU(NK1),RLIQ(NK1), &
                       RICE(NK1),RATIO2(NK1),QNWFRZ,RL,FRC1,EFFQ, &
                       IFLAG,XLV0,XLV1,XLS0,XLS1, &
                       ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
               endif

               ! Call subroutine to calculate ENVIRONMENTAL potential temperature
               ! within glaciation interval.  THETAE must be calculated with
               ! respect to the same degree of glaciation for all entraining air.
               call ENVIRTHT(PP0(I,NK1),TT0(I,NK1),Q00(I,NK1), &
                    THETEE(NK1),RATIO2(NK1),RL, &
                    ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

               ! Trace the level of minimum eq. potential temperature for the downdraft.
               if(NK1.eq.KLCL)THTMIN=THTES(NK1)
               THTMIN=AMIN1(THTES(NK1),THTMIN)
               if(THTMIN.eq.THTES(NK1))KMIN=NK1

               ! BEGINNING OF ENTRAINMENT/DETRAINMENT CALCULATIONS

               ! Compute basic entrainment properties
               REI=VMFLCL*DPP(i,NK1)*0.03/RAD   ! rate of environmental inflow
               TVQU(NK1)=TU(NK1)*(1.+0.608*QU(NK1)-RLIQ(NK1)-RICE(NK1))

               ! If cloud parcels are virtually colder than the environment,
               ! no entrainment is allowed at this level.
               need_mixing = .true.
               if(TVQU(NK1).le.TV00(I,NK1))then
                  UER(NK1)=0.0          ! updraft entrainment rate
                  UDR(NK1)=REI          ! updraft detrainment rate
                  EE2=0.
                  UD2=1.
                  EQFRC(NK1)=0.         ! fraction of env. air for neutral stability
                  need_mixing = .false.
               endif

               ! Compute mixing if entrainment is allowed at this level
               ENV_MIXING: if (need_mixing) then

                  ! Set up level and properties for mixing
                  LET=NK1
                  TTMP=TVQU(NK1)

                  ! Determine the critical mixed fraction of updraft and environmental air.
                  ! First suppose a subparcel with 5% updraft air and 95% environmental air.
                  F1=0.95
                  F2=1.-F1

                  ! Properties of the subparcel
                  THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)     ! equiv. pot. temperature
                  RTMP=F1*Q00(I,NK1)+F2*QU(NK1)            ! specific humidity
                  TMPLIQ=F2*RLIQ(NK1)                      ! liquid water content
                  TMPICE=F2*RICE(NK1)                      ! ice content

                  ! Find the temperature and new liquid for this temporary subparcel.
                  call TPMIX(PP0(I,NK1),THTTMP,TTMP,RTMP,TMPLIQ,TMPICE,QNEWLQ, &
                       QNEWIC,RATIO2(NK1),RL,XLV0,XLV1,XLS0,XLS1, &
                       ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)

                  ! Compute moist virtual temperature of the subparcel
                  TU95=TTMP*(1.+0.608*RTMP-TMPLIQ-TMPICE)

                  ! If the parcel is buoyant with these mixing proportions, then the
                  ! critical fraction of environmental air is close to 1.  (In fact,
                  ! we set it to 1. in this case). No detrainment in this case,
                  ! only entrainment.
                  need_buoyancy_sorting = .true.
                  if(TU95.gt.TV00(I,NK1))then
                     EE2=1.
                     UD2=0.
                     EQFRC(NK1)=1.0
                     need_buoyancy_sorting = .false.
                  endif

                  ! Calculate the critical fraction of environmental air by creating
                  ! a supbarcel with 10% updraft and 90* environmental air.
                  BUOYANCY_SORTING: if (need_buoyancy_sorting) then
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

                     ! It is possible, from this subparcel, to evaluate the critical fraction EQFRC.
                     TVDIFF = TU10-TVQU(NK1)  ! Can be 0 when updraft and env temps. are very close
                     if(TVDIFF.ge.0.0) then
                        EQFRC(NK1)=0.0
                     else
                        EQFRC(NK1)=(TV00(I,NK1)-TVQU(NK1))*F1/(min(TVDIFF,-eps))
                        EQFRC(NK1)=AMAX1(0.0,EQFRC(NK1))
                        EQFRC(NK1)=AMIN1(1.0,EQFRC(NK1))
                     endif

                     ! In the simple cases in which the subparcel mixing possibilities are
                     ! all buoyant or all non-buoyant, no need to integrate the
                     ! Gaussian distribution.
                     FRACTIONAL_MIXING: if (EQFRC(NK1).eq.1) then
                        ! All buoyant case
                        EE2=1.
                        UD2=0.
                     elseif(EQFRC(NK1).eq.0.)then
                        ! All non-buoyant case
                        EE2=0.
                        UD2=1.
                     else
                        ! Integrate over the Gaussian distribution to determine the
                        ! fractional entrainment and detrainment rates.
                        call PROF5(EQFRC(NK1),EE2,UD2)
                     endif FRACTIONAL_MIXING
                  endif BUOYANCY_SORTING

                  ! Set first-level "previous" entrainment/detrainment rates
                  if(NK.eq.KLCLM1)then
                     EE1=1.        ! entrainment rate from the previous level
                     UD1=0.        ! detrainment rate from the previuos level
                  endif

                  ! Net entrainment and detrainment rates are given by the average
                  ! fractional values in the layer.
                  UER(NK1)=0.5*REI*(EE1+EE2)
                  UDR(NK1)=0.5*REI*(UD1+UD2)

               endif ENV_MIXING

               ! If the calculated updraft detrainment rate is greater than the
               ! total updraft mass flux, all cloud mass detrains.  Exit updraft
               ! calculations.
               if(UMF(NK)-UDR(NK1).lt.10.)then
                  ! If the calculated detrained mass flux is greater than the
                  ! total updraft flux, impose total detrainment of updraft
                  ! mass at the previous model level.
                  if (UDLBE > 0.) ABE=ABE-UDLBE*GRAV
                  LET=NK
                  exit UPDRAFT
               endif

               ! New updraft mass flux (remove the part that was detrained - and
               ! add the part that was entrained).
               EE1=EE2
               UD1=UD2
               UPOLD=UMF(NK)-UDR(NK1)
               UPNEW=UPOLD+UER(NK1)
               UMF(NK1)=UPNEW
               DETLQ(NK1)=RLIQ(NK1)*UDR(NK1)    ! liquid mass detrained from updraft
               DETIC(NK1)=RICE(NK1)*UDR(NK1)    ! ice mass detrained from updraft
               QDT(NK1)=QU(NK1)

               ! New properties of the updraft after entrainment processes
               QU(NK1)=(UPOLD*QU(NK1)+UER(NK1)*Q00(I,NK1))/UPNEW
               THETEU(NK1)=(THETEU(NK1)*UPOLD+THETEE(NK1)*UER(NK1))/UPNEW
               RLIQ(NK1)=RLIQ(NK1)*UPOLD/UPNEW
               RICE(NK1)=RICE(NK1)*UPOLD/UPNEW

               ! END OF ENTRAINMENT/DETRAINMENT CALCULATIONS

               ! Compute precipitation generation (fallout)
               if (abs(RATIO2(NK1)-1.).gt.1.E-6) then
                  KFRZ=NK1                         ! highest level with liq. condensate generation
               endif
               PPTLIQ(NK1)=QLQOUT(NK1)*UMF(NK)     ! generation (fallout) of liquid precipitation
               PPTICE(NK1)=QICOUT(NK1)*UMF(NK)     ! generation (fallout) of solid precipitation
               TRPPT=TRPPT+PPTLIQ(NK1)+PPTICE(NK1) ! total generation rate up to this level
               if (NK1.le.KPBL) UER(NK1)=UER(NK1)+VMFLCL*DPP(i,NK1)/DPTHMXG(I)

               ! End of updraft calculations.
            enddo UPDRAFT

            ! Check cloud depth.  If the cloud is deep enough, estimate the
            ! equilibrium temperature level (ETL) and adjust the mass flux
            ! profile at cloud top so that the mass flux decreases to zero
            ! as a linear function of pressure difference the LET and cloud top.
            LTOP=NK     ! level below the one where vert. vel. is negative

            ! If the cloud top height is less than the specified minimum height,
            ! go back to the next highest mixed layer to see if a bigger
            ! cloud can be obtained from the parcel.
            CLDHGT=Z0G(I,LTOP)-ZLCL
            if(CLDHGT < mid_depth .or. ABE < 1.)then
               do NK=KLCLM1,LTOP
                  UMF(NK)=0.
                  UDR(NK)=0.
                  UER(NK)=0.
                  DETLQ(NK)=0.
                  DETIC(NK)=0.
                  PPTLIQ(NK)=0.
                  PPTICE(NK)=0.
               enddo
               cycle TRIGGER
            endif

            ! If we made it here, the column is active
            active(i) = 1.

         enddo TRIGGER

         ! Confirm that the column is active
         if (nint(active(i)) /= 1) cycle MAIN_I_LOOP

         ! ABOVE THE LET

         ! If the LET and LTOP are the same, detrain all of the updraft
         ! mass at this level.
         need_above_let = .true.
         if(LET.eq.LTOP)then
            UDR(LTOP)=UMF(LTOP)+UDR(LTOP)-UER(LTOP)
            DETLQ(LTOP)=RLIQ(LTOP)*UDR(LTOP)*UPNEW/UPOLD
            DETIC(LTOP)=RICE(LTOP)*UDR(LTOP)*UPNEW/UPOLD
            UER(LTOP)=0.
            UMF(LTOP)=0.
            need_above_let = .false.
         endif

         ! Calculations for layers above the LET
         ABOVE_LET: if (need_above_let) then

            ! Begin total detrainment at the level above the LET.
            DPTT=0.
            do NJ=LET+1,LTOP
               DPTT=DPTT+DPP(I,NJ)
            enddo
            DUMFDP=UMF(LET)/DPTT

            ! Adjust mass flux profiles, detrainment rates, and precipitation
            ! fall rates to reflect the linear decrease in mass flux between
            ! the LET and cloud top.
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
         endif ABOVE_LET

         ! BELOW CLOUD BASE

         ! Initialize updraft properties below cloud base
         do NK=1,LC-1
            TU(NK)=0.
            QU(NK)=0.
            TVU(NK)=0.
            UMF(NK)=0.
            WU(NK)=0.
            UER(NK)=0.
         enddo

         ! Between the departure and cloud base levels create a linear
         ! profile of the updraft properties
         UMF(LC)=VMFLCL*DPP(i,LC)/DPTHMXG(I)
         UER(LC)=UMF(LC)
         do NK=LC+1,KPBL
            UER(NK)=VMFLCL*DPP(i,NK)/DPTHMXG(I)
            UMF(NK)=UMF(NK-1)+UER(NK)
         enddo
         do NK=KPBL+1,KLCLM1
            UMF(NK)=VMFLCL
            UER(NK)=0.
         enddo
         do NK=LC,KLCLM1
            TU(NK)=TMIXG(I)+(Z0G(I,NK)-ZMIX)*GDRY
            QU(NK)=QMIXG(I)
            TVU(NK)= TU(NK)*(1.+0.608*QU(NK))
            WU(NK)=WLCL
         enddo

         ! Between the departure and cloud base levels.
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

         ! Set index values surrounding cloud top
         LTOP1=LTOP+1
         LTOPM1=LTOP-1

         ! Make sure that the updraft vertical motion is not too small so that
         ! the cloud fraction becomes not too large just above the lifting
         ! condensation level (first level in the cloud above which WU is
         ! different from WLCL)
         nk = lc+1
         do while (wu(nk) <= wlcl .and. nk <= ltopm1)
            wu(nk) = max(wlcl,wu(nk))
            nk = nk+1
         enddo

         ! Updraft horizontal motions initially match the environment
         UU(LC)=U00(I,LC)
         VU(LC)=V00(I,LC)
         UU(LC+1)=UU(LC)
         VU(LC+1)=VU(LC)

         ! Adjust updraft horizontal motions above the LET for pressure effects
         do NK = LC+1,LTOPM1
            lambda=2.
            if(nk .gt. let) then
               nudr=(lambda+1.)*UDR(NK)
               nuer=uer(nk)+lambda*UDR(NK)
            else
               nudr=UDR(NK)
               nuer=uer(nk)
            endif
            UPOLD=UMF(NK-1)-nudr
            UPNEW=UPOLD+nuer
            UU(NK+1)=(UU(NK)*UPOLD+U00(I,NK)*nuer)/UPNEW
            VU(NK+1)=(VU(NK)*UPOLD+V00(I,NK)*nuer)/UPNEW
         enddo

         ! Define variables above cloud top.
         umf(ltop1:kx) = 0.
         udr(ltop1:kx) = 0.
         uu(ltop1:kx) = 0.
         vu(ltop1:kx) = 0.
         uer(ltop1:kx) = 0.
         qdt(ltop1:kx) = 0.
         rliq(ltop1:kx) = 0.
         rice(ltop1:kx) = 0.
         qlqout(ltop1:kx) = 0.
         qicout(ltop1:kx) = 0.
         detlq(ltop1:kx) = 0.
         detic(ltop1:kx) = 0.
         pptliq(ltop1:kx) = 0.
         pptice(ltop1:kx) = 0.
         thta0(ltop1:kx) = 0.
         thtau(ltop1:kx) = 0.
         qlg(ltop1:kx) = 0.
         qig(ltop1:kx) = 0.
         dtfm(ltop1:kx) = 0.
         omga(ltop1:kx) = 0.
         thadv(ltop1:kx) = 0.
         qadv(ltop1:kx) = 0.
         omg(ltop1:kx) = 0.
         ems(ltop1:kx) = dpp(i,ltop1:kx)*dxsq/GRAV
         emsd(ltop1:kx) = 1./ems(ltop1:kx)
         ug(ltop1:kx) = u00(i,ltop1:kx)
         vg(ltop1:kx) = v00(i,ltop1:kx)
         tg(ltop1:kx) = tt0(i,ltop1:kx)
         qg(ltop1:kx) = q00(i,ltop1:kx)
         qtg(ltop1:kx) = qg(ltop1:kx) + qlg(ltop1:kx) + qig(ltop1:kx)
         do nk=ltop1,kx
            nup(nk) = nk
            if (nk > ltop1) then
               tu(nk) = 0.
               qu(nk) = 0.
               tvu(nk) = 0.
               wu(nk) = 0.
            endif
         enddo
         omg(kx+1) = 0.
         p165 = pp0(i,klcl) - 1.65E4

         ! Initialize some other variables to be used later in the vertical advection
         ! scheme (from the surface to the cloud top).
         do NK=1,LTOP
            EMS(NK)=DPP(I,NK)*DXSQ/GRAV
            EMSD(NK)=1./EMS(NK)
            DTFM(NK)=0.
            EXN(NK)=(P00/PP0(I,NK))**(0.2854*(1.-0.28*QDT(NK)))
            THTAU(NK)=TU(NK)*EXN(NK)
            EXN(NK)=(P00/PP0(I,NK))**(0.2854*(1.-0.28*Q00(I,NK)))
            THTA0(NK)=TT0(I,NK)*EXN(NK)
            if(PP0(I,NK).gt.P165)LVF=NK
            QTDT(NK) = QDT(NK)+RLIQ(NK)+RICE(NK) !updraft total water
            OMG(NK)=0.
         enddo
         CLVF=0.
         do NK = KLCL,LVF+1
            CLVF=CLVF+PPTLIQ(NK)+PPTICE(NK)
         enddo
         USR=UMF(LVF+1)*QU(LVF+1)+CLVF
         USR=AMIN1(USR,TRPPT)

         ! Precipitation efficiency estimates
         if (mid_peff_const >= 0.) then
            ! Precipitation efficiency described by Market et al. (2003; WAF; Precipitation
            ! Efficiency of Warm-Season Midwestern Mesoscale Convective Systems).  A value
            ! of ~0.25 is suggested by Fig. 12b of Fankhauser (1998; WRF; Estimates of
            ! Thuderstorm Precipitation Efficiency from Field Measurements in CCOPE) for
            ! elevated convection
            peff = mid_peff_const
         else if (mid_peff == 'BLUESTEIN83') then
            ! Precipitation efficiency computed from the precipitable water deficit in the
            ! subcloud layer (Bluestein and Parks 1983; MWR)
            do k=1,kx
               pp0c(:,k) = pp0(i,k)
               q00c(:,k) = q00(i,k) * RGASD*tv00(i,k)/pp0(i,k) / GRAV
               qst1c(:,k) = qst1(i,k) * RGASD*tv00(i,k)/pp0(i,k) / GRAV
            enddo
            if (int_profile(subcloud_pw,q00c,pp0c,(/plclg(i)/),(/ps(i)/)) /= INT_OK .or. &
                 int_profile(subcloud_pws,qst1c,pp0c,(/plclg(i)/),(/ps(i)/)) /= INT_OK) then
               call physeterror('kfmid', 'Cannot compute precipitable water integrals')
               return
            endif
            peff = 0.6 * subcloud_pw(1) / subcloud_pws(1)
         else if (mid_peff == 'CONDSE') then
            ! Precipitation efficiency calculated as the sedimentation efficiency of
            ! the gridscale condensation scheme
            if (any(cond_prflux(i,:) > 0.)) then
               cond_maxlev = maxloc(cond_prflux(i,:))
               peff = PEFFMAX*minval(cond_prflux(i,cond_maxlev(1):kx)) / cond_prflux(i,cond_maxlev(1))
            else
               peff = 0.4  ! Set a default value if the gridscale scheme is inactive
            endif
         else if (mid_peff == 'TLCL') then
            peff = max(min((tenv-200.)/90.,PEFFMAX),0.)
         else if (mid_peff == 'QLCL') then
            peff = max(min(0.2+60.*(qenv-0.006),PEFFMAX),0.)
         endif


         ! 10. DOWNDRAFT PROPERTIES
         ! ========================

         ! Let the representative downdraft originate at the level of minimum saturation
         ! equivalent potential temperature (SEQT) in the cloud layer.  Extend downward
         ! toward the surface, down to the level of downdraft buoyancy (LDB), where SEQT
         ! is less than min(SEQT) in the cloud layer.  Let downdraft detrain over a
         ! layer of specified pressure-depth (i.e., DPDD).
         TDER = 0.

         ! KMIN is the level of minimum environmental saturation equivalent potential
         ! temperature.  If the level KMIN is higher than the cloud top then track again
         ! the level between LCL and LTOP.
         if (KMIN >= LTOP) then
            THTMIN=THTES(KLCL)
            KMIN=KLCL
            do NK=KLCL+1,LTOP-1
               THTMIN=AMIN1(THTMIN,THTES(NK))
               if (THTMIN == THTES(NK)) KMIN=NK
            enddo
         endif
         LFS=KMIN

         ! Compute environmental properties at the LFS
         if(RATIO2(LFS).gt.0.) then
            call ENVIRTHT(PP0(I,LFS),TT0(I,LFS),Q00(I,LFS),THETEE(LFS),0.,RL, &
                 ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
         endif

         ! The EQFRC fraction again represents the fraction
         ! of environmental air in mixed subparcels
         ! (from updraft and environmental) at the LFS.
         EQFRC(LFS)=(THTES(LFS)-THETEU(LFS))/(THETEE(LFS)-THETEU(LFS)+1.E-20)
         EQFRC(LFS)=AMAX1(EQFRC(LFS),0.)
         EQFRC(LFS)=AMIN1(EQFRC(LFS),1.)

         ! Eq. pot. temperature of the downdraft at the LFS (saturated).
         THETED(LFS)=THTES(LFS)

         ! Adjust precipitation efficiency using downdraft properties on request
         if (mid_peff == 'TDD') then
            peff = max(min(0.15+0.015*(theted(lfs)-310.),0.9),0.)
         endif
         mpeff(i) = peff

         ! Environmental temperature at the LFS
         TZ(LFS)=TT0(I,LFS)

         ! Recalculate the eq. pot. temperature using this new temperature.
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

         ! Find the level of downdraft buoyancy (LDB) ==> level where THETAES
         ! for the downdraft at the LFS becomes larger then the environmental
         ! THETAES.
         ldb = 1
         FIND_LDB: do nk=1,lfs-1
            nd = lfs-nk
            if (theted(lfs) > thtes(nd) .or. nd == 1) then
               ldb = nd
               exit FIND_LDB
            endif
         enddo FIND_LDB

         ! Determine the LDT level (level where the downdraft detrains).  This
         ! level cannot be higher than the downdraft source level (LFS).
         DPDDMX=mid_dpdd * dpddmult(i)
         DPT=0.
         NK=LDB
         LDT=-1
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
         TVD(LFS)=TT0(I,LFS)*(1.+0.608*QST1(I,LFS))  ! moist virtual temperature
         RDD=PP0(I,LFS)/(RGASD*TVD(LFS))             ! air density
         A1=(1.-PEFF)*AU0
         DMF(LFS)=-A1*RDD                            ! downdraft mass flux (negative)
         DER(LFS)=EQFRC(LFS)*DMF(LFS)                ! downdraft entrainment rate
         DDR(LFS)=0.                                 ! downdraft detrainment rate

        ! Compute the downdraft mass flux
         DOWNDRAFT: do ND=LFS-1,LDB,-1
            ND1=ND+1
            ! At the lowest level
            if(ND.le.LDT)then
               DER(ND)=0.
               DDR(ND)=-DMF(LDT+1)*DPP(I,ND)*FRC/DPDD
               DMF(ND)=DMF(ND1)+DER(ND)+DDR(ND)
               FRC=1.
               THETED(ND)=THETED(ND1)
               QD(ND)=QD(ND1)
               TZ(ND)=TPDD1(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
                    XLV0,XLV1, &
                    ALIQ,BLIQ,CLIQ,DLIQ)
               EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
               THTAD(ND)=TZ(ND)*EXN(ND)
            else
               ! Before reaching its lowest level, downdraft entrains, but
               ! does not detrain.
               DER(ND)=DMF(LFS)*0.03*DPP(I,ND)/RAD
               DDR(ND)=0.
               if (ND <= klcl) DDR(ND) = -DMF(LFS)*0.03*DPP(I,ND)/RAD
               DMF(ND)=DMF(ND1)+DER(ND)+DDR(ND)

               ! Recalculate THETEE
               if(RATIO2(ND).gt.0.) then
                  call ENVIRTHT(PP0(I,ND),TT0(I,ND),Q00(I,ND), &
                       THETEE(ND),0.,RL, &
                       ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
               endif
               ! Downdraft properties deduced from properties
               ! at the previous level and entrained air.
               THETED(ND)=(THETED(ND1)*DMF(ND1)+THETEE(ND)*(DER(ND)+DDR(ND)))/DMF(ND)
               QD(ND)=(QD(ND1)*DMF(ND1)+Q00(I,ND)*(DER(ND)+DDR(ND)))/DMF(ND)
               TZ(ND)=TPDD1(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
                    XLV0,XLV1, &
                    ALIQ,BLIQ,CLIQ,DLIQ)
               EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
               THTAD(ND)=TZ(ND)*EXN(ND)
            endif
         enddo DOWNDRAFT
         TDER=0.

!!$         ! Compute the downdraft mass flux
!!$         DOWNDRAFT: do ND=LFS-1,LDB,-1
!!$            ND1=ND+1
!!$            ! At the lowest level
!!$            if (nd > ldt) then
!!$               ! Before reaching its lowest level, downdraft entrains, but
!!$               ! does not detrain.
!!$               DER(ND)=DMF(LFS)*0.03*DPP(I,ND)/RAD
!!$               DDR(ND)=0.
!!$               if (ND <= klcl) DDR(ND) = -DMF(LFS)*0.03*DPP(I,ND)/RAD
!!$               DMF(ND)=DMF(ND1)+DER(ND)+DDR(ND)
!!$               ! If all mass has left the downdraft, the base is here
!!$               if (dmf(nd) <= 10.) then
!!$                  ldt = nd
!!$                  FRC=1.
!!$                  THETED(ND)=THETED(ND1)
!!$                  QD(ND)=QD(ND1)
!!$                  TZ(ND)=TPDD(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
!!$                       XLV0,XLV1, &
!!$                       ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!!$                  EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
!!$                  THTAD(ND)=TZ(ND)*EXN(ND)
!!$               else
!!$                  ! Recalculate THETEE
!!$                  if(RATIO2(ND).gt.0.) then
!!$                     call ENVIRTHT(PP0(I,ND),TT0(I,ND),Q00(I,ND), &
!!$                          THETEE(ND),0.,RL, &
!!$                          ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!!$                  endif
!!$                  ! Downdraft properties deduced from properties
!!$                  ! at the previous level and entrained air.
!!$                  THETED(ND)=(THETED(ND1)*DMF(ND1)+THETEE(ND)*(DER(ND)+DDR(ND)))/DMF(ND)
!!$                  QD(ND)=(QD(ND1)*DMF(ND1)+Q00(I,ND)*(DER(ND)+DDR(ND)))/DMF(ND)
!!$                  TZ(ND)=TPDD(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
!!$                       XLV0,XLV1, &
!!$                       ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!!$                  EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
!!$                  THTAD(ND)=TZ(ND)*EXN(ND)
!!$               endif
!!$            else
!!$               DER(ND)=0.
!!$               DDR(ND)=-DMF(LDT+1)*DPP(I,ND)*FRC/DPDD
!!$               DMF(ND)=DMF(ND1)+DER(ND)+DDR(ND)
!!$               FRC=1.
!!$               THETED(ND)=THETED(ND1)
!!$               QD(ND)=QD(ND1)
!!$               TZ(ND)=TPDD(PP0(I,ND),THETED(ND),TT0(I,ND),QS,QD(ND),1.0, &
!!$                    XLV0,XLV1, &
!!$                    ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!!$               EXN(ND)=(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
!!$               THTAD(ND)=TZ(ND)*EXN(ND)
!!$            endif
!!$         enddo DOWNDRAFT
!!$         TDER=0.

         ! Effect of downdraft reaching the LD B(cooling and moistening).  Relative
         ! humidity of 90% at this level.
         DOWNDRAFT_DETRAIN: do ND=LDB,LDT
            TZ(ND)=TPDD1(PP0(I,ND),THETED(LDT),TT0(I,ND),QS,QD(ND),1.0, &
                 XLV0,XLV1, &
                 ALIQ,BLIQ,CLIQ,DLIQ)
            ES = ALIQ*exp((TZ(ND)*BLIQ-CLIQ)/(TZ(ND)-DLIQ))
            QS = 0.622*ES/(PP0(I,ND)-ES)
            DQSDT=(CLIQ-BLIQ*DLIQ)/((TZ(ND)-DLIQ)*(TZ(ND)-DLIQ))

            ! Adjust temperature and humidity for a relative humidity of 90%.
            RL=XLV0-XLV1*TZ(ND)
            DTMP=RL*QS*(1.-RHBC)/(cpd+RL*RHBC*QS*DQSDT)
            T1RH=TZ(ND)+DTMP
            ES=RHBC*ALIQ*exp((BLIQ*T1RH-CLIQ)/(T1RH-DLIQ))
            ES=min( ES , 0.5*PP0(I,ND) )
            QSRH=0.622*ES/(PP0(I,ND)-ES)
            QSRH=max( min( QSRH , 0.050 ) , 1.E-6 )

            ! Check to see if the mixing ratio at specified relative humidity is less
            ! than the actual mixing ratio.  If so, adjust to give zero evaporation.
            if(QSRH.lt.QD(ND))then
               QSRH=QD(ND)
               T1RH=TZ(ND)
            endif

            ! If this is not the case, use the adjusted values for temperature
            ! and humidity
            TZ(ND)=T1RH
            QS=QSRH
            TDER=TDER + (QS-QD(ND))*DDR(ND)
            QD(ND)=QS
            THTAD(ND)=TZ(ND)*(P00/PP0(I,ND))**(0.2854*(1.-0.28*QD(ND)))
         enddo DOWNDRAFT_DETRAIN

         ! Deactivate downdraft on user reqest
         if (.not.USE_DOWNDRAFT) TDER = 0.

         ! Compute the downdraft mass flux at the LFS
         ! PPTFLX = precipitation flux =
         ! updraft supply rate times the
         ! precipitation efficiency
         ! RCED   = Rate of condensate evaporation
         ! in downdraft =
         ! total precipitation rate -
         ! precipitation flux
         ! (this is the precip that evaporates
         ! in downdraft).
         ! PPR    = total fallout of liquid water and
         ! and ice in the updraft.
         ! DMFLFS = adjusted downward mass flux at
         ! the LFS
         ! CNDTNF = condensate in updraft air at the LFS.
         ! DPPTDF = total fallout precip rate in updraft.
         DOWNDRAFT_LFS: if (TDER >= 1.) then
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

            ! No downdraft at the LFS
            if(DMFLFS.gt.0.)then
               TDER=0.
            endif
         endif DOWNDRAFT_LFS

         ! Adjust downdraft mass flux so that evaporation rate in downdraft is
         ! consistent with precipitation efficiency relationship
         DOWNDRAFT_ADJUST: if (tder >= 1.) then
            DDINC=DMFLFS/DMF(LFS)                                    ! downdraft increase
            if(LFS.ge.KLCL)then
               UPDINC=(UMF(LFS)-(1.-EQFRC(LFS))*DMFLFS)/UMF(LFS)     ! updraft increase
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

            ! Adjust upward mass flux, mass detrainment rate, and liquid water
            ! detrainment rates to be consistent with the transfer of the
            ! estimate from the updraft to the downdraft at the LFS.
            do NK=LC,LFS
               UMF(NK)=UMF(NK)*UPDINC
               UDR(NK)=UDR(NK)*UPDINC
               UER(NK)=UER(NK)*UPDINC
               PPTLIQ(NK)=PPTLIQ(NK)*UPDINC
               PPTICE(NK)=PPTICE(NK)*UPDINC
               DETLQ(NK)=DETLQ(NK)*UPDINC
               DETIC(NK)=DETIC(NK)*UPDINC
            enddo

            ! Set values below the downdraft buoyancy level
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

            ! Set values above the LFS
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

            ! Downdraft wind adjustments
            UD(LFS)=EQFRC(LFS)*U00(I,LFS)+(1.-EQFRC(LFS))*UU(LFS+1)
            VD(LFS)=EQFRC(LFS)*V00(I,LFS)+(1.-EQFRC(LFS))*VU(LFS+1)
            do NK=1,LFS-LDB
               NJ=LFS-NK
               if(NJ.gt.LDT)then
                  UD(NJ)=(UD(NJ+1)*DMF(NJ+1)+U00(I,NJ)*DER(NJ))/DMF(NJ)
                  VD(NJ)=(VD(NJ+1)*DMF(NJ+1)+V00(I,NJ)*DER(NJ))/DMF(NJ)
               else
                  UD(NJ)=UD(LDT+1)
                  VD(NJ)=VD(LDT+1)
               endif
            enddo

         else

            ! If there is no evaporation, then create no downdraft
            PPTFLX=TRPPT             ! precipitation flux
            CPR=TRPPT                ! condensate production rate
            TDER=0.                  ! total downdraft evaporation rate (?)
            CNDTNF=0.
            UPDINC=1.
            LDB=LFS
            do NDK=1,LTOP
               UD(NDK)=0.
               VD(NDK)=0.
               DMF(NDK)=0.
               DER(NDK)=0.
               DDR(NDK)=0.
               THTAD(NDK)=0.
               TZ(NDK)=0.
               QD(NDK)=0.
            enddo
            AINCM2=100.
            DMFMIN=0.

         endif DOWNDRAFT_ADJUST

         ! Set limits on the updraft and downdraft mass fluxes so that the influx
         ! into convective drafts from a given layer is no more than is available
         ! in that layer initially.  Also, do not allow updraft detrainment to
         ! exceed twice the mass in a layer initially.  Or downdraft detrainment
         ! exceed one time the initial mass in a layer.
         AINCMX=DXDY(I)/(PI*RAD*RAD)
         LMAX=MAX(KLCL,LFS)
         do NK=LC,LMAX
            if((UER(NK)-DER(NK)).gt.0.) then
               AINCM1=EMS(NK)/((UER(NK)-DER(NK))*TIMEC)  ! note that EMS is the total mass in layer
            endif
            AINCMX=MIN(AINCMX,AINCM1)
         enddo

         ! If the entrainment rates for the updraft are reasonable physically
         ! speaking, then we should have AINCMX > 1, thus AINC=1 most of the time
         ! here.  If, on the other hand, the entrainment rates are so large so that
         ! AINCMX < 1, then AINC is equal to this fraction.  (This correction is done
         ! to avoid removing more air than available).

         ! Prescribed total cloud base mass flux for the cell as closure
         AINC = RHOLCL*DXDY(I)*WLCL / VMFLCL
         AINC = MIN(AINC,AINCMX)

         ! Adjust closure based using cloud properties on request
         if (mid_emffrac == 'TLCL') then
            ainc = max(min((280.-tenv)/80.,1.),0.) * ainc
         else if (mid_emffrac == 'TDD') then
            ainc = max(min((360.-theted(lfs))/50.,1.),0.) * ainc
         else if (mid_emffrac == 'QLCL') then
            ainc = max(min(125.*(0.018-qenv),1.),0.) * ainc
         endif
         mainc(i) = ainc

         ! Save the relevant variables for a unit updraft and downdraft.  They are
         ! adjusted by the factor AINC to satisfy the closure

         ! Adjustment of all the updraft/downdraft characteristics.
         TDER=TDER*AINC
         PPTFLX=PPTFLX*AINC
         do K=1,LTOP
            UMF(K)=UMF(K)*AINC
            DMF(K)=DMF(K)*AINC
            DETLQ(K)=DETLQ(K)*AINC
            DETIC(K)=DETIC(K)*AINC
            UDR(K)=UDR(K)*AINC
            UER(K)=UER(K)*AINC
            DER(K)=DER(K)*AINC
            DDR(K)=DDR(K)*AINC
            PPTLIQ(K)=PPTLIQ(K)*AINC
            PPTICE(K)=PPTICE(K)*AINC
            DTFM(K)=0.
         enddo

         ! If the downdraft overdraws the initial mass available at the LFS,
         ! allow PEFF to change so that updraft mass flux is not the limiting
         ! factor in the total convective mass flux.  This is usually necessary
         ! only when the downdraft is very shallow.
         DMFMIN=UER(LFS)-EMS(LFS)/TIMEC

         ! Handle the case of very large (negative) downdraft entrainment rate
         LARGE_DER: if(DER(LFS).lt.DMFMIN .and. DMFMIN.lt.0.)then

            ! RF is the proportionality constant for the decrease of downdraft
            ! fluxes.  (RF < 1).
            RF=DMFMIN/DER(LFS)
            do K=LDB,LFS
               DER(K)=DER(K)*RF
               DDR(K)=DDR(K)*RF
               DMF(K)=DMF(K)*RF
            enddo
            TDER = TDER*RF

            ! Modify the UPDINC factor using the new downdraft fluxes
            if(LFS.ge.KLCL)then
               UPDIN2=1.-(1.-EQFRC(LFS))*DMF(LFS)*UPDINC/UMF(LFS)
            else
               UPDIN2=1.
            endif

            ! Recalculate the condensate production rate (CPR) and
            ! precipitation efficiency.
            CPR=TRPPT+PPR*(UPDIN2-1.)
            PEFF=1.-(TDER+(1.-EQFRC(LFS))*DMF(LFS)*(RLIQ(LFS)+RICE(LFS)))/ &
                 (CPR*AINC)
            PEFF=MIN(MAX(PEFF,0.),1.)

            ! Adjust precipitation fluxes.
            PPTFLX=PEFF*CPR*AINC

            ! Updraft characteristics have to change accordingly.
            F1=UPDIN2/UPDINC
            do K=LC,LFS
               UMF(K)=UMF(K)*F1
               UDR(K)=UDR(K)*F1
               UER(K)=UER(K)*F1
               DETLQ(K)=DETLQ(K)*F1
               DETIC(K)=DETIC(K)*F1
               PPTLIQ(K)=PPTLIQ(K)*F1
               PPTICE(K)=PPTICE(K)*F1
            enddo

         endif LARGE_DER

         ! Evaluate the vertical velocity OMGA (Pa/s) from the entrainment
         ! and detrainment rates.  (vertical velocity for the environment).
         DTT=TIMEC
         do NK = 1, LTOP
            DOMGDP(NK)=-(UER(NK)-DER(NK)-UDR(NK)-DDR(NK))*EMSD(NK)   ! D(omega)/Dp
            if (NK == 1) then
               OMG(NK)=(UMF(NK)+DMF(NK))*GRAV/DXSQ
            else
               OMG(NK)=OMG(NK-1)-DPP(I,NK-1)*DOMGDP(NK-1)
               DTT1 = 0.75*DPP(I,NK-1)/(abs(OMG(NK))+1.E-10)
               DTT=AMIN1(DTT,DTT1)                                   ! mass flux timestep
            endif
         enddo
         NSTEP=nint(TIMEC/DTT+1)
         DTIME=TIMEC/FLOAT(NSTEP)

         ! Set up state variables for advection
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

         ! Recompute updraft and downdraft mass fluxes to ensure consistency with
         ! entr/detr rates
         do nk=lc+1,ltop-1
            umf(nk) = umf(nk-1) + uer(nk) - udr(nk)
         enddo
         udr(ltop) = umf(ltop-1)
         umf(ltop) = 0.
         do nk=lfs-1,ldb,-1
            dmf(nk) = dmf(nk+1) - (-der(nk)) + ddr(nk)
         enddo

         ! Compute environmental mass flux by mass conservation
         do nk=1,ltop
            emf(nk) = -(umf(nk)+dmf(nk))
         enddo

         ! Implicit solution for mass flux equations
         IMPLICIT_MF: do NTC=1,NSTEP

            ! Mass flux equations below cloud top
            do  NK = 1, LTOPM1

               ! Set up coefficients for matrix solution using upwind advection
               ! for numerical stability
               if (omga(nk) >= 0. .or. nk == 1) then
                  delp = pp0(i,nk+1)-pp0(i,nk)
                  ta(nk) = 0.
                  tc(nk) = -dtime*(omga(nk)/delp)
               else
                  delp = pp0(i,nk-1)-pp0(i,nk)
                  ta(nk) = -dtime*(omga(nk)/delp)
                  tc(nk) = 0.
               endif
               tb(nk) = dtime*(omga(nk)/delp - emsd(nk)*(udr(nk) + ddr(nk)))

               ! Set up variable-specific forcing terms
               td_thta(nk) = dtime*emsd(nk)*(udr(nk)*thtau(nk) + ddr(nk)*thtad(nk))
               td_q(nk) = dtime*emsd(nk)*(udr(nk)*qdt(nk) + ddr(nk)*qd(nk))
               td_qt(nk) = dtime*emsd(nk)*(udr(nk)*qtdt(nk) + ddr(nk)*qd(nk))
               td_ql(nk) = dtime*emsd(nk)*(udr(nk)*rliq(nk))
               td_qi(nk) = dtime*emsd(nk)*(udr(nk)*rice(nk))

            enddo

            ! Coefficients for cloud top
            ta(ltop) = 0.
            tb(ltop) = -dtime*emsd(ltop)*udr(ltop)
            tc(ltop) = 0.

            ! Consider only updraft detrainment at cloud top
            td_thta(ltop) = dtime*emsd(ltop)*(udr(ltop)*thtau(ltop))
            td_q(ltop) = dtime*emsd(ltop)*(udr(ltop)*qdt(ltop))
            td_qt(ltop) = dtime*emsd(ltop)*(udr(ltop)*qtdt(ltop))
            td_ql(ltop) = dtime*emsd(ltop)*(udr(ltop)*rliq(ltop))
            td_qi(ltop) = dtime*emsd(ltop)*(udr(ltop)*rice(ltop))

            ! Set up tri-diagonal (D) matrices for implicit solution of mass flux equations
            call difuvd1(tri_thta,1.,ta,tb,tc,thpa,td_thta,1,1,ltop)
            call difuvd1(tri_q,1.,ta,tb,tc,qpa,td_q,1,1,ltop)
            call difuvd1(tri_qt,1.,ta,tb,tc,qtpa,td_qt,1,1,ltop)
            call difuvd1(tri_ql,1.,ta,tb,tc,qlpa,td_ql,1,1,ltop)
            call difuvd1(tri_qi,1.,ta,tb,tc,qipa,td_qi,1,1,ltop)

            ! Adjust coefficients to (I-D) for tendency solution
            ta(1:ltop) = -ta(1:ltop)
            tb(1:ltop) = 1.-tb(1:ltop)
            tc(1:ltop) = -tc(1:ltop)

            ! Matrix solution for increments to variables predicted by mass flux equations
            call difuvd2(thpai,ta,tb,tc,tri_thta,workk,1,1,ltop)
            call difuvd2(qpai,ta,tb,tc,tri_q,workk,1,1,ltop)
            call difuvd2(qtpai,ta,tb,tc,tri_qt,workk,1,1,ltop)
            call difuvd2(qlpai,ta,tb,tc,tri_ql,workk,1,1,ltop)
            call difuvd2(qipai,ta,tb,tc,tri_qi,workk,1,1,ltop)

            ! Update environment properties
            do nk=1,ltop
               thpa(nk) = thpa(nk) + thpai(nk)
               qpa(nk) = qpa(nk) + qpai(nk)
               qtpa(nk) = qtpa(nk) + qtpai(nk)
               qlpa(nk) = qlpa(nk) + qlpai(nk)
               qipa(nk) = qipa(nk) + qipai(nk)
            enddo

         enddo IMPLICIT_MF

         ! Specify the "G" grid values after convective adjustment
         do NK=1,LTOP
            THTAG(NK)=THPA(NK)
            QLG(NK)=QL0(I,NK)
            QIG(NK)=QI0(I,NK)
            QTG(NK)=Q00(I,NK)+QLG(NK)+QIG(NK)
            QG(NK)=MAX(QPA(NK),0.)
         enddo

         ! Convert THETA to T, freeze all supercooled detrained liquid water
         ! because no supercooled water is supposed to be allowed in the
         ! explicit condensation scheme (anyone of those available from
         ! the subroutine VKUOCON).
         do NK=1,LTOP
            EXN(NK)=(P00/PP0(I,NK))**(0.2854*(1.-0.28*QG(NK)))
            TG(NK)=THTAG(NK)/EXN(NK)
            ! Temperature effect of melting
            if (tt0(i,nk) < TRPL) tg(nk) = tg(nk) + detlq(nk)*CHLF*emsd(nk)/CPD
         enddo

!!$       if (kount == 6) then
!!$          do k=ltop-1,1,-1
!!$             cdmf = ''
!!$             if (k==lfs) then
!!$                cdmf = '\/ \/ \/'
!!$             else if (dmf(k) < 0) then
!!$                cdmf = '   \/   '
!!$             endif
!!$             cddrl=''; cddrr=''
!!$             if (ddr(k) > 0) then
!!$                cddrl = '< '
!!$                cddrr = ' >'
!!$             endif
!!$             cderl=''; cderr=''
!!$             if (der(k) < 0) then
!!$                cderl = '>'
!!$                cderr = '<'
!!$             endif
!!$             if (k==ldt) then
!!$                cddrl = '<-'
!!$                cddrr = '->'
!!$             endif
!!$             if (k<ldt .and. k>ldb) then
!!$                cddrl = '<<'
!!$                cddrr = '>>'
!!$             endif
!!$             if (k==ldb) then
!!$                cdmf = '-- -- --'
!!$                cddrl = '<<'
!!$                cddrr = '>>'
!!$             endif
!!$             write(6,'(i,x,f,3(x,e),x,a)') k,z0g(i,k),thpai(k),umf(k),dmf(k),cddrl//cderl//cdmf//cderr//cddrr
!!$             if (k == klcl+1) print*, '------------------------'
!!$          end do
!!$       endif

         ! Diagnostic cloud and condensate outputs
         do K=1,KX
            NK = KX-K+1
            RLIQOUT(I,NK)  = RLIQ(K)
            RICEOUT(I,NK)  = RICE(K)

            ! Compute updraft area
            if (WU(K).gt.0.1) then
               AREAUP(I,NK) = (UMF(K)*RGASD*TVU(K)) / (PP0(I,K)*WU(K))
            endif

            ! Convert updraft areas into cloud fractions and clip the
            ! result to eliminate negative values.
            CLOUDS (I,NK) = min(1.,max(0.,AREAUP(I,NK)/DXDY(I)))
            if (K < KLCL) CLOUDS (I,NK) = 0.

            ! Normalize condensate mixing ratio
            RLIQOUT(I,NK) = RLIQOUT(I,NK)*CLOUDS(I,NK)
            RICEOUT(I,NK) = RICEOUT(I,NK)*CLOUDS(I,NK)
         end do

         ! Convective momentum transport

         ! Initialize the momentum variables for vertical advection calculations
         do NK=1,LTOP
            UPA(NK)=U00(I,NK)
            VPA(NK)=V00(I,NK)
         enddo

         ! Determine the type of numerical scheme to use for the momentum mass flux
         ! equations, used to compute new values for the prognostic wind variables
         ! due to vertical advection of the environmental properties, as well as
         ! detrainment effects from the updraft and downdraft.
         ! Implicit solution for momentum mass flux equations
         IMPLICIT_MOM_STEP: do NTC=1,NSTEP

            ! Momentum mass flux equations below cloud top
            do  nk=1,ltopm1

               ! Set up coefficients for matrix solution using upwind advection
               ! for numerical stability
               if (omga(nk) >= 0. .or. nk == 1) then
                  delp = pp0(i,nk+1)-pp0(i,nk)
                  ta(nk) = 0.
                  tc(nk) = -dtime*(oneminc*omga(nk)/delp)
               else
                  delp = pp0(i,nk-1)-pp0(i,nk)
                  ta(nk) = -dtime*(oneminc*omga(nk)/delp)
                  tc(nk) = 0.
               endif
               tb(nk) = dtime*(oneminc*omga(nk)/delp - emsd(nk)*(udr(nk) + ddr(nk)))

               ! Set up variable-specific forcing terms
               td_u(nk) = dtime*emsd(nk)*(udr(nk)*uu(nk) + ddr(nk)*ud(nk))
               td_v(nk) = dtime*emsd(nk)*(udr(nk)*vu(nk) + ddr(nk)*vd(nk))

            enddo

            ! Coefficients for cloud top
            ta(ltop) = 0.
            tb(ltop) = -dtime*emsd(ltop)*udr(ltop)
            tc(ltop) = 0.

            ! Consider only updraft detrainment at cloud top
            td_u(ltop) = dtime*emsd(ltop)*(udr(ltop)*uu(ltop))
            td_v(ltop) = dtime*emsd(ltop)*(udr(ltop)*vu(ltop))

            ! Set up tri-diagonal (D) matrices for implicit solution of momentum mass flux equations
            call difuvd1(tri_u,1.,ta,tb,tc,upa,td_u,1,1,ltop)
            call difuvd1(tri_v,1.,ta,tb,tc,vpa,td_v,1,1,ltop)

            ! Adjust coefficients to (I-D) for tendency solution
            ta(1:ltop) = -ta(1:ltop)
            tb(1:ltop) = 1.-tb(1:ltop)
            tc(1:ltop) = -tc(1:ltop)

            ! Matrix solution for increments to variables predicted by momentum mass flux equations
            call difuvd2(upai,ta,tb,tc,tri_u,workk,1,1,ltop)
            call difuvd2(vpai,ta,tb,tc,tri_v,workk,1,1,ltop)

            ! Update environment properties
            upa(1:ltop) = upa(1:ltop) + upai(1:ltop)
            vpa(1:ltop) = vpa(1:ltop) + vpai(1:ltop)

         enddo IMPLICIT_MOM_STEP

         ! Assign updated "G" values for final state
         do NK=1,LTOP
            UG(NK)=UPA(NK)
            VG(NK)=VPA(NK)
         enddo

         ! Compute preliminary momentum tendencies
         do K=1, KX
            NK = KX-K+1
            DUDT(I,NK) = ( UG(K)-U00(I,K) ) / TIMEC
            DVDT(I,NK) = ( VG(K)-V00(I,K) ) / TIMEC
         enddo

         ! Impose CMT conservation
         ! Calculate CMT budget from cloud top to sfc
         INTDUDT = 0.
         INTDVDT = 0.
         INTDP   = 0.
         do K=KX-LTOP+1,KX
            NK = KX-K+1
            INTDUDT  = INTDUDT  + DUDT(I,K)  / GRAV * DPP(I,NK)
            INTDVDT  = INTDVDT  + DVDT(I,K)  / GRAV * DPP(I,NK)
            INTDP = INTDP + DPP(I,NK)
         end do

         ! Apply momentum tendency correction from cloud top to surface
         do K=KX-LTOP+1,KX
            DUDT(I,K) = DUDT(I,K) - INTDUDT*GRAV/INTDP
            DVDT(I,K) = DVDT(I,K) - INTDVDT*GRAV/INTDP
         end do

         ! Compute scalar tendencies
         do K = 1, KX
            NK = KX-K+1
            DTDT(I,NK) = (TG(K)-TT0(I,K))/TIMEC
            DQDT(I,NK) = (QG(K)-Q00(I,K))/TIMEC
            DQCDT(I,NK) = (QLG(K)-QL0(I,K))/TIMEC
            DQIDT(I,NK) = (QIG(K)-QI0(I,K))/TIMEC
            DQCDT(I,NK) = DETLQ(K)*EMSD(K)
            DQIDT(I,NK) = DETIC(K)*EMSD(K)
            DQCDT(I,NK) = max( DQCDT(I,NK) , 0. )
            DQIDT(I,NK) = max( DQIDT(I,NK) , 0. )
         enddo

         ! Extract diagnostics
         mcd(i) = cldhgt
         
         ! Diagnose precipitation fluxes for assimilation purposes
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

         ! Impose total water conservation by adjusting precipitation on request
         ZCRR(I) = PPTFLX/DXSQ
         PRECIP_ADJUSTMENT: if (mid_conserve == 'PRECIP') then

            ! Adjust precipitation to so that there is no total water residual
            DQDTDK = 0.
            DQCDTDK = 0.
            do K=1,KX
               NK = KX-K+1
               DQDTDK = DQDTDK  + DQDT(I,K)  / GRAV * DPP(I,NK)
               DQCDTDK = DQCDTDK + (DQCDT(I,K)+DQIDT(I,K)) / GRAV * DPP(I,NK)
            end do
            ZCRR(I) = max(0.0, -(DQDTDK + DQCDTDK))

            ! Normalize precipitation fluxes with conservative precipitation
            do K=1,KX
               NORM = ( ZCRR(I) / max(RNFLX(I,KX)+SNOFLX(I,KX),1.E-10) )
               RNFLX(I,K)  = NORM * RNFLX(I,K)
               SNOFLX(I,K) = NORM * SNOFLX(I,K)
            end do

         endif PRECIP_ADJUSTMENT

      enddo MAIN_I_LOOP

      ! Convert updraft areas into cloud fractions
      do K=1,KX
         do I=1,IX
            CLOUDS(I,K) = min(1.,max(0.,AREAUP(I,K)/DXDY(I)))
         end do
      end do

      return
   end subroutine kfmid1

end module kfmid

