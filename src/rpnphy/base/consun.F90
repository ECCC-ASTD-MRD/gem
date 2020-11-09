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

subroutine CONSUN3(STT  ,  SQT  ,  SWT  ,  SRR  ,  SSR  ,  SCF , &
     CTT  ,  CQT  ,  CWT  ,  CRR  ,  CSR  ,  CCF , &
     TP   ,  TM   ,  QP   ,  QM   ,  CWP  ,  CWM , &
     PSP  ,  PSM  , ilab  ,  DBDT ,  S    , TAU  , &
     PRFLX, SWFLX ,  F12 ,  FEVP , ICEFRAC, &
     MRK2, ni   ,  nlev)
   use tdpack
   use phy_options
   use ens_perturb, only: ens_nc2d, ens_spp_get
   implicit none
!!!#include <arch_specific.hf>

   integer ni , nlev
   real    STT(ni,nlev)    , SQT(ni,nlev)    , SWT(ni,nlev) , &
        CTT(ni,nlev)    , CQT(ni,nlev)    , CWT(ni,nlev) , &
        SRR(ni)         , SSR(ni)         , SCF(ni,nlev) , &
        CRR(ni)         , CSR(ni)         , CCF(ni,nlev) , &
        PRFLX(ni,nlev+1), SWFLX(ni,nlev+1), &
        TP(ni,nlev)     , TM(ni,nlev)     , &
        QP(ni,nlev)     , QM(ni,nlev)     , &
        CWP(ni,nlev)    , CWM(ni,nlev)    , &
        PSP(ni)         , PSM(ni)         , &
        DBDT(ni)        , S(ni,*)         , &
        TAU             , F12(ni,nlev)    , FEVP(ni,nlev), &
        ICEFRAC(ni,nlev), MRK2(ni,ens_nc2d)

   integer ilab(ni,nlev)

   !@Authors Claude Girard and Gerard Pellerin (1995)
   !@Revision
   ! 001     G.Pellerin (Fev 97) Modify onset of downdraft below 500mb
   ! 002     G.Pellerin (Nov 01) Remove negative specific humidity
   !                             generation in the stratosphere
   ! 003     S. Menard and B. Bilodeau (Feb 2003) - Add diagnostics
   !                             for AURAMS
   ! 004     L. Spacek (May 2003) version IBM
   !               - calls to vsexp routine (from massvp4 library)
   !               - calls to optimized routine mfoqst3
   ! 005     R. McTaggart-Cowan (Jul 2006) - inline function FMROFT
   ! 006     B. Bilodeau (May 2007)
   !               - clip on xpradd to avoid the generation of
   !                 infinitesimal precipitation
   ! 007     S. Gravel (July 2014) set some volatile bus variables to zero
   !@Object
   !  This routine deals with parameterization of condensation
   !  and associated generation of precipitation from liquid water
   !  a)  parameterization of stratiform condensation (SUNDQVIST flavour)
   !  b)  convective temperature and humidity changes, along with a cloud
   !      fraction, are imported from a convection parameterization scheme.
   !  c)  parameterization of precipitation generation (SUNDQVIST flavour)
   !@Arguments
   !          - Outputs -
   ! STT      large scale (stable) temperature tendency
   ! SQT      large scale (stable) specific humidity tendency
   ! SWT      large scale (stable) cloud water tendency
   ! SRR      large scale (stable) rain rate
   ! SSR      large scale (stable) snow rate
   ! SCF      large scale (stable) cloud fraction
   !          - Inputs/Outputs
   ! CTT      convective temperature tendency
   ! CQT      convective specific humidity tendency
   !          - Outputs -
   ! CWT      convective cloud water tendency
   ! CRR      convective rain rate
   ! CSR      convective snow rate
   ! CCF      convective cloud fraction
   !          - Inputs
   ! TP      temperature at (t+dt) before condensation
   ! TM      temperature at (t-dt)
   ! QP      specific humidity at (t+dt) before condensation
   ! QM      specific humidity at (t-dt)
   ! CWP     cloud water content at (t+dt) before condensation
   ! CWM     cloud water content at (t-dt)
   ! PSP     surface pressure at (t+dt)
   ! PSM     surface pressure at (t-dt)
   ! S       sigma level values
   ! PRFLX   flux of liquid precipitation
   ! SWFLX   flux of solid precipitation
   ! FEVP    evaporation of precipitation
   ! F12     cloud to rainwater collection tendency
   ! ICEFRAC  ice fraction
   ! MRK2    Markov chains for stochastic perturbations
   ! ilab     label array from convective scheme
   ! ni      number of grid points in the horizontal
   ! nlev    number of levels
   ! TAU     timestep
   !@Note
   !-----------------------------------------------------------------------
   !  I)       NAMES OF PARAMETERS AND OTHER QUANTITIES
   !           ------------------------------------------------------------
   !       CBFEFF      INCREASES CONVERSION RATE DUE TO DUE TO PRESENCE
   !                   OF ICE IN PRECIPITATION COMING IN FROM ABOVE
   !       CFREEZ      INCREASES CONVERSION RATE BELOW TEMP CTFRZ1
   !       COALES      INCREASES CONVERSION RATE DUE TO PRECIPITATION
   !                   COMING IN FROM ABOVE
   !       CONAE       FACTOR TO TUNE TABLE LOOK-UP FOR TETA-AE
   !       CTFRZ1      TEMP BELOW WHICH CONVERSION RATE IS INCREASED
   !       cumask      SET = 0 IF CONVECTION,   OTHERWISE = 1
   !       DPRG        DELTA-P DIVIDED BY GRAVITY
   !       HDCWAD      TENDENCY OF CLOUD WATER DUE TO EFFECTS OTHER
   !                   THAN CONDENSATION
   !       HDQAD       TENDENCY OF VAPOUR DUE TO EFFECTS OTHER THAN
   !                   CONDENSATION
   !       HDTAD       TENDENCY OF TEMPERATURE DUE TO EFFECTS OTHER THAN
   !                   CONDENSATION
   !       HDPMX       MAXIMUM PRECIPITATION CHANGE DUE TO EVAPORATION
   !       HPK         = P AT SIGMA
   !       HQSAT       SATURATION SPECIFIC HUMIDITY: Qs
   !       HSQ         dQs/dT
   !       HSQ2        dQs/dlnP
   !       HU          RELATIVE HUMIDITY OF THE ENVIRONMENT
   !       PRCPCU      RATE OF CONVECTIVE PRECIPITATION AT LEVEL K
   !       PRCPST      RATE OF STRATIFORM PRECIPITATION AT LEVEL K

   !       HCCU        CONVERSION RATE FROM CLOUD TO PRECIP DROPS IN
   !                   CONVECTIVE CLOUD
   !       HCST        CONVERSION RATE FROM CLOUD TO PRECIP DROPS IN
   !                   STRATIFORM CLOUD
   !       HE273       SATURATION VAPOUR PRESSURE AT T=273K
   !       HKMELT      COEFFICIENT FOR MELTING OF ICE
   !       HMRCU       CLOUD WATER MIXING RATIO AT WHICH CONVERSION BECOMES
   !                   EFFICIENT IN CONVECTIVE CLOUD
   !       HMRST       CLOUD WATER MIXING RATIO AT WHICH CONVERSION BECOMES
   !                   EFFICIENT IN STRATIFORM CLOUD (NAMELIST-CONTROLLED)
   !       HPS         TIME AVERAGED SURFACE PRESSURE
   !       HUZ00       MODIFIED HU00
   !       HU00        THRESHOLD RELATIVE HUMIDITY FOR STRATIFORM
   !                   CONDENSATION
   !       HU0MAX      MAXIMUN ALLOWABLE VALUE OF MODIFIED HU00
   !       HU0MIN      MINIMUN ALLOWABLE VALUE OF MODIFIED HU00
   !       STPEVP      EVAPORATION COEFFICIENT FOR STRATIFORM PRECIPITATION
   !       TABDE       DIFFERENCE IN SATURATION VAPOUR PRESSURE OVER
   !                   WATER AND ICE
   !       TABFBF      BERGERON-FINDEISEN EFFECT FROM (DEWI*TABICE)
   !       TABICE      PROBABILITY FOR ICE CRYSTALS AS A FUNCTION OF
   !                   TEMPERATURE
   !-----------------------------------------------------------------------
   !  II)      DECLARATIONS
   !           ------------------------------------------------------------

   real    XCOND  , XN     , HDPAD  , HSQ2   , HDQAD  , &
        HDQSAD , XDE    , XPRB   , BFMOD  , &
        XK     , HFCOX  , XFT    , HFREZX , HFRCOA , HFMRX  , &
        ZCWP   , XPRADD , DTMELT , DMELT  , EVAPRI , &
        XEVACU , XP     , QINCR  , HP0    , HE273  , HEDR   , &
        HDLDCP , HELDR  , HEDLDR , CONAE  , AECON  , CFREEZ , &
        COALES , SIGMIN
   real    CBFEFF , CTFRZ1 , HCCU   , HMRCU , &
        HKMELT , XDT    , DSNMAX , XSNOW  , &
        rTAU   , SNOW   , PRCP   , CONET  , &
        COEF   , COVER  , HMR    , ZDCW   , XT     , &
        SIGMAX , T0I    , WEIGHT , x      , y      , z      , &
        TCI    , TSCALE , APRI   , TOPEQ0 , TODPMX , temp1  , &
        temp2  , XB     , XBB    , xo

   integer il     , jk,inr
   real xxp_t,xhj_t,xf_t,xfprim_t,hsq,huz00t,hu,hcondt
   real xwrk,HACCES

   integer, dimension(NI,NLEV) :: CUMASK
   real, dimension(NI     ) :: HCST
   real, dimension(NI     ) :: STPEVP
   real, dimension(NI     ) :: HPS
   real, dimension(NI     ) :: COVBAR
   real, dimension(NI     ) :: PRCPST
   real, dimension(NI     ) :: STSNOW
   real, dimension(NI     ) :: PRCPCU
   real, dimension(NI     ) :: CUSNOW
   real, dimension(NI     ) :: HSCT
   real, dimension(NI,NLEV) :: SUBCLD
   real, dimension(NI     ) :: HDQMX
   real, dimension(NI     ) :: HDPMX
   real, dimension(NI     ) :: HCOND
   real, dimension(NI,NLEV) :: HQSAT
   real, dimension(NI,NLEV) :: HQSATP
   real, dimension(NI     ) :: HLDCP
   real, dimension(NI     ) :: HDCWAD
   real, dimension(NI,NLEV) :: DPRG
   real, dimension(NI     ) :: HCIMP
   real, dimension(NI     ) :: HDTAD
   real, dimension(NI,NLEV) :: HPK
   real, dimension(NI,NLEV) :: PRESP
   real, dimension(NI,NLEV) :: PRESM
   real, dimension(NI     ) :: XPRBT
   real, dimension(NI     ) :: XDET
   real, dimension(NI     ) :: XDET1
   real, dimension(NI     ) :: CONETT
   real, dimension(NI     ) :: PRMODT
   real, dimension(NI     ) :: YMT
   real, dimension(NI     ) :: HBMRXT
   real, dimension(NI     ) :: COEFT
   real, dimension(NI     ) :: XFIXT
   real, dimension(NI     ) :: HU0MAX
   real, dimension(NI     ) :: HU0MIN
   real, dimension(NI     ) :: XBHU
   real, dimension(NI     ) :: HMRST
   real, dimension(NI     ) :: ICEACC

   !***********************************************************************
   !-----------------------------------------------------------------------
   !  III)     STATEMENT FUNCTIONS
   !           -----------------------------------------------------------

   real    Z1, Z2, Z3, Z4, Z5

   !           TEMPERATURE FUNCTION TO MULTIPLY HMRCU AND HMRST FOR T<273

   Z1(XT) = min(1.33*exp(-(min(0.,(XT-TRPL))*.066)**2) , 1.0)

   Z2(XT) = abs (XT - 232.) / 18.

   Z3(XT) = Z2(XT) * (1. + Z2(XT) * (1. + 1.333 * Z2(XT)))

   Z4(XT) = Z3(XT) / (1. + Z3(XT)) * sign(1.0,XT-232.)

   Z5(XT) = max(0.5*0.15*(1.07+Z4(XT)),0.03)

   !-----------------------------------------------------------------------
   !  IV)      VALUES OF CONSTANTS - IN SI UNITS - INCLUDING DERIVED ONES
   !           -----------------------------------------------------------

   HE273  = 610.78
   T0I    = 1./TRPL
   TCI    = 232.
   TOPEQ0 = 268.
   TODPMX = 256.
   TSCALE = (TODPMX - TCI)*sqrt(2.)
   APRI = 1./(1.-exp(-((TOPEQ0-TCI)/TSCALE)**2))

   HP0    = 1.E5
   HEDR = EPS1/RGASD
   HEDLDR = EPS1*CHLF/RGASD
   HDLDCP = CHLF/CPD

   !-----------------------------------------------------------------------
   !  V)       PARAMATER VALUES IN SI UNITS
   !           ------------------------------------------------------------

   CONAE  = 0.15
   AECON  = exp(CONAE)
   rTAU = 1. / TAU

   !-----------------------------------------------------------------------
   !  VI)      PREPARATIONS
   !           ------------------------------------------------------------

   do il = 1, ni
      HPS(il) = ( PSP(il) + PSM(il) ) * 0.5
   end do

   do jk = 1, nlev
      do il = 1, ni
         HPK(il,jk) = S(il,jk) * HPS(il)
      end do
   end do

   do il = 1, ni
      DPRG(il,1) = 0.5 * ( HPK(il,2) - HPK(il,1) ) / GRAV
      DPRG(il,nlev) = ( 0.5 * ( HPK(il,nlev) - HPK(il,nlev-1) ) &
           + S(il,nlev+1) * HPS(il) - HPK(il,nlev) ) / GRAV
   end do

   do jk  = 2, nlev-1
      do il  = 1, ni
         DPRG(il,jk) = 0.5 * ( HPK(il,jk+1) - HPK(il,jk-1) ) / GRAV
      end do
   end do

   !           SUNQVIST stratiform condensation scheme
   !           ------------------------------------------------------------


   !  A)       LINKAGE CONDITIONS

   !           FROM KUO ilab DEFINE cumask
   !           FROM cumask CALCULATE subcld

   where (ilab == 2)
      cumask = 0
   elsewhere
      cumask = 1
   end where

   do il =1, ni
      subcld(il,1) = 0.
   end do

!VDIR NOLSTVAL
   do jk = 2, nlev
      do il = 1, ni
         subcld(il,jk) = subcld(il,jk-1)
         if( cumask(il,jk).eq.0 .and. HPK(il,jk).le.40000. ) &
              subcld(il,1) = 1.
         !               cloud top above 400 mb
         !               to indicate a likely downdraft below 500 mb
         if( subcld(il,1).gt.0. .and. HPK(il,jk).gt.50000. ) &
              subcld(il,jk)=min(1.,max(0.,2.5*(s(il,jk)-0.5)))
         if( cumask(il,jk).eq.1 .and. cumask(il,jk-1).eq.0) &
              subcld(il,jk) = 1.
      end do
   end do
   do il =1, ni
      subcld(il,1) = 0.
   end do


   !  B)       PARAMATER VALUES IN SI UNITS

   HU0MIN(:) = ens_spp_get('hu0min', mrk2, default=cond_hu0min)
   HU0MAX(:) = ens_spp_get('hu0max', mrk2, default=cond_hu0max)
   SIGMIN = 0.7
   SIGMAX = 0.9
   XBHU(:) = ( HU0MAX(:) - HU0MIN(:) ) / ( SIGMAX - SIGMIN )
   xo = 1.e-12

   !***********************************************************************

   !           HERE BEGINS THE CALCULATION OF STRATIFORM CONDENSATION

   !***********************************************************************

   !           SUNQVIST cloud water and precipitation scheme
   !           ------------------------------------------------------------


   !  A)       LINKAGE CONDITIONS


   !-----------------------------------------------------------------------

   !  B)       PARAMATER VALUES IN SI UNITS

   CFREEZ = 0.12
   COALES = 300.
   CBFEFF = 4.0
   CTFRZ1 = 263.
   HCCU   = 1.E-4
   HMRCU  = 5.E-4
   HKMELT = 3.E-5
   xo = 1.E-16
   HCST(:)   = ens_spp_get('cond_hcst', mrk2, default=1.E-4)
   STPEVP(:) = 2.*GRAV * ens_spp_get('cond_evap', mrk2, default=cond_evap)
   HMRST(:) = ens_spp_get('cond_hmrst', mrk2, default=cond_hmrst)
   ICEACC(:) = ens_spp_get('cond_iceacc', mrk2, default=cond_iceacc)

   !-----------------------------------------------------------------------

   !  C)       INITIALIZATIONS

   do il = 1, ni
      PRCPST(il) = 0.
      STSNOW(il) = 0.
      PRCPCU(il) = 0.
      CUSNOW(il) = 0.
      COVBAR(il) = 0.
      HSCT(il) = 0.
   end do
   PRFLX = 0.
   SWFLX = 0.
   fevp = 0.
   f12 = 0.

   do jk = 1, nlev
      do il = 1, ni
         HSCT(il) = HSCT(il) + amin1( 0. , CTT(il,jk) ) * DPRG(il,jk)
      end do
   end do
   do jk = 1, nlev
      do il = 1, ni
         PRESP(il,jk)=S(il,jk)*PSP(il)
         PRESM(il,jk)=S(il,jk)*PSM(il)
         HQSAT(il,jk) = foqst(TM(il,jk),PRESM(il,jk))
         HQSATP(il,jk) = foqst(TP(il,jk),PRESP(il,jk))
         HLDCP(il) = exp(-(((max(TM(il,jk),tci)-tci)/tscale)**2))

         xwrk= max(((apri*(HLDCP(il)-1.0))+1.0),0.0)
         HLDCP(il)=(CHLC + (CHLF * xwrk))/CPD
         xprbt(il)=xwrk
         temp1 = HQSAT(il,jk)
         temp2 = TM(il,jk)
         HSQ= FODQS( temp1 , temp2 )

         HSQ2 = - HQSAT(il,jk) * ( 1. + DELTA * HQSAT(il,jk) )
         HU= QP(il,jk)/HQSATP(il,jk)
         HU = amin1( HU, 1. )
         HU = amax1( HU, 0. )

         HUZ00t= HU0MIN(il) + XBHU(il) * ( S(il,jk) - SIGMIN )
         HUZ00t= max( HU0MIN(il), min( HU0MAX(il), HUZ00t ) )

         x = 1.+ 0.15 * max( 0., 238. - TM(il,jk) )
         x = ( HU0MAX(il) - HUZ00t ) * ( 1. - 1. / x )
         HUZ00t= AMIN1( HUZ00t + x , HU0MAX(il) )

         SCF(il,jk) = 1. - sqrt( (1.-HU) / (1.-HUZ00t) )
         SCF(il,jk) = amax1( SCF(il,jk)-CCF(il,jk) , 0. )
!!$         mask = sign(1., -1.*abs(SCF(il,jk)))+1./2.
!!$         HCONDt = mask * (- CWP(il,jk) * rTAU * cumask(il,jk))
         if( SCF(il,jk) .eq. 0. ) then
            HCONDt = - CWP(il,jk) * rTAU * cumask(il,jk)
         else
            HCONDt = 0.
         endif
         XB = SCF(il,jk)
         XN = 2. * HQSAT(il,jk) * (1.-HUZ00t) * XB * ( 1.-XB )
         XK =   ( XB * ( XN - xo ) + xo ) &
              / ( XB * ( XN + CWM(il,jk) ) + xo )

         QM(il,jk) = amin1( QM(il,jk) , HQSAT(il,jk) )

         HDTAD(il) = ( TP(il,jk) - TM(il,jk) ) * rTAU
         HDQAD= ( QP(il,jk) - QM(il,jk) ) * rTAU
         HDCWAD(il) = ( CWP(il,jk) - CWM(il,jk) ) * rTAU

         HDPAD = ( PSP(il) - PSM(il) ) / HPS(il) * rTAU
         HDQSAD = HSQ * HDTAD(il) + HSQ2 * HDPAD
         HDQSAD = max( HDQSAD, - HQSAT(il,jk) * rTAU )
         HCIMP(il) = 1. / ( 1. + HLDCP(il) * HSQ )
         HACCES= ( HDQAD- HU* HDQSAD ) &
              / ( 1. + HU * HLDCP(il)*HSQ )
         HDQMX(il) = ( ( QM(il,jk) - HQSAT(il,jk) ) * rTAU &
              +  HDQAD- HDQSAD ) * HCIMP(il)

         !           ( n.b.important limits: if b=0, then k=1; if b=1, then k=0 )

         XCOND =   ( 1. - XK * ( 1. - XB ) ) * HACCES
         XCOND = amax1( XCOND , - CWP(il,jk) * rTAU )

         HCONDt = HCONDt + XCOND


         !           ------------------------------------------------------------
         !           CORRECT THE above CALCULATIONS OF NET CONDENSATION
         !           a) IN CASES OF RESIDUAL SUPER-SATURATION (because, in cloudy
         !           cases, moistening/cond. may have been over/under-estimated
         !           or super-saturation was present initially)
         !           b) FOR diagnosed (b=0) CLEAR SITUATION LEADING eventually
         !           TO SUPER-SATURATION
         !           ------------------------------------------------------------

         QINCR = HDQMX(il) - HCONDt

         XCOND = amax1( QINCR, 0.0 ) * cumask(il,jk)

         HCOND(il) = HCONDt + XCOND

         HDPMX(il) = amax1 ( 0. , - QINCR ) * DPRG(il,jk)

      end do

      !-----------------------------------------------------------------------

      !  D)       CALCULATIONS

      do il = 1, ni
         HELDR = HEDR * CPD * HLDCP(il)
         xdet(il) = exp(HELDR*(T0I - 1. / TM(il,jk)))
         xdet1(il) = exp(HEDLDR*(T0I - 1. / TM(il,jk)))

         CONET = amax1( 0. , CTT(il,jk) + HSCT(il) / DPRG(il,jk) )
         HSCT(il) = HSCT(il)+(amax1(0.,CTT(il,jk))-CONET)*DPRG(il,jk)
         CONETt(il) = HCOND(il) + CONET/HLDCP(il)
         PRCP = PRCPST(il)+PRCPCU(il)
         SNOW = STSNOW(il)+CUSNOW(il)
         COVER = SCF(il,jk) + CCF(il,jk) + 1.E-2
         WEIGHT = ( SCF(il,jk) + 1.E-2 ) / COVER
         HMR = HMRST(il) * WEIGHT + HMRCU * (1-WEIGHT)
         COEF = HCST(il) * WEIGHT + HCCU * (1-WEIGHT)

         !           ------------------------------------------------------------
         !           Factors for coalescence HFCOX, freezing HFREZX.
         !           Reduction of HMR at low temperatures, HFMRX
         !           Modified probability of ICE
         !           resulting from ICE in PRECIP from above
         !           ------------------------------------------------------------

         temp1 = TM(il,jk)
         xde = max(0.0,min(((HE273/temp1*xdet(il)* &
              (1. - xdet1(il)))*9.248487), 1.0))
         XPRB = xprbt(il)
         if (temp1 > 250.) then
            XFT = Z1(temp1)
         else
            XFT = Z5(temp1)
         endif

         PRMODt(il) = XPRB + ( 1. - XPRB ) * SNOW / ( PRCP + 1.E-7 )
         BFMOD = PRMODt(il) * ( 1. - XPRB ) * XDE

         HFCOX = 1. + COALES * sqrt( PRCP )

         HFREZX = 1. + CBFEFF * BFMOD
         HFREZX = HFREZX * ( 1. + CFREEZ * ( 1. - XFT ) / XFT )

         HFRCOA = HFCOX * HFREZX
         HFMRX = HMR * XFT / HFCOX
         HBMRXt(il) = COVER * HFMRX

         !           ------------------------------------------------------------
         !           Special treatment for T.LT.236K
         !           ------------------------------------------------------------

         temp1 = amax1(0.,amin1(1.,0.25*(TM(il,jk)-232.)))
         HFRCOA = temp1 * HFRCOA + ( 1.- temp1 ) * ICEACC(il)
         !           ------------------------------------------------------------
         !           Fixed part of the equation normalized by 2.*b*Mr
         !           ------------------------------------------------------------

         XFIXt(il) = ( 2. * CWM(il,jk) + TAU * &
              (HDCWAD(il) + CONETt(il)))/( 2. * HBMRXt(il) )


         !           ------------------------------------------------------------
         !           Conversion rate times 2*dt
         !           ------------------------------------------------------------

         COEFt(il) = 0.5 * COEF * HFRCOA * TAU

         !           ------------------------------------------------------------
         !           First guess YM is M(t-dt) normalized by b*Mr
         !           ------------------------------------------------------------

         YMt(il) = CWM(il,jk) / HBMRXt(il)

         !           ------------------------------------------------------------
         !           To make M(t+dt).ge.0, YM has to be .ge. M(t-dt)/(2*b*Mr)
         !           ------------------------------------------------------------

         xdet1(il)=0.5 * YMt(il)
      enddo

      !           ------------------------------------------------------------
      !           5 NEWTON - RAPHSON ITERATIONS
      !           ------------------------------------------------------------
      do inr=1,5
         do il=1,ni
            xdet(il) = exp(-min(ymt(il)*ymt(il), 25.0))
            xxp_t = xdet(il)
            xhj_t=1 + coeft(il)*(1-xxp_t)
            xf_t = xhj_t*ymt(il) -xfixt(il)
            xfprim_t = xhj_t + 2.*coeft(il)*ymt(il)*ymt(il)*xxp_t
            ymt(il)  = amax1 ( ymt(il) - xf_t/xfprim_t,xdet1(il))
         enddo
      enddo

      !           ------------------------------------------------------------
      !           Rate of change of cloud water content
      !           Generation of precipitation
      !           ------------------------------------------------------------

      do il = 1, ni
         xdet(il) = exp(-(((max(TP(il,jk), tci) - tci)/tscale)**2))
         ZCWP = amax1( 2. * HBMRXt(il) * YMt(il) - CWM(il,jk) , 0. )

         ZDCW = ( ZCWP - CWP(il,jk) ) * rTAU

         XPRADD = DPRG(il,jk) * amax1( CONETt(il) - ZDCW , 0. )

         !           we make sure that no infinitesimal precipitation is generated
         if (abs(conett(il)-zdcw).le.abs(spacing(zdcw))) xpradd = 0.

         SWT(il,jk) = ZDCW * cumask(il,jk)

         !           Diagnostics for AURAMS

         F12(il,jk) = amax1( CONETt(il) - ZDCW , 0. )
         if (ZCWP.lt.1.0e-09) F12(il,jk)=0.0
         ICEFRAC(il,jk) = max(((apri*(xdet(il)-1.0))+1.0),0.0)


         PRCPST(il) = PRCPST(il) +       XPRADD * cumask(il,jk)
         STSNOW(il) = STSNOW(il) + PRMODt(il)*XPRADD * cumask(il,jk)

         CWT(il,jk) = ZDCW * (1-cumask(il,jk))
         PRCPCU(il) = PRCPCU(il) +       XPRADD * (1-cumask(il,jk))
         CUSNOW(il) = CUSNOW(il) + PRMODt(il)*XPRADD * (1-cumask(il,jk))


         !           ------------------------------------------------------------
         !           Melting of stratiform snow
         !           ------------------------------------------------------------


         XB = SCF(il,jk)
         XBB = COVBAR(il)
         XDT = TM(il,jk) - TRPL &
              + TAU * ( HDTAD(il) + HLDCP(il) * CONETt(il) )
         XDT = amax1( XDT , 0. )
         DSNMAX = XDT * DPRG(il,jk) / ( TAU * HDLDCP )
         XN = HKMELT * ( TAU * HDLDCP )
         z = sqrt ( amax1( xo , XBB ) )

         SNOW = STSNOW(il)
         x = amax1( SNOW , 1.E-16 )
         x = sqrt( x )
         y = 0.5 * XN * DSNMAX / x
         y = XBB * y / ( z + 0.5 * XN * x )
         y = amin1( y , 1. )
         XSNOW = SNOW * ( 1. - y ) ** 2

         DMELT = amin1( SNOW - XSNOW , XBB * DSNMAX )
         DTMELT = HDLDCP * DMELT / DPRG(il,jk)

         STSNOW(il) = amax1( 0. , STSNOW(il) - DMELT )

         !           ------------------------------------------------------------
         !           Melting of convective snow
         !           ------------------------------------------------------------

         SNOW = CUSNOW(il)
         x = amax1( SNOW , 1.E-16 )
         x = sqrt( x )
         y = 0.5 * XN * DSNMAX / x
         y = y / ( 1. + 0.5 * XN * x )
         y = amin1( y , 1. )
         XSNOW = SNOW * ( 1. - y ) ** 2

         DMELT = amin1( SNOW - XSNOW , DSNMAX )

         CUSNOW(il) = amax1( 0. , CUSNOW(il) - DMELT )


         !           ------------------------------------------------------------
         !           Evaporation of stratiform precipitation
         !           ------------------------------------------------------------


         PRCP = PRCPST(il)
         XN = STPEVP(il) * TAU / HCIMP(il)
         x = amax1( PRCP , 1.E-16 )
         x = sqrt ( x )
         y = 0.5 * XN * HDPMX(il) / x
         y = XBB * y / ( z + 0.5 * XN * x )
         y = amin1( y , 1. )
         XP = PRCP * ( XB + ( 1. - XB ) * ( 1. -  y ) ** 2 )

         EVAPRI = amin1( PRCP - XP ,  XBB * HDPMX(il) )

         if (PRCPST(il).gt.0.) FEVP(il,jk) = EVAPRI/(PRCPST(il)+1.0E-28)

         PRCPST(il) = amax1( 0. , PRCPST(il) - EVAPRI )
         STSNOW(il) = amax1( 0. , STSNOW(il) - EVAPRI )
         HCOND(il) = HCOND(il) - EVAPRI / DPRG(il,jk)

         COVBAR(il) = XBB * ( 1. - XB ) + XB

         !           ------------------------------------------------------------
         !           Evaporation of precipitation under the convective cloud
         !           ------------------------------------------------------------

         EVAPRI =  - DBDT(il) * TAU * HDQMX(il) * DPRG(il,jk)
         EVAPRI = EVAPRI * subcld(il,jk)
         EVAPRI = amin1( PRCPCU(il) , amax1( 0. , EVAPRI ) )
         XEVACU = EVAPRI / DPRG(il,jk)

         PRCPCU(il) = PRCPCU(il) - EVAPRI
         CUSNOW(il) = amax1( 0., (CUSNOW(il) - EVAPRI) )

         !           ------------------------------------------------------------
         !           TEMPERATURE AND MOISTURE TENDENCIES, PRECIPITATION FLUXES
         !           ------------------------------------------------------------

         STT(il,jk) = - DTMELT + HCOND(il) * HLDCP(il)
         SQT(il,jk) = - HCOND(il)

         CTT(il,jk) = CTT(il,jk) - XEVACU * HLDCP(il)
         CQT(il,jk) = CQT(il,jk) + XEVACU

         PRFLX(il,jk+1) =  PRCPST(il) + PRCPCU(il)
         SWFLX(il,jk+1) =  STSNOW(il) + CUSNOW(il)

      end do
   end do

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !  9)       SAVE THE STRATIFORM AND CONVECTIVE PRECIPITATION RATES
   !           AND THE LIQUID AND SOLID PRECIPITATION FLUXES
   !           ------------------------------------------------------------

   do il = 1, ni
      SRR(il) = PRCPST(il) - STSNOW(il)
      SSR(il) = STSNOW(il)
      CRR(il) = PRCPCU(il) - CUSNOW(il)
      CSR(il) = CUSNOW(il)
   end do

   do jk = 1, nlev+1
      do il = 1, ni
         PRFLX(il,jk) =  PRFLX(il,jk) - SWFLX(il,jk)
      end do
   end do
   return

   !-----------------------------------------------------------------------

   !***********************************************************************

   !           HERE ENDS THE CALCULATION OF CLOUD WATER CONTENT
   !                         AND PRECIPITATION

   !***********************************************************************

end subroutine CONSUN3
