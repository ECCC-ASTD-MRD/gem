!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

subroutine skocon5(ZTE   , ZQE   , ZCWE, ZSRR, &
     ZSSR  , ZSTCOV, TP    , TM    , QP  , &
     QM    , CWP    , CWM   , PSP   , &
     PSM   , ILAB  , S      , NI    , NLEV  , &
     DT    , SATUCO,  CONVEC, prflx , swflx)
   use tdpack
   implicit none
!!!#include <arch_specific.hf>

   integer ni , nlev
   integer ilab(ni,nlev)
   real    ZTE(NI,nlev)    , ZQE(NI,nlev), ZCWE(NI,nlev), &
        ZSRR(NI)    , ZSSR(NI)  , &
        ZSTCOV(NI,nlev), &
        TP(NI,nlev)   , TM(NI,nlev) , QP(NI,nlev) ,  QM(NI,nlev), &
        CWP(NI,nlev), CWM(NI,nlev), &
        PSP(NI)       , PSM(NI)     , S(ni,*)     , &
        PRFLX(ni,nlev+1), SWFLX(ni,nlev+1)
   real    DT
   character(len=*) :: CONVEC
   logical SATUCO

   !@Author Ashu Dastoor (1991)

   !@Revision
   ! 001     Wei Yu (January 1994) - Code vectorization
   ! 002     Wei Yu (June 1994) - Adaptation to MC2, use local sigma
   ! 003     Bernard Bilodeau - New physics interface
   ! 004     G. Pellerin (June 95) - Different Evaporation Coefficient
   !                                 for convection
   !                               - Bugfix for stable and convective
   !                                 precipitation
   ! 005     R. Sarrazin (June 95) modified stratiform part
   ! 006     R. Sarrazin (Nov  95) put very small precip rates to 0 at the end
   !                               to avoid excessive trace of precip
   ! 007     B. Bilodeau (Jan 01) - Automatic arrays
   ! 008     L. Spacek   (May 03) - IBM conversion
   !             - eliminate NEWKUO code
   !             - eliminate useless calculations
   !             - calls to vsexp routine (from massvp4 library)
   ! 009     R. McTaggart-Cowan (Jul 2006) - eliminate function FMROFT
   ! 010     L. Spacek   (Now 11) - copy prflx,swflx into the lowest level


   !@Object parameterization of stratiform and convective condensation

   !@Arguments
   ! ZTE     tendency of temperature due to convection
   ! ZQE     tendency of specific humidity due to convection
   ! ZCWE    tendency of cloud water due to convection
   ! ZSRR    rate of liquid precipitation due to stratiform condensation
   ! ZSSR    rate of solid precipitation due to stratiform condensation
   ! ZSTCOV  stratiform cloud cover
   ! TP      temperature
   ! TM      environmental T at (T-DT)
   ! QP      specific humidity
   ! QM      environmental Q at (T-DT)
   ! CWP     cloud water content at (T+DT)
   ! CWM     cloud water content at (T-DT)
   ! PSP     surface pressure at (T+DT)
   ! PSM     surface pressure at time (T-DT)
   ! ILAB    flag array from subroutine KUO
   ! S       sigma levels
   ! NI      1st horizontal dimension
   ! NLEV    number of model levels
   ! DT      timestep
   ! SATUCO  parameter to indicate
   ! CONVEC  convection switch
   ! PRFLX   flux of liquid precipitation
   ! SWFLX   flux of solid precipitation

   !@Note

   !  HDCWAD modified, DO NOT PASS to the main routine (1990)


   !  This routine deals with parameterization of condensation along
   !  a constant J - LINE
   !  Points 1) - 8) Handle parameterization of convection by use of
   !  a modified KUO scheme such that it includes prediction
   !  of cloud water content and diagnostic setting of fractional
   !  cloud cover.
   !  Points 11) - 17) Handle parameterization of stratiform condensa-
   !  tion including prediction of cloud water content and
   !  diagnostic setting of fractional cloud cover.

   !  1)   Find lifting condensation level, LCL, HPB of surface air
   !       and corresponding temperature, TB and Humidity, QB. Calculate
   !       THETAE of HPB
   !  2)   Calculate TC(K) and QC(K) along THETAE = CONSTANT
   !  3)   Calculate (TC-T) and (TC-T)-MAX and (QC-Q). Set lowest level
   !       of convection, KB, and top level, KT, and level of (TC-T)-MAX,
   !       KTMAX.  Calculate (HPB - P(KB)).
   !  4)   Determine top of the convection. If the result is KT.GE.KB,
   !       then further calculations are terminated
   !  5)   Calculate accession of vapour, CQ, and basic proportionality
   !       parameter, KSIO.  If CQ.LE.O  or KSIO.LE.O, KT is set = KB
   !       and jump to POINT 5)
   !  6)   Calculate KSI(K), heating and moistening functions and cloud
   !       cover.
   !  7)   Prediction of cloud water content, CW, by an implicit scheme,
   !       requiring an iterative procedure.
   !  8)   Calculate accumulated precipitation.
   !-----------------------------------------------------------------------

   ! 11)   Modification of basic threshold relative humidity, HU00
   ! 12)   General requirement for condensation, REQCON = 1 or 0
   ! 13)   Setting of specific HU00 and REQCON and EVAP of cloud water
   ! 14)   Cloud cover and evaporation of precipitation
   ! 15)   Accession of vapour, and final setting of heating/cooling
   ! 16)   Prediction of cloud water content, CW, by an implicit scheme,
   !       requiring an iterative procedure.
   ! 17)   Calculate accumulated precipitation.

   !*
   !-----------------------------------------------------------------------
   !  I)       NAMES OF PARAMETERS AND OTHER QUANTITIES
   !           ------------------------------------------------------------

   !       TABFBF      BERGERON-FINDEISEN EFFECT FROM (DEWI*TABICE)
   !       CUMASK(K)   SET = 0 IF CONVECTION,   OTHERWISE = 1
   !       DLNPDT(K)   = 1/HPK*DP/DT AT LEVEL K
   !       DPK(K)      DELTA-P AROUND LEVEL K
   !       FSCFLD(I,J) FIELD TO MODIFY SC COVER AT LEVEL KSCFLD
   !                   DUE TO ENTRAINMENT
   !       HDCWAD(K)   TENDENCY OF CLOUD WATER DUE TO OTHER EFFECTS
   !                   BUT CONDENSATION
   !       HDCWST(K)   TENDENCY OF CLOUD WATER DUE TO STRATIFORM
   !       HDQAD(K)    TENDENCY OF VAPOUR DUE TO EFFECTS OTHER THAN
   !                   CONDENSATION
   !       HDQST(K)    TENDENCY OF VAPOUR DUE TO STRATIFORM CONDENSATION
   !       HDTAD(K)    TENDENCY OF TEMPERATURE DUE TO EFFECTS OTHER THAN
   !                   CONDENSATION
   !       HDTST(K)    TENDENCY OF TEMPERATURE DUE TO STRATIFORM
   !                   CONDENSATION
   !       TABDE       DIFFERENCE IN SATURATION VAPOUR PRESSURE OVER
   !                   WATER AND ICE
   !       HEVAC(K)    EVAPORATION RATE OF CLOUD WATER
   !       HEVAPR(K)   EVAPORATION RATE OF PRECIPITATION
   !       HKSIZ(K)    KSI OF KUO SCHEME AS A FUNCTION OF HEIGHT
   !       FMROFT      TEMPERATURE FUNCTION TO MULTIPLY MR WITH
   !                   FOR T<273
   !       HPK(K)      = P AT SIGMA(K)
   !       HPKAP(K)    = (HP0/HPK(JK)) ** KAPPA
   !       HPLOGK(K)   = LOG(P)
   !       HQC(K)      SATURATION MIXING RATIO ALONG THE MOIST ADIABAT
   !                   THROUGH THE LCL
   !       HQCMQ(K)    HQC MINUS ENVIRONMENTAL Q
   !       QM(K)       ENVIRONMENTAL Q AT (T-DT)
   !       HQP1(K)     ENVIRONMENTAL Q AT (T+DT)
   !       HQSAT(K)    SATURATION MIXING RATIO WITH RESECT TO
   !                   ENVIRONMENTAL TEMPERATURE
   !       HSQ(K)      EPS*L**2*QSAT/(R*CP*T**2)
   !       ZSTCOV(K)   CONTAINS STRATIFORM CLOUD COVER
   !       HTC(K)      TEMPERATURE ALONG THE MOIST ADIABAT
   !                   THROUGH THE LCL
   !       HTCMT(K)    HTC MINUS ENVIRONMENTAL T
   !       TM(K)       ENVIRONMENTAL T AT (T-DT)
   !       HTP1(K)     ENVIRONMENTAL T AT (T+DT)
   !       HU(K)       RELATIVE HUMIDITY OF THE ENVIRONMENT
   !       TABICE      PROBABILITY FOR ICE CRYSTALS AS A FUNCTION OF
   !                   TEMPERATURE
   !       PRCPCU(I,J) RATE OF CONVECTIVE PRECIPITATION AT LEVEL K, SUMMED
   !                   FROM KHT
   !       PRCPST(I,J) RATE OF STRATIFORM PRECIPITATION AT LEVEL K, SUMMED
   !                   FROM KHT
   !       CCWM        WORKING COPY OF CWM
   !       ACPRST      ACCUMULATED PRECIPITATION FROM STRATIFORM CLOUDS
   !       ACPRCU      ACCUMULATED PRECIPITATION FROM CONVECTIVE CLOUDS
   !       CFREEZ      INCREASES CONVERSION RATE BELOW TEMP CTFRZ1
   !       COALES      INCREASES CONVERSION RATE DUE TO PRECIPITATION
   !                   COMING IN FROM ABOVE
   !       CBFEFF      INCREASES CONVERSION RATE DUE TO DUE TO PRESENCE
   !                   OF ICE IN PRECIPITATION COMING IN FROM ABOVE
   !       CONAE       FACTOR TO TUNE TABLE LOOK-UP FOR TETA-AE
   !       CTFRZ1      TEMP BELOW WHICH CONVERSION RATE IS INCREASED
   !       DT          TIME-STEP
   !       HCP         SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
   !       HCCU        CONVERSION RATE OF CLOUD TO PRECIP DROPS IN
   !                   CONVECTIVE CLOUD
   !       HCST        CONVERSION RATE OF CLOUD TO PRECIP DROPS IN
   !                   STRATIFORM CLOUD
   !       HCUNRM      CLOUD COVER DEPENDS ON CLOUD DEPTH COMPARED TO
   !                   HCUNRM
   !       HDPB        = (HPB - HPK(KHB) + 0.5 * DPK(KHB))
   !       HELDCP      = HEPS * HLDCP
   !       HELDR       = HEPS * HLV/HRD
   !       HEPS        MOL WEIGHT OF WATER/MOL WEIGHT OF AIR = 0.622
   !       HEPSE0      = HEPS * HE273
   !       HE273       SATURATION VAPOUR PRESSURE AT T=273K
   !       HG          GRAVITY
   !       HKAP        = HRD/HCP
   !       HKMELT      COEFFICIENT FOR MELTING OF ICE IN PRECIPITATION
   !       HKPI        = 1./HKAP
   !       HLDCP       = HLV/HCP
   !       HLV         LATENT HEAT OF VAPORIZATION
   !       HMRCU       CLOUD WATER MIXING RATIO AT WHICH CONVERSION BECOMES
   !                   EFFICIENT IN CONVECTIVE CLOUD
   !       HMRST       CLOUD WATER MIXING RATIO AT WHICH CONVERSION BECOMES
   !                   EFFICIENT IN STRATIFORM CLOUD
   !       HPB         P VALUE AT LCL
   !       HPNORM      = HPS - HPT
   !       HPS         CURRENT SURFACE PRESSURE
   !       HPSKAP      (HP0/HPS) ** KAPPA
   !       HPT         TOP PRESSURE OF THE MODEL ATMOSPHERE (SIGMA-DOT=0)
   !       HP0         REFERENCE PRESSURE (=100 KPA)
   !       HQB         MIXING RATIO AT LCL
   !       HRD         GAS CONSTANT OF DRY AIR
   !       HRTDEL      = HRD * HT273/(HEPS * HLV)
   !       HTAUCU      CHARACTERISTIC TIME USED IN CONVECTIVE CLOUD COVER
   !                   SCHEME
   !       HTB         TEMPERATURE AT LCL
   !       HTD         DEW POINT TEMPERATURE OF SURFACE AIR
   !       HTETAE      THETAE THROUGH LCL
   !       HT273       = 273 K
   !       HUZ00       MODIFIED HU00
   !       HU00        THRESHOLD RELATIVE HUMIDITY FOR STRATIFORM
   !                   CONDENSATION
   !       HVSNOW      TERMINAL VELOCITY OF ICE/SNOW PRECIPITATION
   !       HVTERM      TERMINAL VELOCITY OF PRECIPITATION
   !       ICLDDT      OPTIONAL NUMBER OF STEPS BETWEEN CONDENSATION CALCU-
   !                   LATIONS; 1 FOR EVERY AND 2 FOR EVERY SECOND STEP ETC
   !                   THE VALUE IS SET IN ROUTINE STEP
   !       KCBTOP      NORMALLY = 1, BUT IF THETAE IS SMALL, KCBTOP IS SET
   !                   TO .GT. 1, IMPLYING LOWER ALTITUDE IN ORDER TO AVOID
   !                   TC .LT. TVBEG FOR THE VAPSAT TABLE
   !       KHB         FIRST MODEL LEVEL ABOVE LCL
   !       KHT         TOP LEVEL OF CONVECTION, I.E.,GOING UPWARD IT IS
   !                   THE LAST LEVEL WITH (TC-T).GT.0
   !       NLEV        NUMBER OF MODEL LEVELS.
   !       PCUACC      ACCUMULATED PRECIPITATION FROM CONVECTIVE CLOUDS
   !       PSTACC      ACCUMULATED PRECIPITATION FROM STRATIFORM CLOUDS
   !       REQCON      INDICATOR = 1 IF STRATIFORM CONDENSATION IS ALLOWED
   !                   TO TAKE PLACE,  ELSE REQCON = 0
   !       SEVAPR      ACCUMULATED  EVAPORATION FROM PRECIPITATION
   !       SEVAPC      ACCUMULATED EVAPORATION OF CLOUD WATER
   !       SUPSAT      INDICATOR = 1 IF SUPERSATURATION IS AT HAND, I.E.,
   !                   U.GT.1.02,  ELSE SUPSAT = 0
   !       STPEVP      COEFFICIENT EVAPORATION FOR STRATIFORM PRECIPITATION
   !       CVPEVP      COEFFICIENT EVAPORATION FOR CONVECTION PRECIPITATION
   !       TANVIL      IF T(KHT) .LE. TANVIL A STRATIFORM ANVIL CLOUD IS
   !                   CALCULATED IN THE TOP CONVECTIVE LAYER
   !       TCMTMX      (TC-T) - MAX
   !       DT          TIMESTEP
   !       RTDT        = 1. / DT
   !       U00MAX      MAXIMUN ALLOWABLE VALUE OF MODIFIED HU00
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !  II)      DECLARATIONS
   !           ------------------------------------------------------------


   real    hfis
   real    YMN    , XXP    , XHJ    , XF     , XFPRIM , &
        XCWP1  , XPRADD , DTMELT , XROV   , XRO    , &
        XMIND  , DMELT  , SQMELT , XMELT  , DSNOW  , &
        EVAPRI , XEVAP  , XPSQ   , xp
   real    DU238  , XCLEAR , &
        QPREL  , QINCR  , &
        XSCDIF , ICLDDT
   integer KSCFLD , KSCTOP , IT0    , IT0P1
   real    FSCFLD , FSCFKZ , HCP    , HEPS   , HG     , &
        rHG    , HLV    , HDL    , HP0    , HRD    , &
        HT273  , HE273  , HEPSE0 , HKAP   , HKPI   , &
        HCPI   , HEDR   , HRTDE  , HLDCP  , HDLDCP , &
        HRTDEL , HELDR  , HELDCP , HEDLDR , HELDRI , &
        SHCUM  , CONAE  , AECON  , CFREEZ , COALES
   real    CBFEFF , CTFRZ1 , HCCU   , PMOIST , HCUNRM , &
        HMRCU  , HMRST  , HTAUCU , TANVIL , HCST   , &
        HVTERM , CVPEVP , STPEVP , HVSNOW , SQVSNO , HKMELT , &
        U00MAX , HU00   , RDT , &
        ICLDDT2dt       , ICLDDTdt , &
        rICLDDT2dt      , rICLDDTdt
   real    yp, yt, yq, s1, s2, xt, &
        DUOORL
   real    tabice, tabde, dpksf, tetae
   integer nlevm1, nlevp1

   integer liftst

   integer il, jk, &
        JNR
   integer lqonof
   real    temp1         , temp2

   real    TVBEG  , REDUES , TVIRTC
   real    T0I    , TCI    , TSCALE , APRI  , &
        TOPEQ0 , TODPMX

   logical logic  , logic3

   integer IVBEG

   real     XACCES , &
        XN     , XNPM   , XBNPM  , TPREL  , &
        QSPREL , XDE    , XPRB   , BFMOD  , &
        HFCOX  , XFT    , HFREZX , &
        ANVCOV , X

   integer, dimension(NI,NLEV) :: IPRTCO
   integer, dimension(NI     ) :: KCBTOP
   integer, dimension(NI     ) :: INDCU

   real, dimension(NI,NLEV) :: ACCK
   !LS      REAL, dimension(NI,NLEV) :: HPLOGK
   real, dimension(NI,NLEV) :: CCWM
   real, dimension(NI,NLEV) :: HDCWST
   real, dimension(NI,NLEV) :: HEVAC
   real, dimension(NI,NLEV) :: HEVAPR
   real, dimension(NI,NLEV) :: HDTST
   real, dimension(NI,NLEV) :: HDQST
   real, dimension(NI,NLEV) :: STHEAT
   real, dimension(NI,NLEV) :: HFOORL
   real, dimension(NI,NLEV) :: HQSAT
   real, dimension(NI,NLEV) :: HSQ
   real, dimension(NI,NLEV) :: HU
   real, dimension(NI,NLEV) :: CUMASK
   real, dimension(NI,NLEV) :: HDCWAD
   real, dimension(NI,NLEV) :: HTP1
   real, dimension(NI,NLEV) :: HQP1
   real, dimension(NI,NLEV) :: DPK
   real, dimension(NI,NLEV) :: DLNPDT
   real, dimension(NI,NLEV) :: HDTAD
   real, dimension(NI,NLEV) :: ELOFT
   real, dimension(NI,NLEV) :: HDQAD
   real, dimension(NI,NLEV) :: HPK
   real, dimension(NI     ) :: HPS
   real, dimension(NI     ) :: COVVAP
   real, dimension(NI     ) :: COVBAR
   real, dimension(NI     ) :: DONESC
   real, dimension(NI     ) :: XSCEV
   real, dimension(NI     ) :: XKCOVB
   real, dimension(NI     ) :: DUSTAB
   real, dimension(NI     ) :: XDIFSC
   real, dimension(NI     ) :: SUPSAT
   real, dimension(NI     ) :: XFRCOA
   real, dimension(NI     ) :: XFIX
   real, dimension(NI     ) :: XCST
   real, dimension(NI     ) :: YM
   real, dimension(NI     ) :: YMMIN
   real, dimension(NI     ) :: PRCPCU
   real, dimension(NI     ) :: HFMRX
   real, dimension(NI     ) :: PRCPST
   real, dimension(NI     ) :: CUSNOW
   real, dimension(NI     ) :: STSNOW
   real, dimension(NI     ) :: HEVACU
   real, dimension(NI     ) :: HEVRST
   real, dimension(NI     ) :: HEVCST
   real, dimension(NI     ) :: REQCON
   real, dimension(NI     ) :: HUZ00
   real, dimension(NI     ) :: PRBMOD
   real, dimension(NI,NLEV) :: WORK

   !***********************************************************************

   !  III)     STATEMENT FUNCTIONS
   !           -----------------------------------------------------------

   real Z1, Z2, Z3, Z4, Z5

   !           FORMULA FOR THETHAE FROM A SERIES EXPANSION

   !           TAE = PI*T*EXP(L*Q/(CP*T)) = EXP(CONAE)*PI*T*(1+(X*(1+.5*X))
   !           WHERE X = L*Q/(CP*T)-CONAE  AND CONAE IS A CONST APPR = 0.15

   TETAE(YP,YT,YQ) = AECON * YP * YT * ( 1. + ( HLDCP * YQ / YT &
        - CONAE ) * ( 1. + 0.5 * ( HLDCP * YQ / YT &
        - CONAE ) ) )

   !           FORMULA FOR DPK/P

   DPKSF(S1,S2) = 2. * (S1-S2) / (S1+S2)

   !           DIFFERENCE OF SATURATION PRESSURE BETWEEN WATER AND ICE
   !           PROBABILITY FOR EXISTENCE OF ICE AND THE PRODUCT OF THOSE

   TABDE(XT)  = min((( HE273/XT * exp(HELDR*(T0I - 1./XT)) * &
        (1. - exp(HEDLDR*(T0I - 1./XT))) ) *9.248487), 1.0)

   !     TABDE(XT)  = HDEWI( IT0P1-INT(XT) ) * ( XT-INT(XT) ) +
   !    &             HDEWI( IT0P1-INT(XT)+1 ) * ( 1.-XT+INT(XT) )

   TABICE(XT) = max(APRI*(exp(-(((max(XT,TCI)-TCI)/TSCALE)**2)) &
        -1.)+1. , 0.0)

   !     TABICE(XT) = PRBICE( IT0P1-INT(XT) ) * ( XT-INT(XT) ) +
   !    &             PRBICE( IT0P1-INT(XT)+1 ) * (1.-XT+INT(XT) )

   !      TABFBF(XT) = TABICE(XT) * (1.-TABICE(XT)) * TABDE(XT)

   !      TABFBF(XT) = BFEFF(IT0P1-INT(XT))*(XT-INT(XT)) +
   !     1             BFEFF(IT0P1-INT(XT)+1)*(1.-XT+INT(XT))

   !           TEMPERATURE FUNCTION TO MULTIPLY HMRCU AND HMRST FOT T<273

   Z1(XT) = min(1.33*exp(-((XT-HT273)*.066)**2) , 1.0)

   Z2(XT) = abs (XT - 232.) / 18.

   Z3(XT) = Z2(XT) * (1. + Z2(XT) * (1. + 1.333 * Z2(XT)))

   Z4(XT) = Z3(XT) / (1. + Z3(XT)) * sign(1.0,XT-232.)

   Z5(XT) = 0.5*0.15*(1.07+Z4(XT))

   !      FMROFT(XT) = CVMGT( Z1(XT), Z5(XT), XT.GT.250. )

   !     FMROFT(XT) = HMROFT( IT0P1-INT(XT) ) * (XT-INT(XT) ) +
   !    &             HMROFT( 1+IT0P1-INT(XT) ) * (1.-XT+INT(XT) )

   !-----------------------------------------------------------------------




   !  IV)      VALUES OF CONSTANTS - IN SI UNITS - INCLUDING DERIVED ONES
   !           -----------------------------------------------------------

   TVBEG  = 100.
   REDUES = 0.983
   TVIRTC = 0.61
   IVBEG  = int(TVBEG) - 1

   HT273  = 273.
   HE273  = 611.
   T0I    = 1./HT273
   TCI    = 232.
   TOPEQ0 = 268.
   TODPMX = 256.
   TSCALE = (TODPMX - TCI)*sqrt(2.)

   APRI = 1./(1.-exp(-((TOPEQ0-TCI)/TSCALE)**2))


   ICLDDT = 1.
   KSCFLD = 0
   FSCFLD = 0.
   FSCFKZ = 0.
   KSCTOP = NLEV + 3
   HCP    = CPD
   HEPS   = EPS1
   HG     = GRAV
   rHG    = 1. / hg
   HLV    = CHLC
   !     HDL   =  CHLF
   HDL   = 0.
   HP0    = 1.E5
   HRD    = RGASD
   HEPSE0 = HEPS * HE273
   HKAP   = HRD/HCP
   IT0 = int(HT273)
   IT0P1 = IT0 + 1
   HKPI   = 1./HKAP
   HCPI = 1./HCP
   HEDR = HEPS/HRD
   HRTDE = HRD*HT273/HEPS
   HLDCP = HLV*HCPI
   HDLDCP = HDL*HCPI
   HRTDEL = HRTDE/HLV
   HELDR = HEDR*HLV
   HELDCP = HEPS * HLDCP
   HEDLDR = HEPS*HDL/HRD
   HELDRI = 1./HELDR
   SHCUM = 0.

   !-----------------------------------------------------------------------

   !  V)       PARAMATER VALUES IN SI UNITS
   !           ------------------------------------------------------------

   CONAE  = 0.15
   AECON  = exp(CONAE)
   CFREEZ = 0.12
   COALES = 300.
   CBFEFF = 7.0
   CTFRZ1 = 263.
   HCCU   = 5.E-4
   LQONOF = 0
   PMOIST = 3.
   HCUNRM = 0.3
   HMRCU  = 5.E-4
   HMRST  = 2.E-4
   HTAUCU = 3600.
   TANVIL = 253.
   HCST   = 1.8E-4
   HVTERM = 5.
   STPEVP = 1.6E-4
   CVPEVP = 0.6E-4
   HVSNOW = 2.
   SQVSNO = sqrt(HVSNOW)
   HKMELT = 3.E-5
   U00MAX = 0.99
   HU00 = 0.85
   RDT = 1. / DT        !AJOUT (WEI YU, 1994)
   ICLDDT2dt = ICLDDT * DT
   rICLDDT2dt = 1. / ICLDDT2dt
   ICLDDTdt = ICLDDT * DT
   rICLDDTdt = 1. / ICLDDTdt

   !-----------------------------------------------------------------------

   !  VI)      PREPARATIONS

   do JK = 1 , nlev
      do IL = 1 , NI
         PRFLX(IL,JK) = 0.0
         SWFLX(IL,JK) = 0.0
         CCWM(IL,JK)  = CWM(IL,JK)
      enddo
   enddo

   NLEVM1 = NLEV - 1
   NLEVP1 = NLEV + 1
   !     ------------------------------------------------------------
   do il = 1, ni
      HPS(il) = ( PSP(IL) + PSM(IL) ) * 0.5
      HUZ00(il) = HU00
   enddo

   do jk = 1, nlev
      do il = 1, ni
         HTP1(il,JK) = TP(IL,JK)
         HQP1(il,JK) = QP(IL,JK)
         HPK(il,JK) = S(il,JK) * HPS(il)   !Sigma local (wei YU)
      enddo
   enddo

   do il = 1, ni
      DPK(il,1) = 0.5 * ( HPK(il,2) - HPK(il,1) )
   enddo

   do JK  = 2, NLEVM1
      do il  = 1, ni
         DPK(il,JK) = 0.5 * ( HPK(il,JK+1) - HPK(il,JK-1) )
      enddo
   enddo

   do il = 1, ni
      DPK(il,NLEV) = 0.5 * ( HPK(il,NLEV) - HPK(il,NLEV-1) ) &
           + S(il,nlev+1) * HPS(il) - HPK(il,NLEV)
   enddo

   do JK = 1, NLEV
!vdir noloopchg
      do il = 1, ni
         DLNPDT(il,JK) = ( 1. / HPK(il,JK) ) * S(il,JK)  &!Sigma local
              * ( PSP(IL) - PSM(IL) ) * rDT
         HDTAD(il,JK) = ( HTP1(il,JK) - TM(il,JK) ) * rDT
         HDQAD(il,JK) = ( HQP1(il,JK) - QM(il,JK) ) * rDT
         HDCWAD(il,JK) = ( CWP(il,JK) - CCWM(il,JK) ) * rDT
      enddo
   enddo

   do il = 1, ni
      KCBTOP(il) = 3
      COVVAP(il) = 0.
      COVBAR(il) = 0.
      PRCPCU(il) = 0.
      PRCPST(il) = 0.
      CUSNOW(il) = 0.
      STSNOW(il) = 0.
      HEVACU(il) = 0.
      HEVRST(il) = 0.
      HEVCST(il) = 0.
      indcu(il) = 0
   enddo
   !      KCBTOP = 3
   hfis = 0.
   ANVCOV = 0.
   liftst = 3 + 2  !ne pas varier en fonction de maille

   !     CALCULATE SATURATED PRESSURES

   if (SATUCO) then
      call mfoqst3(HQSAT, TM,HPK,NI,NLEV,NI)
   else
      call mfoqsa3(HQSAT, TM,HPK,NI,NLEV,NI)
   endif



   do jk = 1, nlev
      do il = 1, ni
         !            temp1 = TM(il,jk)
         work(il,jk)=-(((max(TM(il,jk),tci)-tci)/tscale)**2)
      enddo
   enddo

   call vsexp(work,work,ni*nlev)

   if (hdl.ne.0.)then
      do jk = 1, nlev
         do il = 1, ni
            if (TM(il,jk) .le. ht273) then
               temp1 = TM(il,jk)
               eloft(il,jk) = hlv + hdl * &
                    max(((apri*(work(il,jk)-1.0))+1.0),0.0)
            else
               eloft(il,jk) = hlv
            endif
         enddo
      enddo
   else
      do jk = 1, nlev
         do il = 1, ni
            eloft(il,jk) = hlv
         enddo
      enddo
   endif



   !           UPDATE THE VERTICAL ARRAYS

   do JK = 1, NLEV
      do il = 1, ni

         ACCK(il,JK)  = 0.
         !LS         HTCMT(il,JK) = 0.
         !LS         HTCMTb(il,JK) = 0.
         !LS         HQCMQ(il,JK) = 0.
         !LS         HQCMQb(il,JK) = 0.
         !LS         HKSIZ(il,JK) = 0.
         CUMASK(il,JK) = 1.

         HDCWST(il,JK) = 0.
         HEVAC(il,JK) = 0.
         HEVAPR(il,JK) = 0.
         HDTST(il,JK) = 0.
         HDQST(il,JK) = 0.
         STHEAT(il,JK) = 0.
         ZSTCOV(il,JK) = 0.
         HFOORL(il,JK) = 1.
         !LS         HPLOGK(il,JK) = ALOG(HPK(il,JK))
         !LS         HPKAP(il,JK) = (HP0/HPK(il,JK)) ** HKAP
         !LS         if ( TM(il,jk) .le. ht273 ) then
         !LS            temp1 = TM(il,jk)
         !LS            eloft(il,jk) = hlv + hdl * tabice( temp1 )
         !LS         else
         !LS            eloft(il,jk) = hlv
         !LS         endif
         !***
         IPRTCO(il,JK) = 0

         HELDR = HEDR*ELOFT(il,JK)
         HLDCP = HCPI*ELOFT(il,JK)

         temp1 = TM(il,JK)
         temp2 = HPK(il,JK)
         !***
         !LS          if ( SATUCO ) then
         !LS             HQSAT(il,jk) = FOQST( temp1 , temp2 )
         !LS          else
         !LS             HQSAT(il,jk) = FOQSA( temp1 , temp2 )
         !LS          endif
         !***
         HQSAT(il,JK) = HQSAT(il,JK) * ( 1. + TVIRTC * HQSAT(il,JK) )
         HSQ(il,JK) = HELDR * HLDCP * HQSAT(il,JK) / ( TM(il,JK) ** 2)
         HU(il,JK) = AMAX1( QM(il,JK) , HQP1(il,JK) ) / HQSAT(il,JK)
         HU(il,JK) = AMIN1( HU(il,JK) , 1. )
         HU(il,JK) = AMAX1( HU(il,JK) , 0. )
         !LS         HUQSM1(il,JK) = HU(il,JK) * HQSAT(il,JK)
         !LS         MOISTN(il,JK) = ( 1. - HU(il,JK) ) ** PMOIST
         !LS         HDTCU1(il,JK) = 0.
      enddo
   enddo

   !-----------------------------------------------------------------------


   !-----------------------------------------------------------------------

   !           HERE BEGINS THE CALCULATION OF STRATIFORM CONDENSATION

   !-----------------------------------------------------------------------

   if (convec == 'OLDKUO') then
      !           TRANSVIDAGE DE ILAB DANS CUMASK
      do JK=1,NLEV
         do IL=1,NI
            CUMASK (IL,JK)                     = 1
            if(ILAB(IL,JK).eq.2) CUMASK(IL,JK) = 0
         end do
      end do
   endif


   !  11)      MODIFICATION OF BASIC THRESHOLD RELATIVE HUMIDITY, HU00
   !           ------------------------------------------------------------

   DUOORL = 0.
   if (HFIS.gt.0.) DUOORL = 0.1
   !-----------------------------------------------------------------------
   !           BEGIN OF MAIN VERTICAL LOOP
   !           ------------------------------------------------------------
   do il = 1, ni
      DONESC(il) = -1.
      COVBAR(il) =  0.
      XKCOVB(il) = -1.
   enddo

   DO1670: do JK = 1, NLEV
      do il = 1, ni
         !-----------------------------------------------------------------------



         XDIFSC(il) = -1.

         !  12)      GENERAL REQUIREMENT FOR CONDENSATION, REQCON = 1 OR 0
         !           ------------------------------------------------------------
         REQCON(il) = CUMASK(il,JK)
         DUSTAB(il) = 0.
         SUPSAT(il) = 0.
!!!!      SUPCU0 = 0.
         HUZ00(il) = HU00 - DUOORL * HFOORL(il,JK)
      enddo
      !-----------------------------------------------------------------------
      !  13)      SETTING OF SPECIFIC HU00 AND REQCON
      !           SET REQCON = 0. IN WELL MIXED BOUNDARY LAYER OR INCREASE U0
      !           WITH STABILITY IN THE LOWEST LAYERS, I.E., (JK.GE.JMIXTK)
      !           ------------------------------------------------------------

      if( jk .ge. nlev-7 ) then

         !           LET HU00 INCREASE LINEARLY WITH P FOR P/PS .GT. 0.85
         !           ------------------------------------------------------------
         do il = 1, ni
            if ( REQCON(il).eq.1. .and. HPK(il,JK)/HPS(il).gt.0.85 ) then
               HUZ00(il) = HUZ00(il) + ( U00MAX - HUZ00(il) ) &
                    * ( 1. - 6.66 * ( 1. - HPK(il,JK) / HPk(il,nlev) ) )
            endif
         enddo

      endif

      !      IF ( JK .LT. KSCTOP ) GO TO 1310
      !c
      !      if ( KSCFLD .GT. 0 .AND. JK .GT. KSCFLD ) then
      !         do 240 il = 1, ni
      !            if ( SUPSAT(il) .EQ. 0. ) then
      !               REQCON(il) = 0.
      !            endif
      !240      continue
      !      endif
      !c***
      !C
      !C           LET HU00 INCREASE WITH STABILITY
      !C           ------------------------------------------------------------
      !      IF ( JK .LT. NLEVM1 ) GO TO 1310
      !c
      !      do 250 il = 1, ni
      !         logic = REQCON(il) .ne. 0.
      !c***
      !         if ( logic .and. JK .EQ. NLEV ) then
      !            xtp = TSS(il)
      !            xsp = 1.
      !         else
      !            xtp = TM(il,jk+1)
      !c***
      !c   Cette ligne est remplacee le 27 Juin (Wei Yu)
      !c            xsp = s(jk+1)
      !            xsp = s(il,jk+1)     !sigma local
      !c***
      !         endif

      !         if ( logic ) then
      !            X = U00MAX - HUZ00(il)
      !            DUSTAB(il) = X * (1.85-0.019*(TM(il,JK)-277.)-1.67E-2
      !c***
      !c   Cette ligne est remplacee le 27 Juin (wei Yu)
      !c
      !c     &                  *(XTP-TM(il,JK-1)) / (XSP-S(JK-1)))
      !     &                  *(XTP-TM(il,JK-1)) / (XSP-S(il,JK-1)))
      !c***
      !            DUSTAB(il) = AMAX1( DUSTAB(il), 0. )
      !         endif
      !250   continue
      !1310  CONTINUE

      !      do 260 il = 1, ni
      !         HUZ00(il) = AMIN1( (HUZ00(il) + DUSTAB(il)), U00MAX )
      !260   continue


      !           LET U00 INCREASE FOR T.LT.238K
      !           ------------------------------------------------------------
      do il = 1, ni
         if ( TM(il,JK) .lt. 238. ) then
            X = 1.+ 0.15 * ( 238. - TM(il,JK) )
            DU238 = ( U00MAX - HUZ00(il) ) * ( 1. - 1. / X )
            HUZ00(il) = AMIN1( (HUZ00(il)+DU238), U00MAX )
         endif
      enddo
      !-----------------------------------------------------------------------
      !           IF U.LT.U0,SET REQCON = 0
      !           ------------------------------------------------------------
      do il = 1, ni
         if ( HU(il,JK).le.HUZ00(il) .and. SUPSAT(il).eq.0. ) then
            REQCON(il) = 0.
         endif
      enddo

      !-----------------------------------------------------------------------
      !           CONDENSATION CALCULATIONS
      !-----------------------------------------------------------------------
      !  14)      CLOUD COVER AND EVAPORATION OF PRECIPITATION
      !           ------------------------------------------------------------
      do il = 1, ni
         ZSTCOV(il,JK) = ( 1. - sqrt( ( 1.-HU(il,JK) ) &
              / ( 1.-HUZ00(il) ) ) ) * REQCON(il)
         ZSTCOV(il,JK) = AMAX1( ZSTCOV(il,JK), 0. )
      enddo

      !      IF ( KSCTOP .GT. NLEV .OR. JK .LT. KSCTOP ) GO TO 1412
      !c
      !      do 300 il = 1, ni
      !C
      !C           MODIFICATION OF SC COVER DUE TO ENTRAINMENT
      !C
      !         if ( DONESC(il) .le. 0. .and. FSCFLD .GT. 0.
      !     &         .AND. KSCFLD .EQ. JK .and. ZSTCOV(il,JK-1) .EQ. 0. ) then
      !c***
      !            ZSTCOV(il,JK) = ZSTCOV(il,JK) / ( 1.+ZSTCOV(il,JK)*FSCFLD )
      !            XDIFSC(il) = 1.
      !            DONESC(il) = 1.
      !         endif
      !300   continue
      !1412  CONTINUE

!vdir nodep
      do il = 1, ni

         PRFLX(il,JK) =  PRFLX(il,JK)+PRCPST(il)-STSNOW(il)
         SWFLX(il,JK) =  SWFLX(il,JK)+STSNOW(il)

         logic = PRCPST(il) .ne. 0. .and. SUPSAT(il) .ne. 1.
         !***
         COVVAP(il) = 0.
         if ( logic .and. jk .ne. 1 ) then
            COVVAP(il) = ( (JK-2.) * COVVAP(il) + (JK-1.) &
                 * ZSTCOV(il,JK-1) ) / ( JK - 1. + 1.E-6 )
         endif
         !***

         !         XCLEAR = 1.
         !***
         EVAPRI = 0.
         if ( logic ) then
            !***
            XRO  = HPK(il,JK) / ( HRD * TM(il,JK) )
            XROV = XRO * HVTERM
            XMIND = AMIN1( (DPK(il,JK) / (HG*XROV) ), ICLDDT2dt )
            XEVAP = 0.5*STPEVP * ( 1.-HUZ00(il)*REQCON(il)-HU(il,JK) &
                 *(1.-REQCON(il))) * XMIND * XROV
            !***
            XCLEAR = ( 1. - ZSTCOV(il,JK) ) * COVVAP(il)
            XEVAP = XEVAP * XCLEAR
            XPSQ = sqrt( PRCPST(il) ) - XEVAP
            XP = AMAX1( XPSQ, 0. )
            xp = xp * xp
            !***
            EVAPRI = PRCPST(il) - XP
            !     HEVAPR(JK) = EVAPRI/(ICLDDT*DT*XRO)
            HEVAPR(il,JK) = EVAPRI / ( XMIND * XROV )
         endif
         !***
         PRCPST(il) = PRCPST(il) - EVAPRI
         STSNOW(il) = AMAX1( 0., (STSNOW(il)-EVAPRI) )
      enddo

      !      IF ( KSCTOP .GT. NLEV .OR. JK .LT. KSCTOP ) GO TO 1430
      !c
      !      do 320 il = 1, ni
      !C
      !C           EVAPORATION OF CLOUDWATER COMING UP FROM INVERSION CLOUD
      !C           DUE TO ENTRAINMENT
      !C
      !         if ( DONESC(il) .le. 0. .and. FSCFLD .GT. 0.
      !     &         .AND. JK .EQ. KSCFLD-1 .AND. ZSTCOV(il,JK) .EQ. 0. ) then

      !            XCEVSC = 0.5 * ICLDDT2dt * FSCFKZ * DPK(il,JK+1)
      !     &              / DPK(il,JK)
      !            XSCEV(il) = ( CWM(il,JK+1) * (1.-XCEVSC) + ICLDDT2dt
      !     &                  * HDCWAD(il,JK+1) ) / (1.+XCEVSC)

      !            XSCEV(il) = XCEVSC * ( XSCEV(il)+CWM(il,JK+1) )
      !     &                 * rICLDDT2dt
      !c      XSCEV1 = XSCEV
      !            HEVAC(il,JK) = HEVAC(il,JK) + XSCEV(il)
      !            XSCEV(il) = XSCEV(il) * DPK(il,JK) / DPK(il,JK+1)
      !         endif

      !320   continue
      !1430  CONTINUE

      !-----------------------------------------------------------------------

      !  15)      ACCESSION OF VAPOUR, AND FINAL SETTING OF HEATING/COOLING
      !           ------------------------------------------------------------

      do il = 1, ni
         if ( REQCON(il) .ne. 0. ) then

            HLDCP = HCPI * ELOFT(il,JK)
            XACCES = HDQAD(il,JK) - HU(il,JK) * HSQ(il,JK) &
                 * HDTAD(il,JK) / HLDCP + QM(il,JK) * DLNPDT(il,JK)
            !-----------------------------------------------------------------------
            !           CALCULATION OF HEATING RATE COMPOSED OF RELEASE OF LATENT
            !           HEAT AND COOLING DUE TO EVAPORATION OF PRECIPITATION
            !           ------------------------------------------------------------
            XN = 2. * ( 1.-HUZ00(il) ) * ( 1.-ZSTCOV(il,JK)+1.E-3 )
            XNPM = XN + CCWM(il,JK) / ( ZSTCOV(il,JK) * HQSAT(il,JK) &
                 + 1.E-6 )
            XBNPM = XNPM - XN + ZSTCOV(il,JK) * XN
            STHEAT(il,JK) = (XBNPM * XACCES - XN *HEVAPR(il,JK)) / &
                 ((1. + HU(il,JK) * HSQ(il,JK)) * XNPM)
         endif

         if ( SUPSAT(il) .eq. 0. .and. &
              (STHEAT(il,JK) + CCWM(il,JK)*rICLDDT2dt+ HDCWAD(il,JK)) &
              .lt. 0. ) then
            REQCON(il) = 0.
         endif
      enddo
      !-----------------------------------------------------------------------
      !           FINAL SETTING OF HEATING/COOLING AND CLOUD COVER
      !           ------------------------------------------------------------
      do il = 1, ni
         STHEAT(il,JK) = STHEAT(il,JK) * REQCON(il) - &
              (1.-REQCON(il)) * (HEVAC(il,JK) + HEVAPR(il,JK))
         ZSTCOV(il,JK) = ZSTCOV(il,JK) * REQCON(il)
      enddo


      !           CHECK IF SUPER SATURATION MAY RESULT AND IF SO, ADJUST

!vdir nodep
      do il = 1, ni
         !***
         QPREL = 0.                               !pour la vectorisation
         QSPREL = 0.                              !pour la vectorisation
         !***
         if ( CUMASK(il,JK) .eq. 1. ) then
            HLDCP = HCPI * ELOFT(il,JK)

            QPREL = QM(il,JK) + DT * ( HDQAD(il,JK)-STHEAT(il,JK) )
            TPREL = TM(il,JK) + DT * ( HDTAD(il,JK) + STHEAT(il,JK) &
                 * HLDCP )
            QSPREL = HQSAT(il,JK) + HSQ(il,JK) * ( TPREL - TM(il,JK) ) &
                 / HLDCP - QM(il,JK) * DLNPDT(il,JK) * DT
            QSPREL = AMAX1( QSPREL, 0.0 )
         endif

         !***
         QINCR = 0.                               !pour la vectorisation
         !***
         if ( CUMASK(il,JK) .eq. 1. .and. QPREL .gt. QSPREL ) then
            QINCR = ( QPREL - QSPREL ) / ( DT * ( 1. + HSQ(il,JK) ) )
            STHEAT(il,JK) = STHEAT(il,JK) + CUMASK(il,JK) * QINCR
         endif

         if ( CUMASK(il,JK) .eq. 1. .and. QPREL .gt. QSPREL &
              .and. STHEAT(il,JK) .gt. 0. ) then
            REQCON(il) = 1.
            CCWM(il,JK) = AMAX1( CCWM(il,JK) - QINCR, 0. )
            ZSTCOV(il,JK) = 1.
         endif
      enddo

      do il = 1, ni
         if ( REQCON(il) .eq. 0. &
              .and. CWP(il,JK) * CUMASK(il,JK) .gt. 0. ) then

            STHEAT(il,JK) = STHEAT(il,JK) - CWP(il,JK) * rDT
            HDCWST(il,JK) = - CWP(il,JK) * rDT
         endif
      enddo
      !-----------------------------------------------------------------------

      !  16)      PREDICTION OF CLOUD WATER CONTENT, CW, BY AN IMPLICIT SCHEME
      !           REQUIRING AN ITERATIVE PROCEDURE.
      !           NEWTON - RAPHSON ITERATIVE METHOD IS USED
      !           ------------------------------------------------------------
      DO370: do il = 1, ni
         logic3 = REQCON(il) .ne. 0.

         !           CALCULATE FACTORS FOR COALESCENCE HFCOX, FREEZING HFREZX,
         !           REDUCTION OF HMRST AT LOW TEMPS, HFMRX
         !           ------------------------------------------------------------
         if ( logic3 ) then
            IPRTCO(il,JK) = 3
            XKCOVB(il) = XKCOVB(il) + 1.
            COVBAR(il) = ( XKCOVB(il) * COVBAR(il) + ZSTCOV(il,JK) ) &
                 / ( XKCOVB(il) + 1. )
         endif

         !***

         !           MODIFIED PROBABILITY OF ICE RESULTING FROM ICE IN PRECIP
         !           COMING FROM ABOVE
         !***
         if ( logic3 .and. TM(il,JK) .le. HT273 ) then
            temp1 = TM(il,JK)
            XDE = TABDE(temp1)
            XPRB = TABICE(temp1)
            !            XPRB = max(((apri*(work(il,jk)-1.0))+1.0),0.0)

            if (temp1.gt.250.) then
               XFT = Z1(temp1)
            else
               XFT = Z5(temp1)
            endif

         else
            XDE = 0.
            XPRB = 0.
            xft = 1.
         endif

         if ( logic3 ) then
            PRBMOD(il) = XPRB + (1.-XPRB) * STSNOW(il) &
                 / ( PRCPST(il)+1.E-7 )
            BFMOD = PRBMOD(il) * (1.-XPRB) * XDE

            STSNOW(il) = AMAX1( 0., STSNOW(il) )
            HFCOX = 1. + COALES * sqrt( PRCPST(il) + ACCK(il,JK) &
                 * rICLDDTdt )

            XFT = AMAX1( XFT, 3.E-2 )

            HFREZX = 1. + CBFEFF * BFMOD
            HFREZX = HFREZX * ( 1. + CFREEZ * (1.-XFT) / XFT )

            XFRCOA(il) = HFCOX * HFREZX
            HFMRX(il) = HMRST * XFT / HFCOX
         endif


         !-----------------------------------------------------------------------
         !           SPECIAL TREATMENT FOR T.LT.236K
         !           -------------------------------
         if ( logic3 .and. TM(il,JK) .lt. 236. &
              .and. TM(il,JK) .gt. 232. ) then
            XFRCOA(il) = 0.25 * ( XFRCOA(il) * ( TM(il,JK) - 232. ) &
                 + 5. * ( 236. - TM(il,JK) ) )
         endif
         !***
         if ( logic3 .and. TM(il,JK) .le. 232. ) then
            XFRCOA(il) = 5.
         endif
         !-----------------------------------------------------------------------
         !           CALCULATE THE FIXED PART OF EQU AND NORMALIZE BY B * MR
         !           ------------------------------------------------------------
         if ( logic3 .and. XDIFSC(il).gt.0. ) then
            XSCDIF = XSCEV(il)
         else
            XSCDIF = 0.
         endif

         if ( logic3 ) then
            XFIX(il) = ((cwp(il,jk)+CCWM(il,JK)) + &
                 ICLDDT2dt * ( HDCWAD(il,JK) &
                 - XSCDIF + STHEAT(il,JK) ) ) &
                 / ( 2. * ( ZSTCOV(il,JK)+1.E-2 ) * HFMRX(il) )

            !           THE CONVERSION RATE TIMES 2*DT

            XCST(il) = 0.5 * HCST * XFRCOA(il) * ICLDDT2dt

            !           FIRST GUESS IS THE NORMALIZED M(T-DT)

            YM(il) = 0.5*(ccwm(il,JK)+cwp(il,JK)) / &
                 ( ( ZSTCOV(il,JK)+1.E-2 ) * HFMRX(il) )

            !           TO MAKE M(T+DT)GE.0, THE FINAL SOLUTION
            !           OF YM HAS TO BE YM.GE.M(T-DT)/(2.*B*MR)

            YMMIN(il) = 0.5 * YM(il)
         endif
      enddo DO370

      !           NEWTON - RAPHSON LOOP
      !           ------------------------------------------------------------
      do JNR = 1,2
         do il = 1, ni
            if ( REQCON(il) .ne. 0. ) then
               YMN = AMIN1( YM(il), 5. )
               XXP = exp( -YMN * YMN )
               XHJ = 1.+ XCST(il) * ( 1. - XXP )
               XF = YM(il) * XHJ - XFIX(il)
               XFPRIM = XHJ + 2. * XCST(il) * YM(il) * YM(il) * XXP
               YM(il) = AMAX1( (YM(il)-XF/XFPRIM), YMMIN(il) )
            endif
            !***
            !                ca se converge au bout de 5 iterations
            !      XERMAX = ABS(XF)
            !      IF (XERMAX.LT.1.E-5) GO TO 1635
            !***
         enddo
      enddo

      !1635 CONTINUE
      !***
      do il = 1, ni
         logic3 = REQCON(il) .ne. 0.

         if ( logic3 ) then
            XCWP1 = 2. * YM(il) * ( ZSTCOV(il,JK)+1.E-2 ) * HFMRX(il) &
                 - CWp(il,JK)
            XCWP1 = AMAX1( XCWP1, 0. )
            HDCWST(il,JK) = ( XCWP1 - CCWM(il,JK) ) * rICLDDT2dt &
                 - HDCWAD(il,JK)
         endif

         if ( logic3 .and. abs( HDCWST(il,JK)) .lt. 1.E-16 ) then
            HDCWST(il,JK) = 0.
         endif

         !           RATE OF PRECIPITATION
         !           ------------------------------------------------------------

         if ( logic3 ) then

            XPRADD = DPK(il,JK) * rHG &
                 * AMAX1( (STHEAT(il,JK)-HDCWST(il,JK)), 0. )
            PRCPST(il) = PRCPST(il) + XPRADD
            STSNOW(il) = AMIN1( (STSNOW(il) + PRBMOD(il) * XPRADD), &
                 PRCPST(il) )
         endif
         !***
         if ( logic3 .and. PRCPST(il) .eq. 0. ) then
            COVBAR(il) = 0.
         endif
      enddo

      !           ADJUST THE HEATING FOR THE THETA EQU AND UPDATE
      !           THE MOISTURE TENDENCY
      !           ------------------------------------------------------------

      do il = 1, ni
         HLDCP = HCPI * ELOFT(il,JK)
         HDTST(il,JK) = STHEAT(il,JK) * HLDCP
         HDQST(il,JK) = -STHEAT(il,JK)
      enddo

      !           CALCULATE POSSIBLE MELTING

      !      DTMELT = 0.
      do il = 1, ni

         if ( TM(il,JK) .gt. HT273 .and. STSNOW(il) .gt. 0. ) then
            XRO = HPK(il,JK) / ( HRD * TM(il,JK) )
            XROV = XRO * HVSNOW
            XMIND = AMIN1( (DPK(il,JK) / (HG*XROV) ), ICLDDT2dt )
            DMELT = 0.5 * HKMELT * ( TM(il,JK)-HT273 ) * XMIND * XROV
            SQMELT = sqrt( STSNOW(il) ) - DMELT
            XMELT = AMAX1( SQMELT, 0. )
            xmelt = xmelt * xmelt
            DSNOW = STSNOW(il) - XMELT
            STSNOW(il) = XMELT
            DMELT = DSNOW / ( ICLDDT2dt * XRO )
            DTMELT = HDLDCP*DMELT
         else
            DTMELT = 0.
         endif

         HDTST(il,JK) = HDTST(il,JK) - DTMELT
      enddo

      !           ACCUMULATE EVAPORATED STRATIFORM PRECIPITATION AND
      !           CLOUD WATER FOR TABLE OUTPUT
      !           ------------------------------------------------------------
      do il = 1, ni
         HEVRST(il) = HEVRST(il) + HEVAPR(il,JK) * DPK(il,JK) * rHG
         HEVCST(il) = HEVCST(il) + HEVAC(il,JK) * DPK(il,JK) * rHG
      enddo

   enddo DO1670
   !           END OF MAIN VERTICAL LOOP
   !-----------------------------------------------------------------------

   !  17)      ACCUMULATE STRATIFORM PRECIPITATION.
   !           ------------------------------------------------------------
   !           ------------------------------------------------------------

   !           RETURN CALCULATED TENDENCIES TO THE MAIN ARRAYS, INCLUDING
   !           ADJUSTMENT OF THE HEATING FOR THE THETA EQU.
   !           ------------------------------------------------------------

   do il = 1, ni
      !******************************
      ! small precip rates --> 0.0  *
      !******************************
      if(prcpst(il) .lt. 1.0E-05) then
         zsrr(il) = 0.0
         zssr(il) = 0.0
      else
         ZSRR(IL) = PRCPST(il) - STSNOW(il)
         ZSSR(IL) = STSNOW(il)
      endif
   enddo

   do JK = 1, NLEV
      do il = 1, ni
         ZTE    (il,JK) = HDTST (il,JK)
         ZQE    (il,JK) = HDQST (il,JK)
         ZCWE   (il,JK) = HDCWST(il,JK)
      enddo
   enddo

   prflx(:,nlev+1)=prflx(:,nlev)
   swflx(:,nlev+1)=swflx(:,nlev)

   !           HERE ENDS THE CALCULATION OF STRATIFORM CONDENSATION
   !***********************************************************************

   return
end subroutine skocon5
