!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

!/@*
subroutine cslm_main(bus, bussiz, ptsurf, ptsurfsiz, lcl_indx, trnch, kount, n, m, nk )
   !use sfclayer_mod, only: sl_prelim,sl_sfclayer,SL_OK
   use sfclayer, only: sl_prelim,sl_sfclayer,SL_OK
   use cpl_itf     , only: cpl_update
   use sfc_options
   use sfcbus_mod
   use mu_jdate_mod, only: jdate_day_of_year, mu_js2ymdhms
   implicit none
#include <arch_specific.hf>
   !@Object    CANADIAN SMALL LAKE MODEL
   !             Computes fluxes, lake surface temperatures, and other lake properties
   !             (e.g. temperature profile, mixed layer depth, mixed layer TKE etc)
   !             Fluxes are weighted averages based on: open water, bare ice cover, and snow covered
   !             ice fractions.  Snow physics routines taken from CLASS Ver. 3.6.1
   !@Arguments
   !             - Input/Output -
   ! BUS         Bus for the WATER surface scheme
   !             - Input -
   ! BUSSIZ      dimension of bus
   ! PTSURF      surface pointers
   ! PTSURFSIZ   dimension of ptsurf
   ! TRNCH       row number
   ! DELT        timestep
   ! KOUNT       timestep number
   ! N           horizontal dimension (row length)
   ! M           horizontal dimensions of fields
   !             (not used for the moment)
   ! NK          vertical dimension

   integer bussiz, kount, trnch, N, M, NK
   real, target :: bus(bussiz)
   integer ptsurfsiz
   integer ptsurf(ptsurfsiz), lcl_indx(2,n)

   !@Author M. MacKay (March 2018)
   !@Notes
   !
   !
   !@Revisions
   !
   !*@/
   !include "thermoconsts.inc"
   include "sfcinput.cdk"
   !include "dintern.inc"
   !include "fintern.inc"

   integer surflen
#define x(fptr,fj,fk) ptsurf(vd%fptr%i)+(fk-1)*surflen+fj-1

   integer, parameter ::INDX_SFC = INDX_LAKE

!========================================================================================
!     * FIXED PARAMETERS FOR CSLM
! Sediment properties (fixed for now)
      REAL, PARAMETER :: DSED = 5.0       !THICKNESS OF SEDIMENT LAYER (m) 
      REAL, PARAMETER :: TSEDB = 277.15   !CONSTANT TEMP AT BOTTOM OF SEDIMENT (K)
      REAL, PARAMETER :: HCPSED = 2.13E6  !HEAT CAPACITY OF SEDIMENT (J/m3.K)
!     REAL, PARAMETER :: TCSED = 2.5      !THERMAL COND.SEDIMENT (J/m3.K) - zero for adiabatic lower B.C.
      REAL, PARAMETER :: TCSED = 0.0      !THERMAL COND.SEDIMENT (J/m3.K) - zero for adiabatic lower B.C.
      REAL, PARAMETER :: QSED  = 1.0      !SEDIMENT HEAT FLUX MULTIPLIER
! Other fixed parameters
      include "cslm.cdk"                  !Defines NLAKMAX, lake parms
      INTEGER, PARAMETER :: IGL = 1             !CLASS parameter
      INTEGER, PARAMETER :: NBS = 4             !CLASS parameter
      INTEGER, PARAMETER :: ISNOW = 1           !CLASS parameter
      INTEGER, PARAMETER :: ISLFD = 0           !CLASS parameter
      INTEGER, PARAMETER :: IPCP = 1            !CLASS parameter: rain/snow partitioning
      INTEGER, PARAMETER :: IZREF = 1           !CLASS parameter
      INTEGER, PARAMETER :: ITG = 1           !CLASS parameter
      INTEGER, PARAMETER :: IALS = 0           !CLASS parameter: snow albedo computed, not prescribed
      INTEGER, PARAMETER :: ISNOALB = 0        !CLASS parameter: 4-band radiation scheme disables for now
      REAL, PARAMETER :: SNOLIM  = 0.10      !
      logical, parameter :: WATER_TDIAGLIM = .false.

! CSLM common blocks
      REAL TKECN,TKECF,TKECE,TKECS,HDPTHMIN, DUMAX, QAMIN,              &
           TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,TKECL
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,TKEMIN,DELMAX, &
                      DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,TKECL,DUMAX, QAMIN                                      
      DATA  TKECN,      TKECF,      TKECE,      TKECS,      TKECL      &
      /     1.33,       0.25,       1.15,       0.20,      0.2350/ 
      DATA  HDPTHMIN, TKEMIN,  DELMAX,  DELMIN,  DELZLK,  DELSKIN      &
      /     0.5,      1.0E-12, 5.0,     0.5,     0.5,     0.050 /   
      DATA  DHMAX,    DUMAX,   EMSW,  QAMIN                             &
      /     2.0,      0.1,     0.97,  1.0E-5   /                                       

! CLASS common blocks
!     * THE FOLLOWING COMMON BLOCKS ARE DEFINED SPECIFICALLY FOR USE 
!     * IN CLASS, AND ARE RETAINED HERE TO MAINTAIN COMPATIBILITY
      REAL DELTAT,TFREZ, CLHVAP, CPDzz,                                   &
           RGAS,RGASVzz,GRAVzz,SBC,VKC,CTzz,VMIN,TCW,TCICE,TCSAND,TCCLAY, &
           TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,      &
           HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT
      COMMON /CLASS1/ DELTAT,TFREZ
      COMMON /CLASS2/ RGAS,RGASVzz,GRAVzz,SBC,VKC,CTzz,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,             &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,              &
                      TCGLAC,CLHMLT,CLHVAP
      DATA      CTzz,       TCW,        TCICE,      TCSAND,    TCCLAY,  TCOM     &
             /  1.15E-3,    0.57,       2.24,       2.5,       2.5,     0.25  /
      DATA      TCDRYS,     RHOSOL,     RHOOM                                    &
             /  0.275,      2.65E3,     1.30E3  /    
      DATA      HCPW,       HCPICE,     HCPSOL,    HCPOM                         &
             /  4.187E6,    1.9257E6,   2.25E6,    2.50E6 /
      DATA      HCPSND,     HCPCLY,     SPHW,      SPHICE,  SPHVEG               &
             /  2.13E6,     2.38E6,     4.186E3,   2.10E3,  2.70E3 /
      DATA      RHOW,       RHOICE,     TCGLAC                                   &
             /  1.0E3,      0.917E3,    2.24                       /

! these parameters kept local for now
      REAL CKARM,CGRAV,ASXlocal,ANGMAXlocal,BSlocal,FACTNlocal,HMINlocal
      REAL DELTAlocal, ASlocal, CIlocal, BETAlocal
      COMMON /PHYCON/ DELTAlocal,CGRAV,CKARM,CPDzz
      COMMON /CLASSD2/ ASlocal,ASXlocal,CIlocal,BSlocal,BETAlocal,FACTNlocal,HMINlocal,ANGMAXlocal              
 
      DATA      DELTAlocal,    ASlocal,       ASXlocal,      ANGMAXlocal           &
             /  0.608,         12.0,          4.7,           0.85   /
      DATA      CIlocal,       BSlocal,       BETAlocal,     FACTNlocal, HMINlocal & 
             /  40.0,          1.0,           1.0,           1.2,        40.   /

!-- LOCAL CONSTANTS AND VARIABLES.
      INTEGER I,J,K,JBOT,JTOP,NLEV,JMIX,JMIXP1
      REAL Z,DAY,DECL,HOUR,TKE_OLD,ZBOT,TBOT,WEDB,COSZ,                  &
           Z2,TMIX,ZMIX,HTOP,HBOT,ZTOP,TTOP,                             &
           RHO1,RHO2,TC1,TC2,TBAR,ATTEN1,ATTEN2,ATTEN3,DELTHRM,HCAP,     &
           EAVAIL,EHEAT,ECOOL,EFLUX,ZICE,NEWICE,TC,ICELIM,               &
           FACTM,FACTH,RATIOM,RATIOH
      REAL TSEDT,FFLXSED,ICEMASS
      REAL ICECAP,PORE,MASSL,SNO0,ALPHA,ETA,ZSNOW0,ICE0,BBLFAC
      REAL MASS,TC0,RHO0,RHOIW,ICETOP,ICEBOT,ICETOP0,ICEBOT0,ROFICE,     &
           HFREZ,HWARM,HFLX,HFLXI,HFLXS,HCPRAT,QMLTS,XXX
      INTEGER ILG,IL1,IL2

! ----* Solar Rad variables for COS(ZenithAngle) *----------------------------------------------
!
      integer(INT64), parameter :: MU_JDATE_HALFDAY = 43200
      integer yy, mo, dd, hh, mn, sec
      REAL HZ, HZ0, JULIEN

! ----* INTERNAL ARRAYS *----------------------------------------------
!
      REAL,DIMENSION(N) :: Q0,DISS,BFLX,FQU,FSHEAR,FENTRA,TRAN,          &
                             CQ1A, CQ1B, CQ2A, CQ2B, CQ3A, CQ3B,         &
                             CQ1BI, CQ2BI, CQ3BI, CDHL, CDML,      &
                             HRPCPL,HSPCPL,HCONV, ZPONDL
      REAL,DIMENSION(N) :: TCMIX, VA, G0, F0, USTAR, LKICEH0,            &
                          CSZ, DEGLAT, DEGLON, ZDIAGM, ZDIAGH,           &
                          SUNOTHER1, SUNOTHER2, SUNOTHER3, SUNOTHER4
      REAL,DIMENSION(N) :: ASVDAT, ASIDAT                 ! assigned snow albedos - not used
      REAL,DIMENSION(NLAKMAX) :: TLAK0
      REAL   FLS   (N),  HCPSNO(N),  ZSNOW (N), ALBW  (N),  ALBI  (N)
      REAL   ALVSSN(N),  ALIRSN(N),  ALVSSC(N), ALIRSC(N),               &
             ALVSG (N),  ALIRG (N),  TRSNOWC(N),ALVSL (N),  ALIRL (N)
!                                                                       
      REAL   GCONSTS(N), GCOEFFS(N), CPHCHS(N),  TCSNOW(N),              &
             ZRSLDM(N),  ZRSLDH(N),  ZRSLFM(N),  ZRSLFH(N),              &
             ZDSLM (N),  ZDSLH (N),  ZOSCLM(N),  ZOSCLH(N),              &
             ZOMLNS(N),  ZOELNS(N),  ZOM   (N),  ZOH   (N),              &
             TVIRTA(N),  TPOTA (N),  CRIBS (N),  CEVAP (N),              &
             TSTART(N),  FCOR  (N),  PCP   (N),  FTOT  (N),              &
             TSURF (N),  QSURFL(N),  TSNBOT(N),  FSCALL(N),              &
             RPCN  (N),  TRPCN (N),  SPCN  (N),  TSPCN (N)   
!
      REAL   QMELT (N),  RHOMAX(N),                                &
             RALB  (N),  WLOST (N),  SN    (N),  SL    (N),        &
             RN    (N),  RL    (N),                                &
             TRUNOF(N),  RUNOFF(N),  TOVRFL(N),  OVRFLW(N)         
!
      REAL   HTC(N,IGL), THLIQ(N,IGL), THLMIN(N,IGL), DELZW(N,IGL)
      INTEGER  IWATER(N), ISAND (N,IGL)

      ! For water balance calculations
      REAL,DIMENSION(N) :: HLAKSILM1, SNOM1, WSNOWM1
                                                                                  
!     * TSOLVE OUTPUT ARRAYS.                                                  
      REAL QSWNS (N),    QLWOS (N),    QTRANSL(N),   QSENSS(N),    &
           QEVAPS(N),    EVAPS (N),    TZEROS(N),    QZEROS(N),    &
           GSNOW (N),    GZEROSL(N),   QMELTS(N),    CDHS  (N),    &
           CDMS  (N),    RIBS  (N),    CFLUXS(N),    FTEMPS(N),    &
           FVAPS (N),    ILMOS (N),    UES   (N),    HS    (N)  
      INTEGER  ITERCT(N,6,50),   IEVAP (N)                  
                                                                        
!     * RADIATION BAND-DEPENDENT ARRAYS.                                          
      REAL TRSNOWG(N,NBS), ALSNO(N,NBS), TRTOP (N,NBS),           &
           FSDB   (N,NBS), FSFB (N,NBS), FSSB  (N,NBS) 
                                                                        
!     * SNOALBA arrays not currently used
      REAL REFSN (N), BCSN (N)

!     * INTERNAL WORK ARRAYS FOR TSOLVE.                                           
      REAL TSTEP (N),    TVIRTS(N),    EVBETA(N),    Q0SAT (N),    &
           RESID (N),    DCFLXM(N),    CFLUXM(N),    WZERO (N),    &
           A     (N),    B     (N),                                &
           LZZ0  (N),    LZZ0T (N),    FM    (N),    FH    (N) 
      INTEGER  ITER  (N),  NITER (N),  JEVAP (N), KF   (N)                                  

! ----* LAKE MODEL VARIABLES *----------------------------------------
!
      INTEGER,DIMENSION(N) :: NLAK
      REAL,DIMENSION(N) :: HLAK, LLAK, BLAK, RC1, RC2, RUNFRAC
      REAL,DIMENSION(N) :: VPD, TADP, PADRY, RHOAIR, RHOSNI, R,TR,S,TS,    &
                           RRATE, SRATE, PCPIN
      REAL,DIMENSION(N,NLAKMAX) :: QFLX, FFLX
 
      real,pointer,dimension(:) ::  T0, TKE, HDPTH,LKICEH,SNICEH,        &
           EXPW,DTEMP,DELU,GRED,RHOMIX,TSED,ROFICEH, SNO, RHOSNO, TSNOW, ALBSNO, WSNOW, &
           QSENS, QEVAP, ALVS, QSURF, HLAKSIL, ZGRIDAREA, ZROFINLAK, ZLAKEFR, ZLAKD, FICE, &
           ZLAKEAREA
      real,pointer,dimension(:,:) ::  TLAK
      REAL,POINTER,DIMENSION(:) :: LSTD, LSTF, LFXI, LFXO, zevlak

! ----* INPUT FORCING DATA   *----------------------------------------
      real,pointer,dimension(:) ::  PRES, QA, TA, UWIND, VWIND, QSWIN, QLWIN, RT
      real,pointer,dimension(:) ::  ZDLAT, ZDLON, ZREFM, ZREFH 
! ----* SFC LAYER   *----------------------------------------
      real,pointer,dimension(:) ::  z0h, z0m, cmu, ctu, hst_wat, ilmo_wat, th, zfrv
      real,pointer,dimension(:) ::  zalfaq, zalfat, zftemp, zfvap, zrunofftot
      real,pointer,dimension(:) ::  zqdiag, ztdiag, ztsurf, ztsrad, zudiag, zvdiag
      real,pointer,dimension(:) ::  zqdiagtyp, ztdiagtyp, zudiagtyp, zvdiagtyp
      real,pointer,dimension(:) ::  zsnodp

! ----* DIAGNOSTIC OUTPUT FIELDS *-------------------------------------
      REAL,DIMENSION(N) :: FSGL, FLGL, HFSL, HEVL, HMFL, HTCL, DRAGL,       &
                             FSGS, FLGS, HFSS, HEVS, HMFN, HTCS, DRAGS,     &
                             PCPL, QFL, PCPN, QFN, ROFN,                    &
                             CDH, CDM, TFLUX, EVAP, QFLUX,    &
                             EVPPOT, EVAPB, GT, DRAG,                &
                             ST, SU, SV, SQ, SH, QLWAVG, ALIR,        &
                             ILMO, HBL, UE, FTEMP, FVAP, OSW, CDHNL, vdir
!----------------------------------------------------------------------------------------
! MODULES
!----------------------------------------------------------------------------------------
      EXTERNAL CLASSI, XIT, EQNST, NDRCOEFL, FREECONV, TSOLVL, LKTRANS, MIXLYR
      EXTERNAL SNOALBA, TLSPREP, TSOLVE, TLSPOST, SNOVAP, TMELT, SNINFL, SNOALBW, SNOADD, DRCOEF, FLXSURFZ

!----------------------------------------------------------------------------------------
! CLASS common block variables whose values are taken from thermoconsts.inc or set elsewhere
      TFREZ=0.2731500000000e+03    !thermoconsts TCDK=0.2731500000000e+03  ! passage k a c         C
      SPHAIR=0.1005460000000e+04   !thermoconsts CPD = 0.1005460000000e+04  ! chal. spec. air sec   J kg-1 K-1
      CPDzz=SPHAIR
      CLHMLT=0.3340000000000e+06   !thermoconsts CHLF = 0.3340000000000e+06  ! ch. lat. fusion       J kg-1
      CLHVAP=0.2501000000000e+07   !thermoconsts CHLC = 0.2501000000000e+07  ! ch. lat. condens.(0C) J kg-1
      CGRAV=0.9806160000000e+01    !only for CLASS !thermoconsts GRAV = 0.9806160000000e+01  ! acc. de gravite       m s-2
      GRAVzz=CGRAV                 !only for CLASS
      RGASVzz=0.4615100000000e+03  !thermoconsts RGASV = 0.4615100000000e+03  ! cte gaz - vap eau     J kg-1 K-1
      RGAS=0.2870500000000e+03     !thermoconsts RGASD = 0.2870500000000e+03  ! cte gaz - air sec     J kg-1 K-1
      SBC=0.5669800000000e-07      !thermoconsts STEFAN = 0.5669800000000e-07  ! cte stefan-boltzmann  J m-2 s-1 K-4
      CKARM=0.4000000000000e+00    !only for CLASS !thermoconsts KARMAN = 0.4000000000000e+00  ! cte de von karman
      VKC=CKARM                    !only for CLASS
      DELTAT=DELT
      VMIN=VAMIN          !from sfc_options
! CLASS variables retained for now

   ILG=N
   IL1=1
   IL2=N
   DO I=1,N
     ZDIAGM(I) = ZU
     ZDIAGH(I) = ZT
   END DO
!=============================================================================================
   SURFLEN = M
! Lake variables
!! T0       (1:n) => bus( x(twater,1,1)       : )  !lake sfc temp         (K)
   T0       (1:n) => bus( x(lst,1,1)          : )  !lake sfc temp         (K)
   TKE      (1:n) => bus( x(tke,1,1)          : )  !lake mixed layer tke  (m2/s2)
   HDPTH    (1:n) => bus( x(hdpth,1,1)        : )  !lake mixed layer depth (m)
   LKICEH   (1:n) => bus( x(lkiceh,1,1)       : )  !total lake ice thickness    (m)
   SNICEH   (1:n) => bus( x(sniceh,1,1)       : )  !lake snow-ice thickness (m)
   EXPW     (1:n) => bus( x(expw,1,1)         : )  !volume expansivity of lake water (K^-1)
   HLAKSIL  (1:n) => bus( x(hlaksil,1,1)      : )  !height of lake water above sill at outlet (m)
   DTEMP    (1:n) => bus( x(dtemp,1,1)        : )  !temp jump across bottom of mixed layer (K)
   DELU     (1:n) => bus( x(delu,1,1)         : )  !current jump across bottom of mixed layer (m/s)
   GRED     (1:n) => bus( x(gred,1,1)         : )  !reduced gravity across bottom of mixed layer (m/s2)
   RHOMIX   (1:n) => bus( x(rhomix,1,1)       : )  !mean density of mixed layer (kg/m3)
   TSED     (1:n) => bus( x(tsed,1,1)         : )  !temp of sediment layer (K)
   ROFICEH  (1:n) => bus( x(roficeh,1,1)      : )  !lake ice thickness from snowmelt runoff (m) TURNED OFF FOR NOW
   TLAK     (1:n,1:NLAKMAX) => bus( x(tlak,1,1) : )!lake temp profile
   ZROFINLAK (1:n) => bus( x(rofinlak,1,1)      : ) !runoff input to lake water balance based on lake fraction (mm or kg/m2)
   FICE     (1:n) => bus( x(ficl,1,1)          : ) !fractional lake ice cover
   LSTD     (1:n) => bus( x(lstd,1,1)          : ) !lake initial water storage (kg/m2)
   LSTF     (1:n) => bus( x(lstf,1,1)          : ) !lake final water storage (kg/m2)
   LFXI     (1:n) => bus( x(lfxi,1,1)          : ) !accumulated lake water flux in (kg/m2)
   LFXO     (1:n) => bus( x(lfxo,1,1)          : ) !accumulated lake water flux out (kg/m2)
   zevlak   (1:n) => bus( x(evlak,1,1)         : ) !lake evaporation (kg/m2)
! Snow on lake variables
   SNO      (1:n) => bus( x(snol,1,1)          : )  !snow on ice 
   RHOSNO   (1:n) => bus( x(rhosnol,1,1)       : )  !density of snow on ice 
   TSNOW    (1:n) => bus( x(tsnowl,1,1)        : )  !temp. of snow on ice 
   ALBSNO   (1:n) => bus( x(albsnol,1,1)       : )  !albedo of snow on ice 
   WSNOW    (1:n) => bus( x(wsnowl,1,1)        : )  !liq. water content of snow on ice 
! Atmospheric forcing data
   if (atm_tplus) then
      PRES     (1:n) => bus( x(pplus,1,1)       : )   !sfc pressure previous timestep (Pa)
      QA       (1:n) => bus( x(huplus,1,nk)     : )   !sfc sp. humidity previous timestep (kg/kg)
      TA       (1:n) => bus( x(tplus,1,nk)      : )   !sfc air temp (K)
      UWIND    (1:n) => bus( x(uplus,1,nk)      : )   !sfc wind speed HOPEFULLY IN M/S NOT KNOTS!
      VWIND    (1:n) => bus( x(vplus,1,nk)      : )   !sfc wind speed HOPEFULLY IN M/S NOT KNOTS!
      th       (1:n) => bus( x(thetaap,1,1)     : )   !potential air temp. at first level (K?)     
   else
      PRES     (1:n) => bus( x(pmoins,1,1)       : )   !sfc pressure previous timestep (Pa)
      QA       (1:n) => bus( x(humoins,1,nk)     : )   !sfc sp. humidity previous timestep (kg/kg)
      TA       (1:n) => bus( x(tmoins,1,nk)      : )   !sfc air temp (K)
      UWIND    (1:n) => bus( x(umoins,1,nk)      : )   !sfc wind speed HOPEFULLY IN M/S NOT KNOTS!
      VWIND    (1:n) => bus( x(vmoins,1,nk)      : )   !sfc wind speed HOPEFULLY IN M/S NOT KNOTS!
      th       (1:n) => bus( x(thetaa,1,1)       : )   !potential air temp. at first level (K?)
   endif
   QSWIN    (1:n) => bus( x(flusolis,1,1)     : )   !DOWNWARD SW FLX AT SFC       (W/m2)
   QLWIN    (1:n) => bus( x(fdsi,1,1)         : )  !DOWNWARD LW FLX AT SFC       (W/m2)
   RT       (1:n) => bus( x(rt,1,1)           : )  !total precip rate       (m/s)
! Other input data
   ZDLAT    (1:n) => bus( x(dlat,1,1)         : )  !Latitude
   ZDLON    (1:n) => bus( x(dlon,1,1)         : )  !Longitude
   ZREFM    (1:n) => bus( x(zusl,1,1)         : )
   ZREFH    (1:n) => bus( x(ztsl,1,1)         : )
   ZGRIDAREA(1:n) => bus( x(gridarea,1,1)     : )         !Grid cell surface area (m2)
   ZLAKEAREA(1:n) => bus( x(lakearea,1,1)     : )         !Lake surface area (km2)
   ZLAKD(1:n)     => bus( x(lakd,1,1)         : )         !Lake depth (average in grid cell) (m)
   ZLAKEFR  (1:n) => bus( x(lakefr,1,1)     : )           !Fraction of grid cell covered by lakes
! Output variables
   QSENS    (1:n) => bus( x(fc,1,indx_sfc)    : )  !Lake average sensible heat flux (inc. ice,snow)
   QEVAP    (1:n) => bus( x(fv,1,indx_sfc)    : )  !Lake average latent heat flux (inc. ice,snow)
   ALVS     (1:n) => bus( x(alvis,1,indx_sfc) : )  !Lake average vis. albedo  (inc. ice,snow)
   QSURF    (1:n) => bus( x(qsurf,1,indx_sfc) : )  !Lake average sat. sp. humidity  (inc. ice,snow)
   zrunofftot(1:n) => bus( x(runofftot,1,indx_sfc) : ) !Lake outflow (kg/m2)

! Sfc layer variables
   z0h      (1:n) => bus( x(z0t,1,indx_sfc)   : )
   z0m      (1:n) => bus( x(z0,1,indx_sfc)    : )
   cmu      (1:n) => bus( x(bm,1,1)           : )
   ctu      (1:n) => bus( x(bt,1,indx_sfc)    : )
   hst_wat  (1:n) => bus( x(hst,1,indx_sfc)   : )
   ilmo_wat (1:n) => bus( x(ilmo,1,indx_sfc)  : )
   zalfaq   (1:n) => bus( x(alfaq,1,1)        : )
   zalfat   (1:n) => bus( x(alfat,1,1)        : )
   zftemp   (1:n) => bus( x(ftemp,1,indx_sfc) : )
   zfvap    (1:n) => bus( x(fvap,1,indx_sfc)  : )
   ztsurf   (1:n) => bus( x(tsurf,1,1)        : )
   ztsrad   (1:n) => bus( x(tsrad,1,1)        : )
   zudiag   (1:n) => bus( x(udiag,1,1)        : )
   zvdiag   (1:n) => bus( x(vdiag,1,1)        : )
   ztdiag   (1:n) => bus( x(tdiag,1,1)        : )
   zqdiag   (1:n) => bus( x(qdiag,1,1)        : )
   zudiagtyp(1:n) => bus( x(udiagtyp,1,indx_sfc) : )
   zvdiagtyp(1:n) => bus( x(vdiagtyp,1,indx_sfc) : )
   ztdiagtyp(1:n) => bus( x(tdiagtyp,1,indx_sfc) : )
   zqdiagtyp(1:n) => bus( x(qdiagtyp,1,indx_sfc) : )
   zfrv     (1:n) => bus( x(frv,1,indx_sfc)   : )
   zsnodp   (1:n) => bus( x(snodp,1,indx_sfc) : )  !Lake snow depth (cm)

!========================================================================================
! PHYSICAL LAKE PROPERTIES FIXED FOR NOW.  (based on small boreal lake L239 at ELA)
! TODO: These should be geophysical input data
      DO I=1,N
!        HLAK(I)=HLAKCON
        HLAK(I)=-ZLAKD(I)        ! Read from geophysical data (in negative meters)
        IF (HLAK(I).LE.0.0) THEN ! Replace zero values with constant value for points where
          HLAK(I)=HLAKCON        ! VF3 is greater than zero but there is no data on depth (i.e. openstreet map data)
        ENDIF
!        LLAK(I)=SQRT(ZLAKEFR(I)*ZGRIDAREA(I)) ! Lake fetch length scale (square-root of lake surface area)
        IF (ZLAKEAREA(I).LT.0.01) THEN
          LLAK(I)=SQRT(ZLAKEFR(I)*ZGRIDAREA(I)) ! Lake fetch length scale (square-root of lake fraction * grid area)
        ELSE
          LLAK(I)=SQRT(ZLAKEAREA(I)*1000000.0) ! Lake fetch length scale (square-root of lake surface area in m2)
        ENDIF
        BLAK(I)=BLAKCON
        NLAK(I)=MIN(NLAKMAX, NINT(HLAK(I)/DELZLK))
        NLAK(I)=MAX(NLAKMIN,NLAK(I))
        RC1 (I)=RC1CON
        RC2 (I)=RC2CON
        RUNFRAC(I)=ZLAKEFR(I)
! fixed water extinction coefficients
        CQ1A(I)=CQ1ACON
        CQ1B(I)=BLAK(I)	!input value
        CQ2A(I)=CQ2ACON
        CQ2B(I)=CQ2BCON
        CQ3A(I)=CQ3ACON
        CQ3B(I)=CQ3BCON
! fixed ice values (from Patterson and Hamblin, 1988, L&O)
        CQ1BI(I)=CQ1BICON
        CQ2BI(I)=CQ2BICON
        CQ3BI(I)=CQ3BICON
      END DO 
!-----------------------------------------------------------------------
!    * SOME INITIAL CONDITIONS (KOUNT=0) BASED IN INITIAL TLAK COMPUTED HERE
!    * EQUATION OF STATE FROM FARMER AND CARMACK (1981,JPO), BUT 
!    * NEGLECTS SALINITY AND PRESSURE EFFECTS
!
      IF (KOUNT .EQ. 0) THEN
        DO I=1,N
!!       DO J=1,NLAK(I) 
!!        TLAK(I,J)=TFREZ+4.0 
!!        IF (J.LE.10) tlak(i,j)=tlak(i,j)+2.0
!!       END DO

!!       LKICEH(I)=0.0     !INITIAL ICE COVER SET TO ZERO
!!       SNICEH(I)=0.0     !INITIAL SNOW ICE  SET TO ZERO
!!       T0(I)=TLAK(I,1)   !INITIAL SKIN TEMP SET TO FIRST LAYER TEMP
!!       TKE(I)=TKEMIN
!!       DELU(I)=0.0
!!       SNO   (I) = 0.0
!!       RHOSNO(I) = 0.0
!!       TSNOW (I) = 0.0
!!       ALBSNO(I) = 0.0
!!       WSNOW (I) = 0.0
         ROFICEH(I)=0.0    !NOT USED  !INITIAL runoff ICE  SET TO ZERO
         NLEV=NLAK(I)
!!       TSED(I)=TLAK(I,NLEV)   !INITIAL SEDIMENT TEMP SET TO LOWEST LAYER TEMP
! ---* 
! ---* INITIAL MIXED LAYER DEPTH ESTIMATED BASED ON INITIAL T PROFILE
! ---* 
         DO 50 J=1,NLAK(I)-1
             JMIX=J
             TTOP=TLAK(I,J)
             TBOT=TLAK(I,J+1)
             IF (TTOP-TBOT .GT. 1.0) EXIT
50       CONTINUE
         HDPTH(I)=DELZLK*JMIX
         DTEMP(I)=TLAK(I,JMIX)-TLAK(I,JMIX+1) !see Spigel et al 1986
! ---* 
! ---* INITIAL MIXED LAYER TEMP (CELSIUS), DENSITY, EXPANSIVITY, 
! ---* AND REDUCED GRAVITY
! ---* 
         TCMIX(I)=( (TLAK(I,1)+TLAK(I,JMIX))/2.0) -TFREZ 
         CALL EQNST(EXPW(I),RHOMIX(I),TCMIX(I),0.0)
         GRED(I) = CGRAV*ABS(RHOW-RHOMIX(I))/RHOW
         IF (GRED(I).LT.0.)  CALL XIT('CLASSL',-1)
        END DO
      ENDIF
   
!----------------------------------------------------------------------------------------
   !# In offline mode the t-step 0 is (correctly) not performed
   !if (FLUVERT.eq.'SURFACE'.and.KOUNT.eq.0) return

!
!========================================================================================
!    * CANADIAN SMALL LAKE MODEL Ver. 2.1
!
!========================================================================================
     RHOIW=RHOICE/RHOW
!----------------------------------------------------------------------------------------
! 0.  Compute cosine solar zenith angle for lake albedo
!----------------------------------------------------------------------------------------
!!!     CALCULATE GREENWICH HOUR
!!!   HZ0 = DATE(5) + float(DATE(6))/360000.0 
!!!   HZ = AMOD ( HZ0+(FLOAT(KOUNT)*DELTAT)/3600. , 24. )
!       Determine the current julian day
!!!   JULIEN = JULIAND( DELTAT, KOUNT, DATE )
!
!!!     Get cosinus of solar angle at LOCAL HOUR
!!!   CALL SUNCOS1(CSZ,SUNOTHER1,SUNOTHER2,SUNOTHER3,SUNOTHER4, N, & 
!!!       bus(x (DLAT,1,1)), BUS(x(DLON,1,1)),  &
!!!        HZ, JULIEN, DATE, .false.)
      ! Calculate greenwich hour 
      call mu_js2ymdhms(jdateo, yy, mo, dd, hh, mn, sec)
    ! hz0 = hh + float(mn)/60. + float(sec)/3600.
      hz0 = hh + mn/60. + sec/3600.
    ! hz = amod(hz0+ (float(kount)*DELT)/3600., 24.)
      hz = amod(hz0+ (kount*DELT)/3600., 24.)
      
      !Determine the current julian day
      julien = real(jdate_day_of_year(jdateo + kount*int(DELT) + MU_JDATE_HALFDAY))
      !Get local solar angle
      call suncos2(CSZ,sunother1,sunother2,sunother3,sunother4,n, &
                   bus(x(dlat,1,1)),bus(x(dlon,1,1)),hz,julien,.false.)

!----------------------------------------------------------------------------------------
! 1.  Compute atmospheric variables; radiation input for TSOLVE
!----------------------------------------------------------------------------------------
      DO I=1,N
       PCPIN(I)=RT(I)*RHOW        !in kg/m2.s
       RRATE(I)=0.0               !only needed for IPCP=4
       SRATE(I)=0.0               !only needed for IPCP=4
       FSSB(I,1) = 0.5*QSWIN(I)   ! used in TSOLVE
       FSSB(I,2) = 0.5*QSWIN(I)
!       IF (QA(I) .LT. EPSILON_CSLM) QA(I)=EPSILON_CSLM ! To avoid a surface spe. humi. of 0.0
      END DO
      CALL CLASSI(VPD,TADP,PADRY,RHOAIR,RHOSNI,R,TR,S,TS,   &
                TA,QA,PCPIN,RRATE,SRATE,PRES, IPCP,N,1,N)

!----------------------------------------------------------------------------------------
! 2.  * Initialize snow and precipitation related variables for beginning of timestep
!----------------------------------------------------------------------------------------
      DO 120 I=1,N
          ZPONDL(I)=0.0
          LKICEH0(I)=LKICEH(I)
!         ICELIM=LLAK(I)*5.45E-6		! Garnaud et al 2022
          ICELIM=(1.909E-4)*(LLAK(I))**0.666667 ! VE Mackay GRL 2023
          FICE(I)=MIN(ICELIM,LKICEH(I))/ICELIM   
          FTOT(I)=1.0
          SNOM1(I)=SNO(I)
          WSNOWM1(I)=WSNOW(I)
          IF(SNO(I).GT.0.0)    THEN
              ZSNOW(I)=SNO(I)/RHOSNO(I)
              IF(ZSNOW(I).GE.(SNOLIM-0.00001)) THEN
                  FLS(I)=1.0
              ELSE
                  FLS(I)=ZSNOW(I)/SNOLIM
                  ZSNOW(I)=SNOLIM
                  WSNOW(I)=WSNOW(I)/FLS(I)
              ENDIF
          ELSE
              ZSNOW(I)=0.0
              FLS(I)=0.0
          ENDIF
          IF(S(I).GT.0.0)                                THEN
              IF(LKICEH(I).GE.DELSKIN)          THEN
                  IF(FLS(I).GT.0.0)         THEN
                      SN(I)=S(I)*FICE(I)/FLS(I)
                  ELSE
                      SN(I)=S(I)
                  ENDIF
                  SL(I)=(1.0-FICE(I))*S(I)
              ELSE
                  SN(I)=0.0
                  SL(I)=S(I)
              ENDIF
          ELSE
              SN(I)=0.0
              SL(I)=0.0
          ENDIF
          IF(R(I).GT.0.0)                                THEN
              IF(FLS(I).GT.0.0)         THEN
                  RN(I)=R(I)
                  RL(I)=(1.0-FLS(I))*R(I)
              ELSE
                  RN(I)=0.0
                  RL(I)=R(I)
              ENDIF
          ELSE
              RN(I)=0.0
              RL(I)=0.0
          ENDIF
          IF(FLS(I).GT.0) THEN
              FSCALL(I)=FLS(I)
              TSTART(I)=TSNOW(I)
          ELSE
              FSCALL(I)=FICE(I)
              TSTART(I)=TFREZ
          ENDIF
          PCPN(I) =FLS(I)*RHOW*RN(I)+FSCALL(I)*RHOSNI(I)*SN(I)
          PCPL(I) =RHOW*RL(I)+RHOSNI(I)*SL(I)
          HTCS(I)=0.0
          HMFN(I)=0.0
          QFN(I)=0.0
          ROFN(I)=0.0
          GZEROSL(I)=0.0
          QMELTS(I)=0.0
          HRPCPL(I)=(TR(I)+TFREZ-T0(I))*RL(I)*HCPW
          HSPCPL(I)=(TS(I)+TFREZ-T0(I))*SL(I)*SPHICE*RHOSNI(I)
          HCONV(I)=CLHMLT*SL(I)*RHOSNI(I)
          HTCL(I)=HRPCPL(I)+HSPCPL(I)-HCONV(I)
!
!         * INITIALIZE OTHER SNOW-RELATED VARIABLES.
!
          SPCN  (I)=0.
          TSPCN (I)=0.
          RPCN  (I)=0.
          TRPCN (I)=0.
          CDMS  (I)=0.
          CDHS  (I)=0.
          DRAGS (I)=0.
          TZEROS(I)=0.
          QZEROS(I)=0.
          ALVSSN(I)=0.
          ALIRSN(I)=0.
          QLWOS (I)=0.
          QSENSS(I)=0.
          QEVAPS(I)=0.
          EVAPS (I)=0.
          HCPSNO(I)=0.
          TCSNOW(I)=0.
          GSNOW (I)=0.
120   CONTINUE

!----------------------------------------------------------------------------------------
! 3.  * Wind speed, precipitation, albedo
!----------------------------------------------------------------------------------------
      DO 140 I=1,N
          VA(I)=MAX(VMIN,SQRT(UWIND(I)*UWIND(I)+VWIND(I)*VWIND(I)))
          FCOR(I)=2.0*7.29E-5*SIN(ZDLAT(I))
          PCP(I)=R(I)+S(I)
140   CONTINUE
! 
!     * DEFINE ICE AND WATER ALBEDOS, USED IN TSOLVL.
!     * STARTING "GROUND ALBEDOS" FOR ROUTINE SNOALBA.
!
      DO 145 I=1,N
        ALBW(I)=0.045/MAX(CSZ(I),0.1)          !std value for water
!
        ALBI(I)=0.08+0.44*(LKICEH(I))**0.28    !thin ice albedo (Vavrus et al 1996)
        ALBI(I)=MIN(ALBI(I),0.44)
!
        IF(LKICEH(I).GT.0.0)          THEN
          ALVSG(I)=ALBI(I)
          ALIRG(I)=ALBI(I)
        ELSE
          ALVSG(I)=0.0
          ALIRG(I)=0.0
        ENDIF
145   CONTINUE

!----------------------------------------------------------------------------------------
! 4.  * SNOW PACK ENERGY AND WATER BALANCE
!     * Routines taken from CLASS Ver. 3.6.1
!----------------------------------------------------------------------------------------
      CALL SNOALBA(ALVSSN,ALIRSN,ALVSSC,ALIRSC,ALBSNO,                &
                   TRSNOWC,ALSNO,TRSNOWG,FSDB,FSFB,RHOSNO,            &
                   REFSN,BCSN,SNO,CSZ,ZSNOW,FLS,ASVDAT,ASIDAT,        &
                   ALVSG,ALIRG,                                       &
                   ILG,IGL,IL1,IL2,N,IALS,NBS,ISNOALB)            
                                                                        
      CALL TLSPREP(GCOEFFS,GCONSTS,CPHCHS,TCSNOW,HCPSNO,IWATER,       &
                   ZRSLDM,ZRSLDH,ZRSLFM,ZRSLFH,ZDSLM,ZDSLH,           &
                   ZOSCLM,ZOSCLH,ZOMLNS,ZOELNS,ZOM,ZOH,               &
                   TVIRTA,TPOTA,CRIBS,DRAGS,CEVAP,IEVAP,ISAND,        &
                   FLS,ZSNOW,TSNOW,RHOSNO,WSNOW,ZREFM,ZREFH,          &
                   ZDIAGM,ZDIAGH,TA,QA,VA,IZREF,ILG,IL1,IL2,IGL,N)    
 
      CALL TSOLVE(ISNOW,FLS,                                          &
                  QSWNS,QLWOS,QTRANSL,QSENSS,QEVAPS,EVAPS,            &
                  TZEROS,QZEROS,GSNOW,QMELTS,CDHS,CDMS,RIBS,CFLUXS,   &
                  FTEMPS,FVAPS,ILMOS,UES,HS,                          &
                  QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,                     &
                  ALVSSN,ALIRSN,CRIBS,CPHCHS,CEVAP,TVIRTA,            &
                  ZOSCLH,ZOSCLM,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,           &
                  GCONSTS,GCOEFFS,TSTART,PCP,TRSNOWG,FSSB,ALSNO,      &
                  THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPONDL,             &
                  IWATER,IEVAP,ITERCT,ISAND,                          &
                  ISLFD,ITG,ILG,IGL,IL1,IL2,N,NBS,ISNOALB,           &
                  TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,                    &
                  DCFLXM,CFLUXM,WZERO,TRTOP,A,B,                      &
                  LZZ0,LZZ0T,FM,FH,ITER,NITER,JEVAP,KF)           
 
      CALL TLSPOST(GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZEROSL,           &
                  TSNBOT,HTCS,HMFN,QFN,EVAPS,RPCN,TRPCN,SPCN,TSPCN,   &
                  GCONSTS,GCOEFFS,T0,ZSNOW,TCSNOW,HCPSNO,QTRANSL,     &
                  RN,TR,SN,TS,TZEROS,RHOSNI,                          &
                  FLS,DELSKIN,ILG,IL1,IL2,N      )
                                                                        
      CALL SNOVAP(RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAPS,QFN,QFL,HTCS,       &
                  WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,                  &
                  FLS,RPCN,SPCN,RHOSNI,WSNOW,ILG,IL1,IL2,N)
 
      CALL TMELT (ZSNOW,TSNOW,QMELTS,RPCN,TRPCN,GZEROSL,RALB,         &
                  HMFN,HTCS,HTC,FLS,HCPSNO,RHOSNO,WSNOW,              &
                  ISAND,IGL,ILG,IL1,IL2,N)
 
      CALL SNINFL(RPCN,TRPCN,ZSNOW,TSNOW,RHOSNO,HCPSNO,WSNOW,         &
                  HTCS,HMFN,PCPL,ROFN,FLS,ILG,IL1,IL2,N)   
 
      CALL SNOALBW(ALBSNO,RHOSNO,ZSNOW,HCPSNO,TSNOW,                  &
                   FLS,SPCN,RALB,WSNOW,RHOMAX,ISAND,                  &
                   ILG,IGL,IL1,IL2,N)                             
                                                                                   
      CALL SNOADD(ALBSNO,TSNOW,RHOSNO,ZSNOW,HCPSNO,HTCS,              &
                  FSCALL,SPCN,TSPCN,RHOSNI,WSNOW,ILG,IL1,IL2,N)
 
!----------------------------------------------------------------------------------------
!
      DO I=1,N
          IF(ZSNOW(I).GT.0.0)                           THEN
              TSNOW(I)=TSNOW(I)+TFREZ
              SNO(I)=ZSNOW(I)*RHOSNO(I)*FSCALL(I)
              WSNOW(I)=WSNOW(I)*FLS(I)
          ELSE
              TSNOW(I)=0.0
              RHOSNO(I)=0.0
              SNO(I)=0.0
              WSNOW(I)=0.0
          ENDIF
!-----------------------------------------------------------
! Melt last bit of snow if below threshold or no ice exists
!
          IF(SNO(I).GT.0.0 .AND. (SNO(I).LT.1.0E-6 .OR. LKICEH(I) &
              .LE.0.01))                                 THEN
              ROFN(I)=ROFN(I)+(SNO(I)+WSNOW(I))/DELT
              PCPL(I)=PCPL(I)+(SNO(I)+WSNOW(I))/DELT
              SL(I)=SL(I)+SNO(I)/(RHOSNO(I)*DELT)
              RL(I)=RL(I)+WSNOW(I)/(RHOW*DELT)
              HTCS(I)=HTCS(I)-TSNOW(I)*(SPHICE*SNO(I)+SPHW*WSNOW(I))/DELT
              TSNOW(I)=0.0
              RHOSNO(I)=0.0
              SNO(I)=0.0
              WSNOW(I)=0.0
          ENDIF
!---------------------------------------------------------------
! Heat capacity calculation for use if necessary by: (1) heat 
! added from snowmelt/rain runoff; and (2) snow ice production
! Any heat added goes into layer 1, not the skin layer.
          ICEBOT=RHOIW*LKICEH0(I)   !bottom of floating ice
          ZTOP=DELSKIN             !top of first layer
          ZBOT=DELSKIN + DELZLK    !bottom of first layer
          IF (ICEBOT .GE. ZBOT) THEN 
            HCAP=HCPICE
          ELSE IF (ICEBOT .LE. ZTOP) THEN
            HCAP=HCPW
          ELSE 
            Z=ICEBOT-ZTOP
            HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK
          ENDIF
!-----------------------------------------------------------
! Add heat of snow melt/rain reaching bottom of snowcover to ice
! Heat added to first layer below skin
! ROFN cools then runs off - ie does not add to LKICEH
!
!-----------------------------------------------------------
! RUNOFF ICE PRODUCTION/HEATING TURNED OFF: JUNE 2016, M.MACKAY
!
          IF (ROFN(I).GT.0.0 .AND. LKICEH0(I).GT.0.0) THEN
             ROFICE=0.0             !runoff ice production turned off
!mdm         ROFICE=DELT*ROFN(I)/RHOICE    !ice production this timestep
!mdm         ROFICEH(I)=ROFICEH(I)+ROFICE  !cumulative diagnostic
! Cool ROFN to TFREZ, then add heat to ice layer below skin
             HWARM=HCPW*(TRPCN(I)-0.0)*DELT*ROFN(I)/RHOW
             HFREZ=RHOICE*CLHMLT*ROFICE
             TLAK(I,1)=TLAK(I,1) + (HFREZ+HWARM)/(HCAP*DELZLK)
!mdm         LKICEH(I)=LKICEH(I)+ROFICE
          ENDIF
!----------------------------------------------------------------------------------------
! SNOW ICE
!---------
! Produce snow-ice if weight of snow cover sufficient
! Pore volume in snow fills with water then freezes
!
          IF(SNO(I).GT.0.0)    THEN
              ROFN(I)=ROFN(I)+SNO(I)/DELT
              HTCS(I)=HTCS(I)-HCPSNO(I)*TSNOW(I)*SNO(I)/(RHOSNO(I)*DELT)
          ENDIF
          BBLFAC=1.0                         !non-bubble fraction in snow ice
          ICECAP =LKICEH0(I)*(RHOW-RHOICE)   !snow holding capacity of ice

           ZBOT = MIN(10.0, (NLAK(I)-1.)*DELZLK)   !limit snow ice production to 10m or depth of lake (less 0.5m)
!          IF (SNO(I).GT.ICECAP .AND. LKICEH0(I).GT.0.0) THEN
           IF (SNO(I).GT.ICECAP .AND. LKICEH0(I).GT.0.0  &
              .AND. LKICEH0(I).LE.ZBOT ) THEN
             ICE0=LKICEH0(I)            !initial ice thickness
             SNO0=SNO(I)                     !initial snow mass
             ZSNOW0=SNO(I)/RHOSNO(I)         !initial mean snow depth 
             PORE=(RHOICE-RHOSNO(I))/RHOICE  !pore volume in snow
             ALPHA=(RHOW*PORE*BBLFAC + RHOICE*(1.0-PORE))/RHOICE
             ETA=(RHOW-RHOICE)/RHOSNO(I)
! Final snow cover mass equals new ice holding capacity
             SNO(I)=RHOSNO(I)*ETA*(ICE0+ALPHA*ZSNOW0)/(1.0+ALPHA*ETA)
             LKICEH(I)=SNO(I)/(RHOSNO(I)*ETA)
! Mass of liquid water that freezes in slush layer
             MASSL=RHOW*PORE*BBLFAC*(SNO0-SNO(I))/RHOSNO(I)
! Cumulative snow-ice diagnostic
             SNICEH(I)=SNICEH(I)+LKICEH(I)-ICE0
! Add heat of fusion to snow pack
! First warm snow to TFREZ, then freeze lake water
! Lake water is assumed to be at TFREZ already (ie comes from layer that
! contains ICEBOT, which is at TFREZ)
            HWARM=HCPSNO(I)*(TFREZ-TSNOW(I))*(SNO0-SNO(I))/RHOSNO(I)
            HFREZ=CLHMLT*MASSL
            HFLX=HFREZ-HWARM
! Latent heat goes goes into snowpack unless snow too thin
           IF (SNO(I) .GE. 5.0E-3) THEN
! Partition latent heat between ice and snow based on heat capacities
!mdm test   HCPRAT=HCPSNO(I)/HCPICE
!mdm test   HFLXS=HFLX*(HCPRAT/(1.0+HCPRAT))
!mdm test   HFLXS=HFLX/2.0
!mdm test   HFLX=0.0   !mdm test
            HFLXS=HFLX           !all latent heat goes into snow cover

            HFLXI=HFLX-HFLXS
            TLAK(I,1)=TLAK(I,1) + HFLXI/(HCAP*DELZLK) !mdm test
            TSNOW(I)=TSNOW(I)+HFLXS/(HCPSNO(I)*SNO(I)/RHOSNO(I)) 
            IF(TSNOW(I).GT.TFREZ) THEN
              QMLTS=(TSNOW(I)-TFREZ)*HCPSNO(I)*SNO(I)/RHOSNO(I) !J/m2
              IF (QMLTS .GT. SNO(I)*CLHMLT) THEN
! All snow melts. Excess heat put in lake
                HFLX=QMLTS-SNO(I)*CLHMLT
                TLAK(I,1)=TLAK(I,1) + HFLX/(HCAP*DELZLK)
                SNO(I)=0.0
                TSNOW(I)=0.0
                RHOSNO(I)=0.0
                WSNOW (I)=0.0
                ZSNOW (I)=0.0
              ELSE
! Some snow melts
                SNO(I)=SNO(I)-QMLTS/CLHMLT
                TSNOW(I)=TFREZ
              ENDIF
            ENDIF                                                  
           ELSE
            TLAK(I,1)=TLAK(I,1) + (HFREZ-HWARM)/(HCAP*DELZLK)
           ENDIF
          ENDIF
          IF(SNO(I).GT.0.0)    THEN
              ROFN(I)=ROFN(I)-SNO(I)/DELT
              HTCS(I)=HTCS(I)+HCPSNO(I)*TSNOW(I)*SNO(I)/(RHOSNO(I)*DELT)
          ENDIF
      ENDDO

!----------------------------------------------------------------------------------------
! 5.  * Drag coefficients and water friction velocity
!       Neutral drag coefficients now computed in NDRCOEFL.F90 and used to estimate
!       roughness lengths for SL_SFCLAYER
!       DRAG is weighted average of snow and bare ice/open water values
!----------------------------------------------------------------------------------------
!!    CALL DRCOEFL (CDML,CDHL,DRAGL,VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,  &
!!                  FLS,ILG,IL1,IL2)
      CALL NDRCOEFL (DRAGL,CDHNL,VA,ZREFM,ZREFH,FICE,FLS,ILG,IL1,IL2)
      DO 150 I=IL1,IL2
          DRAG(I)=FLS(I)*DRAGS(I)+(1.0-FLS(I))*DRAGL(I)
          z0m(I)=ZREFM(I)/EXP(VKC/SQRT(DRAG(I)))
          z0h(I)=z0m(I)/3.0
150   CONTINUE

!     COMPUTE transfer coef. and BOUNDARY LAYER DIAGNOSTICS for water/bare ice (at temp T0)
!
   i = sl_prelim(TA,QA,UWIND,VWIND,PRES,ZREFM,dir_air=vdir,min_wind_speed=VAMIN)
   if (i /= SL_OK) then
      call physeterror('cslm', 'error returned by sl_prelim()')
      return
      !print*, 'Aborting in cslm() because of error returned by sl_prelim()'
      !stop
   endif

   i = sl_sfclayer(th,QA,VA,vdir,ZREFM,ZREFH,T0,QSURF,z0m,z0h,ZDLAT,FCOR, &
        hghtm_diag=zu,hghtt_diag=zt,coefm=cmu,coeft=ctu,flux_t=zftemp, &
        flux_q=zfvap,ilmo=ilmo_wat,ue=zfrv,h=hst_wat,t_diag=ztdiag,    &
        q_diag=zqdiag,u_diag=zudiag,v_diag=zvdiag,tdiaglim=WATER_TDIAGLIM)


   if (i /= SL_OK) then
      call physeterror('cslm', 'error returned by sl_sfclayer()')
      return
      !print*, 'Aborting in cslm() because of error returned by sl_sfclayer()'
      !stop
   endif


   !# Fill surface type-specific diagnostic values
   zqdiagtyp = zqdiag
   ztdiagtyp = ztdiag
   zudiagtyp = zudiag
   zvdiagtyp = zvdiag
   
!
!     COMPUTE water friction velocity, stability corrected drag coef.
      DO 155 I=IL1,IL2
          CDML(I)=zfrv(I)*zfrv(I)/(VA(I)*VA(I))
          CDHL(I)=ctu(I)/VA(I)
          IF (LKICEH(I) .LE. 0) THEN
              USTAR(I)=VA(I)*SQRT(CDML(I)*RHOAIR(I)/RHOW)
          ELSE
              USTAR(I)=0.0            !no wind stress through ice
          ENDIF
155   CONTINUE
 
!----------------------------------------------------------------------------------------
! 6.  * REMOVE STATIC INSTABILITY IF PRESENT BY MIXING WITH ADJACENT LAYER
!     * COMPUTE WATER EXTINCTION COEFFICIENTS
!----------------------------------------------------------------------------------------
!  - include skin layer if no ice is present       (Sept 2, 2011 M MacKay)
      CALL FREECONV(LKICEH,T0,TLAK,RHOIW,NLAK,NLAKMAX,ILG,IL1,IL2)
!
      CALL LKTRANS (CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,BLAK,IL1,IL2,ILG,   &
                    CQ1BI,CQ2BI,CQ3BI)
 
!----------------------------------------------------------------------------------------
! 7.  * COMPUTE SURFACE ENERGY BALANCE FOR CURRENT TIMESTEP AND STEP FORWARD SKIN TEMP T0.
!     * COMPUTE ICE COVER IN SKIN LAYER.
!----------------------------------------------------------------------------------------
      CALL TSOLVL(TLAK,T0,LKICEH,Q0,F0,QSURFL,                   &
                  FSGL,FLGL,HFSL,HEVL,QFL,ALVSL,ALIRL,           &
                  QSWIN,QLWIN,CSZ,TA,QA,VA,PRES,RHOAIR,CDHL,     &
                  GZEROSL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,       &
                  CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,                 &
                  CQ1BI,CQ2BI,CQ3BI,G0,                          &
                  NLAKMAX,ILG,IL1,IL2   )
 
!----------------------------------------------------------------------------------------
! 8.  * COMPUTE SOLAR FLUX (QFLX) FOR CURRENT TIMESTEP
!        include attenuation through ice (including leads) if present
!     * COMPUTE CONDUCTIVE HEAT FLUX (FFLX) FOR CURRENT TIMESTEP
!        thermal conductivity is weighted average of ice and water if ice is present in layer
!----------------------------------------------------------------------------------------
      DO 200 I=IL1,IL2
          ICEBOT=RHOIW*LKICEH(I)
          ICETOP=LKICEH(I)-ICEBOT
          DO 180 J=1,NLAK(I)
              Z=DELSKIN + DELZLK*J 
              IF ( LKICEH(I) .LE. 0.0) THEN        !NO ICE 
                  ATTEN1=CQ1B(I)*Z
                  ATTEN2=CQ2B(I)*Z
                  ATTEN3=CQ3B(I)*Z
              ELSE 
                 IF (ICEBOT .GT. Z) THEN        !Z inside ice
                     ATTEN1=FICE(I)*(CQ1BI(I)*(Z+ICETOP)) + (1.-FICE(I))*CQ1B(I)*Z
                     ATTEN2=FICE(I)*(CQ2BI(I)*(Z+ICETOP)) + (1.-FICE(I))*CQ2B(I)*Z
                     ATTEN3=FICE(I)*(CQ3BI(I)*(Z+ICETOP)) + (1.-FICE(I))*CQ3B(I)*Z
                 ELSE
                  ATTEN1=FICE(I)*(CQ1BI(I)*LKICEH(I)+CQ1B(I)*(Z-ICEBOT))  &
                       + (1.-FICE(I))*CQ1B(I)*Z
                  ATTEN2=FICE(I)*(CQ2BI(I)*LKICEH(I)+CQ2B(I)*(Z-ICEBOT))  &
                       + (1.-FICE(I))*CQ2B(I)*Z
                  ATTEN3=FICE(I)*(CQ3BI(I)*LKICEH(I)+CQ3B(I)*(Z-ICEBOT))  &
                       + (1.-FICE(I))*CQ3B(I)*Z
                 ENDIF
              ENDIF
              QFLX(I,J)=FSGL(I)*(CQ1A(I)*EXP(-ATTEN1) +   &
                              CQ2A(I)*EXP(-ATTEN2) +      &
                              CQ3A(I)*EXP(-ATTEN3) )
180       CONTINUE
200   CONTINUE

!----------------------------------------------------------------------------------------
!
      DO 300 I=IL1,IL2
        ICEBOT=RHOIW*LKICEH(I)
        ICETOP=LKICEH(I)-ICEBOT
        DO 250, J=1,NLAK(I)-1
          ZTOP=DELSKIN + DELZLK*(J -1)
          ZBOT=DELSKIN + DELZLK*J
          IF (ICEBOT .GE. ZBOT) THEN 
            FFLX(I,J)=(-TCICE/DELZLK)*(TLAK(I,J+1)-TLAK(I,J))
          ELSE IF (ICEBOT .LT. ZBOT .AND. ICEBOT .GT. ZTOP) THEN
            TC=((ICEBOT-ZTOP)*TCICE + (ZBOT-ICEBOT)*TCW)/DELZLK
            FFLX(I,J)=(-TC/DELZLK)*(TLAK(I,J+1)-TLAK(I,J))
          ELSE
            FFLX(I,J)=(-TCW/DELZLK)*(TLAK(I,J+1)-TLAK(I,J))
          ENDIF
250     CONTINUE
        NLEV=NLAK(I)
! ---* COMPUTE THERMAL FLUX AT LAKE BOTTOM INTO SEDIMENT
! ---* ASSUMES LAKE DOES NOT FREEZE TO BOTTOM
! ---* (ADIABATIC BC RECOVERED WHEN TCSED=0)
        TSEDT=( (TCSED*TSED(I)/DSED)+(TCW*TLAK(I,NLEV)/DELZLK) )/  &
              ( (TCSED/DSED)+(TCW/DELZLK) )
        FFLX(I,NLEV)= (-2.0*TCW/DELZLK)*(TSEDT-TLAK(I,NLEV))
300   CONTINUE
!
!----------------------------------------------------------------------------------------
! 9.  * COMPUTE TEMPERATURE PROFILE (TLAK) BEFORE TURBULENT MIXING
!        hEat capacity of layer is weighted average of water and ice if necessary
!----------------------------------------------------------------------------------------
!
      DO 400 I=IL1,IL2
        ICEBOT=RHOIW*LKICEH(I)
        ICETOP=LKICEH(I)-ICEBOT
        ICEBOT0=RHOIW*LKICEH0(I)
        ICETOP0=LKICEH0(I)-ICEBOT
!
! --COLUMN  --- TEMP BEFORE MELT/FREEZE
        DO 310, J=1,NLAK(I)
          TLAK0(J)=TLAK(I,J)
          ZTOP=DELSKIN + DELZLK*(J -1)
          ZBOT=DELSKIN + DELZLK*J
          IF (ICEBOT .GE. ZBOT) THEN 
            HCAP=HCPICE
          ELSE IF (ICEBOT .LE. ZTOP) THEN
            HCAP=HCPW
          ELSE 
            Z=ICEBOT-ZTOP
            HCAP=(Z*HCPICE + (DELZLK-Z)*HCPW)/DELZLK
          ENDIF
          IF (J .EQ. 1) THEN
            EFLUX=F0(I)-FFLX(I,1)+Q0(I)-QFLX(I,1)
          ELSE
            EFLUX=FFLX(I,J-1)-FFLX(I,J)+QFLX(I,J-1)-QFLX(I,J)
          ENDIF
          TLAK(I,J)=TLAK(I,J) + (DELT/(DELZLK*HCAP))*EFLUX
310     CONTINUE
! --UPDATE SEDIMENT TEMPERATURE
        NLEV=NLAK(I)
        FFLXSED=(-2.0*TCSED/DSED)*(TSEDB-TSED(I))
        TSED(I)=TSED(I) + (DELT/(DSED*HCPSED))*(FFLX(I,NLEV)-FFLXSED+QSED*QFLX(I,NLEV))
!  
!   ICE GROWTH OR DECAY
! 
        DO 320, J=1,NLAK(I)
          ZTOP=DELSKIN + DELZLK*(J -1)
          ZBOT=DELSKIN + DELZLK*J
          IF (J .EQ. 1) THEN
            EFLUX=F0(I)-FFLX(I,1)+Q0(I)-QFLX(I,1)
          ELSE
            EFLUX=FFLX(I,J-1)-FFLX(I,J)+QFLX(I,J-1)-QFLX(I,J)
          ENDIF
!
! --- FREEZING --------------------
!
          IF (EFLUX .LT. 0.0 .AND. ICEBOT0 .LT. ZBOT .AND. ICEBOT0 .GE. ZTOP) THEN
! ----Net energy flux used to lower T to TFREZ
            IF (TLAK0(J) .GT. TFREZ) THEN
              ECOOL=(ZBOT-ICEBOT0)*HCPW*(TLAK0(J)-TFREZ)/DELT
            ELSE
              ECOOL=0.0
            ENDIF
! ----Remaining energy flux (if any) used to freeze ice
            EAVAIL=EFLUX+ECOOL
            IF (EAVAIL .LT. 0.0) THEN
              NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL
              LKICEH(I)=LKICEH(I)+NEWICE
              ICEBOT=RHOIW*LKICEH(I)
              TLAK(I,J)=TFREZ
!-----LIMIT ICE GROWTH TO THE CURRENT LAYER
              IF (ICEBOT .GT. ZBOT) THEN
                EHEAT=(RHOICE*CLHMLT*(ICEBOT-ZBOT))/DELT
                TLAK(I,J)=TLAK(I,J) - (EHEAT*DELT)/(DELZLK*HCPICE)
                LKICEH(I)=ZBOT/RHOIW
              ENDIF
            ENDIF
          ENDIF
!
! --- MELTING --------------------
!
          IF (EFLUX .GT. 0.0 .AND. ICEBOT0 .GT. ZTOP) THEN
! ----Net energy flux used to raise T to TFREZ
            ZICE=MIN(DELZLK, ICEBOT0-ZTOP)
            EHEAT=ZICE*HCPICE*(TFREZ-TLAK0(J))/DELT
! ----Remaining energy flux (if any) used to melt ice
            EAVAIL=EFLUX-EHEAT
            IF (EAVAIL .GT. 0.0) THEN
              NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL
              LKICEH(I)=LKICEH(I)+NEWICE
              ICEBOT=RHOIW*LKICEH(I)
              TLAK(I,J)=TFREZ
!-----LIMIT ICE MELT TO THE CURRENT LAYER
              IF (ICEBOT .LT. ZTOP) THEN
                EHEAT=RHOICE*CLHMLT*(ZTOP-ICEBOT)/DELT
                TLAK(I,J)=TFREZ + (EHEAT*DELT)/(DELZLK*HCPW)
                LKICEH(I)=ZTOP/RHOIW
              ENDIF
            ENDIF
          ENDIF
320     CONTINUE
400   CONTINUE

!----------------------------------------------------------------------------------------
! 10.  * COMPUTE MIXED LAYER DEPTH (HDPTH), TKE, SHEAR FLOW (DELU)
!      * MIX TEMP OVER MIXED LAYER (now mass weighted)
!----------------------------------------------------------------------------------------
      DO 220 I=IL1,IL2
        IF (TKE(I) .LT. TKEMIN) TKE(I)=TKEMIN
220   CONTINUE
       
      CALL MIXLYR(DTEMP,Q0,NLAK,USTAR,IL1,IL2,ILG,   &
                    HDPTH,TKE,DELU,FQU,BFLX,DISS,EXPW,FSGL,  &
                    FSHEAR,FENTRA,HLAK,LLAK,GRED,TRAN,       &
                    CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,RHOMIX,    &
                    FLGL,HFSL,HEVL,LKICEH)
!
!----------------------------------------------------------------------------------------
!
      DO 600 I=IL1,IL2
       ICEBOT=RHOIW*LKICEH(I)
!
       Z2=DELZLK+DELSKIN
       ZBOT=DELSKIN+DELZLK*NLAK(I)
       IF (HDPTH(I) .LT. Z2) THEN	!fully stratified: no mixing
         JMIX=1
         ZMIX=DELZLK+DELSKIN
       ELSE IF (HDPTH(I) .GE. ZBOT) THEN	!fully mixed
         JMIX=NLAK(I)
         ZMIX=ZBOT
       ELSE
         DO 510, J=2,NLAK(I)
           Z=DELSKIN + DELZLK*J
           IF (HDPTH(I) .LT. Z) EXIT
510      CONTINUE
         JMIX=J-1
         ZMIX=Z-DELZLK
       ENDIF
! ---------------------------------------------------------------------------------------
! -- MIXING UNDER ICE ALLOWED (ie temp mixed between bottom of mixed layer
!                              and bottom of ice cover)
! ---------------------------------------------------------------------------------------
        IF (LKICEH(I) .LE. 0.0) THEN
          TC1=T0(I)-TFREZ
          CALL EQNST(XXX,RHO1,TC1,0.05)
          TBAR=DELSKIN*RHO1*T0(I)
          MASS=DELSKIN*RHO1
        ELSE
          TBAR=0.0
          MASS=0.0
        ENDIF
        DO 520, J=1,JMIX
          ZTOP=DELSKIN + (J-1)*DELZLK
          ZBOT=ZTOP+DELZLK
          IF (ICEBOT .LE. ZTOP) THEN
           TC1=TLAK(I,J)-TFREZ
           CALL EQNST(XXX,RHO1,TC1,ZBOT)
           TBAR=TBAR + RHO1*TLAK(I,J)*DELZLK
           MASS=MASS + RHO1*DELZLK
          ENDIF
520     CONTINUE
        IF ((JMIX .GE. 2) .AND. ((ZMIX-ICEBOT).GE.DELZLK) ) THEN
          TMIX=TBAR/MASS
        ELSE
          TMIX=TLAK(I,1)
        ENDIF  
!-----------------------------------------------------------------------------------------
! MIX TEMPERATURES: include skin layer if no ice 
!-----------------------------------------------------------------------------------------
!--Mix skin with first layer temperature
!mdm    IF ( (LKICEH(I).LE.0.0) .AND. (JMIX.LT.2) ) THEN
        IF ( (LKICEH(I).LE.0.0) .AND. (JMIX.EQ.1) ) THEN
          TC1=T0(I)-TFREZ
          TC2=TLAK(I,1)-TFREZ
          CALL EQNST(XXX,RHO1,TC1,0.05)
          CALL EQNST(XXX,RHO2,TC2,0.5)
          T0(I) = (RHO1*DELSKIN*T0(I) + RHO2*DELZLK*TLAK(I,1))/  &
                   (RHO1*DELSKIN + RHO2*DELZLK)
          TLAK(I,1)=T0(I)
        ELSE IF ( (LKICEH(I).LE.0.0) .AND. (JMIX.GE.2) ) THEN
          T0(I)=TMIX
          DO 525, J=1,JMIX
             TLAK(I,J)=TMIX
525       CONTINUE
!
!--- MIXING UNDER ICE
        ELSE IF ((JMIX .GE. 2) .AND. ((ZMIX-ICEBOT).GE.DELZLK) ) THEN
          DO 530, J=1,JMIX
            ZTOP=DELSKIN + (J-1)*DELZLK
            IF (ICEBOT .LE. ZTOP) THEN
             TLAK(I,J)=TMIX
            ENDIF
530       CONTINUE
        ENDIF
!-----------------------------------------------------------------------
! ---* COMPUTE TEMPERATURE, DENSITY, EXPANSIVITY OF MIXED LAYER WATER
!      EQN OF STATE FROM FARMER AND CARMACK, 1982
!
       TCMIX(I)=TMIX-TFREZ
       JMIXP1=MIN(NLAK(I),JMIX+1)
       CALL EQNST(EXPW(I),RHOMIX(I),TCMIX(I),HDPTH(I))
       CALL EQNST(XXX,RHO1,TLAK(I,JMIXP1)-TFREZ,HDPTH(I))
       GRED(I) = CGRAV*ABS(RHO1-RHOMIX(I))/RHOW
!mdm   GRED(I) = GRAV*ABS(RHO1-RHOMIX(I))/RHO1
       IF (GRED(I).LT.0.) 		CALL XIT('CLASSL',-2)
!
!-----------------------------------------------------------------------
! ---* MIX TEMPERATURE IN THERMOCLINE LAYER 
!          DISABLED FOR NOW: JAN 2013 (M.MACKAY)
!
!--Compute thermocline layer thickness (DELTHRM)
!--Limit allowed values 
!
      IF (LKICEH(I) .LE. 0.0) THEN 
       IF (GRED(I) .GT. 0.0) THEN
        DELTHRM=0.3*DELU(I)*DELU(I)/GRED(I)
       ELSE
        DELTHRM=0.0
       ENDIF
       IF (DELTHRM.GT.DELMAX) DELTHRM=DELMAX 
       IF (DELTHRM.LT.DELMIN) DELTHRM=DELMIN 

       HTOP=HDPTH(I)-(DELTHRM/2.)
       HBOT=HDPTH(I)+(DELTHRM/2.)
       DO 540, J=1,NLAK(I)
         Z=DELZLK*J
         IF (HTOP .LT. Z) EXIT
540    CONTINUE
       JTOP=J-1
       IF (JTOP.LT.1) JTOP=1
       ZTOP=DELZLK*JTOP
       TTOP=TLAK(I,JTOP)

       DO 550, J=1,NLAK(I)
         Z=DELZLK*J
         IF (HBOT .LT. Z) EXIT
550    CONTINUE
       JBOT=J-1
       IF (JBOT.LT.1) JBOT=1
       ZBOT=DELZLK*JBOT
       TBOT=TLAK(I,JBOT)
       IF (JBOT.LT.JTOP) 		CALL XIT('CLASSL',-3)

!-DISABLED
!-       IF (JBOT.GT.JTOP) THEN
!-        DO 560, J=JTOP,JBOT
!-          Z=DELZLK*J
!-          TLAK(I,J)=TTOP+((Z-ZTOP)*(TBOT-TTOP)/(ZBOT-ZTOP))
!-560     CONTINUE
!-       ENDIF

!--Compute temperature jump beneath mixed layer
!--
       IF (JMIX .LT. NLAK(I)) THEN
          DTEMP(I)=TMIX-TLAK(I,JBOT) 	!see Spigel et al 1986
       ELSE 
          DTEMP(I)=0.
       ENDIF
      ELSE 
        DTEMP(I)=0.
        DELTHRM=0.
      ENDIF
600   CONTINUE

!----------------------------------------------------------------------------------------
!--Compute Wedderburn number and other diagnostics (not used for now)
!
      DO 605 I=IL1,IL2
       IF  (LKICEH(I) .LE. 0.0 .AND. USTAR(I) .GT. 0.0 ) THEN
         WEDB=GRED(I)*HDPTH(I)*HDPTH(I)/(LLAK(I)*USTAR(I)*USTAR(I))
       ELSE
         WEDB=-999    !because USTAR=0
       ENDIF
!--Compute heat of melting/freezing for energy balance diagnostic
       HMFL(I) = (LKICEH(I)-LKICEH0(I))*RHOICE*CLHMLT/DELT
       OSW(I) = QSWIN(I) - FLS(I)*(QSWNS(I)-QTRANSL(I)) - FSGL(I)
! --RESET SNICEH IF ICE HAS COMPLETELY MELTED; ENSURE CONSISTENCY BETWEEN
!   FICE AND LKICEH
       IF  (LKICEH(I) .LE. 0.0 ) SNICEH(I)=0.0
!      ICELIM=LLAK(I)*5.45E-6		! Garnaud et al 2022
       ICELIM=(1.909E-4)*(LLAK(I))**0.666667 ! VE Mackay GRL 2023
       FICE(I)=MIN(ICELIM,LKICEH(I))/ICELIM

605   CONTINUE

!----------------------------------------------------------------------------------------
! 11.  * SCREEN LEVEL DIAGNOSTICS AND FINAL FLUXES
!----------------------------------------------------------------------------------------
      DO I=IL1,IL2
          CDM(I)=FLS(I)*CDMS(I)+(1.0-FLS(I))*CDML(I)
          CDH(I)=FLS(I)*CDHS(I)+(1.0-FLS(I))*CDHL(I)
          ZOM(I)=ZREFM(I)/EXP(VKC/SQRT(DRAG(I)))
          ZOH(I)=ZOM(I)/3.0
          TSURF(I)=FLS(I)*TZEROS(I)+(1.0-FLS(I))*T0(I)
          QSURF(I)=FLS(I)*QZEROS(I)+(1.0-FLS(I))*QSURFL(I)
          ALVS(I)=FLS(I)*ALVSSN(I)+(1.0-FLS(I))*ALVSL(I)
          ALIR(I)=FLS(I)*ALIRSN(I)+(1.0-FLS(I))*ALIRL(I)
      ENDDO
!
      DO I=IL1,IL2                                            
          FACTM=ZDIAGM(I)+ZOM(I)  
          FACTH=ZDIAGH(I)+ZOM(I)   
          RATIOM=SQRT(CDM(I))*LOG(FACTM/ZOM(I))/VKC              
          RATIOM=MIN(RATIOM,1.)                                   
          RATIOH=SQRT(CDM(I))*LOG(FACTH/ZOH(I))/VKC              
          RATIOH=MIN(RATIOH,1.)                                   
          IF(TSURF(I).GT.TA(I))                THEN
              RATIOH=RATIOH*CDH(I)/CDM(I)                         
              RATIOH=MIN(RATIOH,(FACTH/ZREFH(I))**(1./3.))         
          ENDIF                                                   
          ST(I)=TSURF(I)-(MIN(RATIOH,1.))*(TSURF(I)-TA(I))       
          SQ(I)=QSURF(I)-(MIN(RATIOH,1.))*(QSURF(I)-QA(I))       
          SU(I)=RATIOM*UWIND(I)                                  
          SV(I)=RATIOM*VWIND(I)                                  
      ENDDO
!                                                                       
!  * Total fluxes are weighted averages of snow-covered and non-snow covered fluxes
      DO I=IL1,IL2                                              
          FSGS(I) =FLS(I)*(QSWNS(I)-QTRANSL(I))           
          FLGS(I) =FLS(I)*(QLWIN(I)-QLWOS(I))            
          HFSS(I) =FLS(I)*QSENSS(I)                     
          HEVS(I) =FLS(I)*QEVAPS(I)                     
          HTCS(I) =HTCS(I)-(FLS(I)*GZEROSL(I))

          HTCL(I) =HTCL(I)+(FLS(I)*GZEROSL(I))
          QLWAVG(I)=FLS(I)*QLWOS(I)+(1.0-FLS(I))*QLWIN(I)-FLGL(I)
          QSENS(I)=HFSL(I)+HFSS(I)
          TFLUX(I)=-QSENS(I)/(RHOAIR(I)*SPHAIR)
          QEVAP(I)=HEVL(I)+HEVS(I)
          EVAP(I)=QFL(I)+QFN(I)
	  zevlak(I) = zevlak(I) + QFL(I)/adj_cslm_liquevap * DELT + QFN(I) * DELT
          QFLUX(I)=-EVAP(I)/RHOAIR(I)
          EVPPOT(I)=EVAP(I)
          EVAPB(I)=1.0
! The lake average surface temp is GT, not T0.  T0 does not include temp of snow covered part of lake.
          GT(I)=FLS(I)*TZEROS(I)+(1.0-FLS(I))*T0(I)
      ENDDO
!                                                                       
!  * FINAL FLUXES FOR GEM
!VDIR NODEP
      DO I=IL1,IL2    
          ZTSURF (I) = GT(I)
          ZTSRAD (I) = GT(I)
          ZALFAT (I) = TFLUX(I)
          ZALFAQ (I) = QFLUX(I)
          IF (.not. IMPFLX) CTU(I) = 0.0
          IF (IMPFLX) then
            ZALFAT (I) = - CTU(I) * GT(I)
            ZALFAQ (I) = - CTU(I) * QSURF(I)
          ENDIF
      ENDDO
!                                                                       
!  * WATER BALANCE
      DO I=IL1,IL2    
          IF (HLAKSIL(I) .GT. 0.0) THEN
!           ZRUNOFFTOT(I) = DELT*RC1(I)*(HLAKSIL(I))**RC2(I)  !Rating curve for lake outflow (kg/m2/dt)

! EG_MOD 
!           ZRUNOFFTOT(I) = (DELT*RHOW/(ZLAKEFR(I)*ZGRIDAREA(I)))* &
!                 0.001*LLAK(I)*RC1(I)*(HLAKSIL(I)**RC2(I))

! check if the grid-cell belongs to a big lake by comparing the area of the grid cell
! to the actual full area of the lake.
	    IF ( (ZLAKEAREA(I)*1000000.0) .GT. ZGRIDAREA(I) ) THEN
	       ! Grid-cell belongs to a large lake spanning over multiple grid cells;
	       ! don't slowdown lake outflow, let the routing scheme handle this.
                ZRUNOFFTOT(I) = HLAKSIL(I) * RHOW
            ELSE
	       ! Rating curve for rectangular weir (CMS) converted to kg/m^2
                ZRUNOFFTOT(I) = (DELT*RHOW/(ZLAKEFR(I)*ZGRIDAREA(I)))* &
                 0.001*SQRT(ZLAKEFR(I)*ZGRIDAREA(I))*RC1(I)*(HLAKSIL(I)**RC2(I))
            ENDIF

! END_EG_MOD

          ELSE
           ZRUNOFFTOT(I) = 0.0
          ENDIF

          HLAKSILM1(I)=HLAKSIL(I) ! from previous timestep

!         Includes the weighted runoff from land elements (soil, urban, glaciers)   
! EG_MOD2:
! do not remove evap. from liquid water from the lake in CSLM, but put it into lake outlet runoff
! to still close the water balance, and let the routing scheme remove the evaporation from lakes.
!          HLAKSIL(I) = HLAKSIL(I) + (DELT*(PCPIN(I)-PCPN(I)-EVAP(I)+QFN(I)+ROFN(I)) &
!                       + ZROFINLAK(I) - ZRUNOFFTOT(I))/RHOW
	  HLAKSIL(I) = HLAKSIL(I) + (DELT*(PCPIN(I)-PCPN(I)+ROFN(I)) &
                       + ZROFINLAK(I) - ZRUNOFFTOT(I))/RHOW
! END_EG_MOD2

          ! Variables outputted to check the water balance of the lake tile
          ! All are in kg/m2
          LFXI(I) = LFXI(I)           + DELT*PCPIN(I) + ZROFINLAK(I)

! EG_MOD3: in order to close the water balance due to not removing evap. from lake, change 
! the way LFXO is computed:
!          LFXO(I) = LFXO(I)           + DELT*EVAP(I)  + ZRUNOFFTOT(I)
	 LFXO(I) = LFXO(I)           + DELT*QFN(I)  + ZRUNOFFTOT(I)
! END_EG_MOD3

          LSTD(I) = RHOW*HLAKSILM1(I) + SNOM1(I)      + WSNOWM1(I)
          LSTF(I) = RHOW*HLAKSIL(I)   + SNO(I)        + WSNOW(I)

!          IF (I .EQ. 59 .AND. TRNCH .EQ. 59) THEN
!            print *, "----------------------", KOUNT
!            print *, "FLUX IN  =", LFXI(I)
!            print *, "FLUX OUT =", LFXO(I)
!            print *, "STOR DEB =", LSTD(I)
!            print *, "STOR FIN =", LSTF(I)
!          ENDIF

!         Include CSLM simulated snow depth in snodpl(indx_lake)
          ZSNODP(I)  = ZSNOW(I)

! EG_MOD4: 
! Remove liquid water evap. from CSLM output runoff in order
! to let Watroute remove the evaporation from lakes (and prevent CSLM lakes from drying out).
! and tune evaporation from liquid water by dividing it by constant "adj_cslm_liquevap":
          ZRUNOFFTOT(I) = ZRUNOFFTOT(I) - (QFL(I)*DELT/adj_cslm_liquevap)
! END_EG_MOD4

      ENDDO
!========================================================================================
!  * FILL THE ARRAYS TO BE AGGREGATED LATER IN S/R AGREGE
   call FILLAGG ( BUS, BUSSIZ, PTSURF, PTSURFSIZ, INDX_LAKE, SURFLEN )

   return

 end subroutine cslm_main
