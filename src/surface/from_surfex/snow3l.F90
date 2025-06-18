!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE SNOW3L(HSNOWRES, OMEB, HIMPLICIT_WIND,                           &
                  PPEW_A_COEF, PPEW_B_COEF,                                 &
                  PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                  PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                  PSNOWAGE, PTSTEP, PPS, PSR, PUNLOAD, PRR, PPSN3L,PRESA_SV,&
                  PTAR,PTAC,PTG,PSW_RAD,PQA,PVMOD,PWIND_DRIFT,              &
                  PLW_RAD, PRHOA,                                           &
                  PUREF,PEXNS,PEXNA,PDIRCOSZW,                              &
                  PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
                  PSOILCOND,PD_G,PLVTT,PLSTT,                               &
                  PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
                  PTHRUFAL,PSNOWMELT,PSNREFREEZ,                            &
                  PGRNDFLUX,PEVAPCOR,PSOILCOR,                              &
                  PGFLXCOR,PSNOWSFCH,PDELHEATN,PDELHEATN_SFC,               &
                  PDELPHASEN, PDELPHASEN_SFC,                               &
                  PSWNETSNOW,PSWNETSNOWS,PLWNETSNOW,PRESTOREN,              &
                  PRNSNOW,PHSNOW,PGFLUXSNOW,                                &
                  PHPSNOW,PLES3L,PLEL3L,PEVAP,PSNDRIFT,PRI,                 &
                  PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,PQS,           &
                  PPERMSNOWFRAC,PFORESTFRAC,PZENITH,                        &
                  HSNOWDRIFT,OSNOWDRIFT_SUBLIM                              )
!     ##########################################################################
!
!!****  *SNOW3L*
!!
!!    PURPOSE
!!    -------
!
!     3(N)-Layer snow scheme option (Boone and Etchevers, J Hydrometeor., 2000)
!
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    ISBA-ES: Boone and Etchevers (2001)
!!    ISBA: Belair (1995)
!!    ISBA: Noilhan and Planton (1989)
!!    ISBA: Noilhan and Mahfouf (1996)
!!
!!    AUTHOR
!!    ------
!!      A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    7/99
!!      Modified by A.Boone 05/02 (code, not physics)
!!      Modified by A.Boone 11/04 i) maximum density limit imposed (although
!!                                rarely if ever reached), ii) check to
!!                                see if upermost layer completely sublimates
!!                                during a timestep (as snowpack becomes vanishly
!!                                thin), iii) impose maximum grain size limit
!!                                in radiation transmission computation.
!!
!!      Modified by B. Decharme  (03/2009): Consistency with Arpege permanent
!!                                          snow/ice treatment (LGLACIER for alb)
!!      Modified by A. Boone     (04/2010): Implicit coupling and replace Qsat and DQsat
!!                                          by Qsati and DQsati, respectively.
!!      Modified by E. Brun      (08/2012): Mass conservation in SNOW3LEVAPGONE
!!      Modified by B. Decharme  (08/2012): Loop optimization
!!      Modified by B. Decharme  (09/2012): New wind implicitation
!!      Modified by E. Brun      (10/2012): Bug in vertical snow redistribution
!!      Modified by B. Decharme  (08/2013): Qsat as argument (needed for coupling with atm)
!!                                          Add snow drift effect like in CROCUS
!!                                          Change snow albedo (CROCUS spectral albedo)
!!                                          Change snow viscosity like in CROCUS
!!      Modified by P. Samuelsson(10/2014): MEB added an optional snow albedo calculation
!!                                          for litter
!!      Modified by A. Boone     (10/2014): MEB modifs: removed above since spectral albedo 
!!                                          method from CROCUS added...can eventually add above
!!                                          back in a manner consistent with spectral bands...
!!      Modified by A. Boone     (10/2014): MEB modifs: permit option to impose fluxes at sfc  
!!      Modified by A. Boone     (10/2014): SNOW3LREFRZ and SNOW3LEVAPN edited to give
!!                                          better enthalpy conservation.
!!      Modified by B. Decharme  (03/2016): No snowdrift under forest
!!      Modified by B. Decharme  (03/2020): Diagnostics
!!      Modified by A. Boone     (09/2022): Consistent with modifs for CROCUS by A. Bouchet,
!!                                          add unloading rate from explicit canopy (MEB) 
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT, XCL, XDAY
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_SNOW_METAMO, ONLY : XSNOWDZMIN
USE MODD_SNOW_PAR,    ONLY : XSNOWDMIN, NSPEC_BAND_SNOW
!
USE MODE_SNOW3L
USE MODE_THERMOS,  ONLY : QSATI
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration

!
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES
!                                      HSNOWRES  = ISBA-SNOW3L turbulant exchange option
!                                      'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                                      'RIL' = Limit Richarson number under very stable
!                                              conditions (currently testing)
!                                      'RI1' = Limit Richarson number under very stable
!                                              conditions to 0.1
!                                      'RI2' = Limit Richarson number under very stable
!                                              conditions to 0.026
!
LOGICAL, INTENT(IN)                 :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
!                                                 ! as an upper boundary condition to the explicit snow schemes. 
!                                                 ! If = False, then energy
!                                                 ! budget and fluxes are computed herein.
!
CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)    :: PPS, PTAR, PTAC, PSW_RAD, PQA,                       &
                                     PVMOD, PWIND_DRIFT, PLW_RAD, PSR, PRR, PUNLOAD  
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = snow rate (SWE) [kg/(m2 s)]
!                                      PUNLOAD = snow unloading rate [kg/(m2 s)]
!                                      PTAC    = atmospheric temperature at level za (K)
!                                                NOTE, when MEB used, it corresponds to the canopy air T
!                                      PTAR    = atmospheric temperature above surface & canopy (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PWIND_DRIFT   = modulus of the wind under canopy if PVMOD is above canopy, otherwise PWIND_DRIFT = PVMOD (m/s)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
!
REAL, DIMENSION(:), INTENT(IN)    :: PSOILCOND, PD_G, PPSN3L
!                                      PSOILCOND = soil thermal conductivity [W/(m K)]
!                                      PD_G      = Assumed first soil layer thickness (m)
!                                                  Used to calculate ground/snow heat flux
!                                      PPSN3L    = snow fraction
REAL, DIMENSION(:)                :: PRESA_SV     
!                                 Aerodynamic  resistance computed externally (0 if not used)
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PRHOA, PZ0, PZ0EFF, &
                                       PALB, PZ0H, PPERMSNOWFRAC, PFORESTFRAC 
!                                      PZ0EFF    = roughness length for momentum
!                                      PZ0       = grid box average roughness length
!                                      PZ0H      = grid box average roughness length for heat
!                                      PZREF     = reference height of the first
!                                                  atmospheric level
!                                      PUREF     = reference height of the wind
!                                      PRHOA     = air density
!                                      PEXNS     = Exner function at surface
!                                      PEXNA     = Exner function at lowest atmos level
!                                      PDIRCOSZW = Cosinus of the angle between the
!                                                  normal to the surface and the vertical
!                                      PALB      = soil/vegetation albedo
!                                      PPERMSNOWFRAC  = fraction of permanet snow/ice
!                                      PFORESTFRAC = fraction of forest
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                         PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)    :: PTG
!                                      PTG       = Surface soil temperature (effective
!                                                  temperature the of layer lying below snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLVTT, PLSTT ! = latent heats for hydrology
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWALB
!                                      PSNOWALB = Prognostic surface snow albedo
!                                                 (does not include anything but
!                                                 the actual snow cover)
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PSNOWHEAT, PSNOWRHO, PSNOWSWE
!                                      PSNOWHEAT = Snow layer(s) heat content (J/m2)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
!

REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE  ! Snow grain age
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PRNSNOW, PHSNOW, PLES3L, PLEL3L, &
                                       PHPSNOW, PEVAP,  PGRNDFLUX, PEMISNOW
!                                      PLES3L      = evaporation heat flux from snow (W/m2)
!                                      PLEL3L      = sublimation (W/m2)
!                                      PHPSNOW     = heat release from rainfall (W/m2)
!                                      PRNSNOW     = net radiative flux from snow (W/m2)
!                                      PHSNOW      = sensible heat flux from snow (W/m2)
!                                      PEVAP       = total evaporative flux (kg/m2/s)
!                                      PGRNDFLUX   = soil/snow interface heat flux (W/m2)
!                                      PEMISNOW    = snow surface emissivity
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGFLUXSNOW
!                                      PGFLUXSNOW  = net heat flux from snow (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
!                                      PSWNETSNOW = net shortwave radiation entering top of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                      PSWNETSNOWS= net shortwave radiation in uppermost layer of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                                   Used for surface energy budget diagnostics
!                                      PLWNETSNOW = net longwave radiation entering top of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PUSTAR, PCDSNOW, PCHSNOW, PRI
!                                      PCDSNOW    = drag coefficient for momentum over snow (-)
!                                      PUSTAR     = friction velocity over snow (m/s)
!                                      PCHSNOW    = drag coefficient for heat over snow (-)
!                                      PRI        = Richardson number (-)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWLIQ, PSNOWDZ
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)
!                                      PSNOWTEMP = Snow layer(s) temperature (m)
!                                      PSNOWDZ   = Snow layer(s) thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL, PEVAPCOR, PSOILCOR, PGFLXCOR, &
                                       PRESTOREN, PSNOWSFCH, PDELHEATN, PDELHEATN_SFC, &
                                       PDELPHASEN, PDELPHASEN_SFC, PSNOWMELT, PSNREFREEZ
!                                      PTHRUFAL  = rate that liquid water leaves snowpack :
!                                                  paritioned into soil infiltration/runoff
!                                                  by ISBA [kg/(m2 s)]
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!                                      PSOILCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!                                      PGFLXCOR  = flux correction to underlying soil for vanishing snowpack
!                                                  (to put any energy excess from snow to soil) (W/m2)
!                                      PRESTOREN = heat flux between the surface and sub-surface 
!                                                  snow layers (W/m2)
!                                      PSNOWSFCH = snow surface layer pseudo-heating term owing to
!                                                  changes in grid thickness            (W m-2)
!                                      PDELHEATN = total snow heat content change  (W m-2)
!                                      PDELHEATN_SFC = surface layer snow heat content change during the timestep (W m-2)
!                                      PDELPHASEN = latent heating due to snow melt/freeze  (W m-2)
!                                      PDELPHASEN_SFC = latent heating due to surface layer snow melt/freeze  (W m-2)
!                                      PSNOWMELT = snowmelt in the snowpack  [kg/(m2 s)]
!                                      PSNREFREEZ = refreezing of water in the snowpack  [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNDRIFT
!                                      PSNDRIFT    = blowing snow sublimation (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT)   ::   PSNOWHMASS
!                                      PSNOWHMASS  = heat content change due to mass
!                                                    changes in snowpack (J/m2): for budget
!                                                    calculations only.
!
REAL, DIMENSION(:), INTENT(OUT)   :: PQS
!                                    PQS = surface humidity
!
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle
!
!
CHARACTER(4), INTENT(IN)          :: HSNOWDRIFT  ! Snowdrift scheme :
                                                 ! 'NONE': No snowdrift scheme
                                                 !  'DFLT':  Snowdrift scheme activated
                                                 !  Other options are available in Crocus

LOGICAL, INTENT(IN)               ::  OSNOWDRIFT_SUBLIM ! activate snowdrift, sublimation during drift
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                     :: ZEPSILON = 1.0E-12
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWTEMP, ZSCAP, ZSNOWDZN, ZSCOND,    &
                                                      ZRADSINK, ZWORK2D, ZSNOWTEMP0, ZSNOWLIQ0 
!                                      ZSNOWTEMP  = Snow layer(s) averaged temperature (K)
!                                      ZSCAP      = Snow layer(s) heat capacity [J/(K m3)]
!                                      ZSNOWDZN   = Updated snow layer thicknesses (m)
!                                      ZSCOND     = Snow layer(s) thermal conducivity [W/(m K)]
!                                      ZRADSINK   = Snow solar Radiation source terms (W/m2)
!                                      ZWORK2D    = working variable (*)
!                                      ZSNOWTEMP0 = value of snow layer temperature before time integration (K)
!                                      ZSNOWLIQ0  = value of liquid water before time integration (m)
!
REAL, DIMENSION(SIZE(PTAR))         :: ZSNOW, ZSFCFRZ, ZTSTERM1, ZTSTERM2,                   &
                                       ZCT, ZRA, ZSNOWTEMP01  
!                                      ZSNOW      = Total snow depth (m)
!                                      ZCT        = inverse of the product of snow heat capacity
!                                                   and layer thickness [(m2 K)/J]
!                                      ZRA        = Surface aerodynamic resistance
!                                      ZTSTERM1,ZTSTERM2 = Surface energy budget coefficients
!                                      ZSNOWTEMP01= value of uppermost snow temperature
!                                                   before time integration (K)
!
LOGICAL, DIMENSION(SIZE(PTAR))      :: GSFCMELT
!                                      GSFCMELT   = FLAG if surface melt is occurring, used
!                                                   for surface albedo calculation.
!
REAL, DIMENSION(SIZE(PTAR))         :: ZRSRA, ZDQSAT, ZQSAT, ZRADXS, ZMELTXS, ZLIQHEATXS, &
                                       ZLWUPSNOW, ZGRNDFLUX, ZGRNDFLUXO, ZGRNDFLUXI,      &
                                       ZPSN3L, ZSNOWHMASS1, ZSNOWMASSNEW, ZDELPHASE,       &
                                       ZDELPHASE_SFC, ZSNOWGONE_DELTA
!                                      ZRSRA    = air density over aerodynamic resistance
!                                      ZDQSAT   = derrivative of saturation specific humidity
!                                      ZQSAT    = saturation specific humidity
!                                      ZRADXS   = shortwave radiation absorbed by soil surface
!                                                 (for thin snow sover) (W m-2)
!                                      ZMELTXS  = excess energy for snowmelt for vanishingly
!                                                 thin snowpacks: used to heat underlying surface
!                                                 in order to maintain energy conservation (W m-2)
!                                      ZLIQHEATXS = excess snowpack heating for vanishingly thin
!                                                 snow cover: add energy to snow/ground heat
!                                                 flux (W m-2)
!                                      ZLWUPSNOW = upwelling longwave radiative flux (W m-2)
!                                      ZGRNDFLUXO= snow-ground flux before correction (W m-2)
!                                      ZGRNDFLUX = snow-ground flux after correction (W m-2)
!                                      ZGRNDFLUXI= for the case where the ground flux is imposed,
!                                                  this is the actual imposed value.
!                                      ZPSN3L    = snow fraction: different use if MEB "on".
!                                                  In this case, it is only used for Tg update
!                                                  since only this variable has a sub-grid relevance.
!                                      ZSNOWMASSNEW = combined mass flux downward at surface from snowfall and 
!                                                     possible unloading from the explicit canopy (kg m-2 s-1)
!                                      ZDELPHASE = Total snow liquid phase change (kg/m2/s)
!                                      ZDELPHASE_SFC= uppermost snow layer liquid phase change (kg/m2/s)
!                                      ZSNOWGONE_DELTA = snow melts during a time step (-)
!
REAL, DIMENSION(SIZE(PTAR))         :: ZUSTAR2_IC, ZTA_IC, ZQA_IC, ZWORK, ZWORK2, ZWORK3,                 &
                                       ZPET_A_COEF_T, ZPEQ_A_COEF_T, ZPET_B_COEF_T, ZPEQ_B_COEF_T  
!                                      ZUSTAR2_IC    = implicit lowest atmospheric level friction (m2/s2)
!                                      ZTA_IC        = implicit lowest atmospheric level air temperature
!                                      ZQA_IC        = implicit lowest atmospheric level specific humidity
!                                      ZPET_A_COEF_T = transformed A-air temperature coefficient
!                                      ZPET_B_COEF_T = transformed B-air temperature coefficient
!                                      ZPEQ_A_COEF_T = transformed A-air specific humidity coefficient
!                                      ZPEQ_B_COEF_T = transformed B-air specific humidity coefficient
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),NSPEC_BAND_SNOW)  :: ZSPECTRALALBEDO, ZSPECTRALWORK
!                                                     ZSPECTRALALBEDO=spectral albedo 
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEAT0
!
INTEGER                            :: JJ, JI     ! Loop control
!
INTEGER                            :: INI        ! number of point
INTEGER                            :: INLVLS     ! number of snow layers
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! - - ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3L',0,ZHOOK_HANDLE)
!
!       0.     Initialization
!               --------------
!
PGFLUXSNOW    (:) = 0.0
PTHRUFAL      (:) = 0.0 
PEVAPCOR      (:) = 0.0
PSOILCOR      (:) = 0.0
PGFLXCOR      (:) = 0.0
PRESTOREN     (:) = 0.0
PSNOWSFCH     (:) = 0.0
PDELHEATN     (:) = 0.0
PDELHEATN_SFC (:) = 0.0
PDELPHASEN    (:) = 0.0
PDELPHASEN_SFC(:) = 0.0
PSNDRIFT      (:) = 0.0
PSNOWHMASS    (:) = 0.0
PSNOWMELT     (:) = 0.0
PSNREFREEZ    (:) = 0.0
!   
! NOTE that snow layer thickness is used throughout this code: SWE
! is only used to diagnose the thickness at the beginning of this routine
! and it is updated at the end of this routine.
!
PSNOWDZ(:,:) = PSNOWSWE(:,:)/PSNOWRHO(:,:)
!
INI          = SIZE(PSNOWSWE(:,:),1)
INLVLS       = SIZE(PSNOWSWE(:,:),2)    ! total snow layers
!
ZUSTAR2_IC = 0.0
ZTA_IC     = 0.0
ZQA_IC     = 0.0
!
ZGRNDFLUXO = 0.0
ZGRNDFLUX  = PGRNDFLUX
!
ZMELTXS     (:) = 0.0
ZLIQHEATXS  (:) = 0.0
ZRADXS      (:) = 0.0
!
WHERE(PSNOWHEAT(:,:)==XUNDEF)
     ZSNOWHEAT0(:,:) = 0.0
ELSEWHERE
     ZSNOWHEAT0(:,:) = PSNOWHEAT(:,:) ! save initial heat content
ENDWHERE
!
ZSNOWTEMP0(:,:) = PSNOWTEMP(:,:) ! for MEB, this is the updated T profile
!
!*       1.     Snow total depth
!               ----------------
!
ZSNOW(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSNOW(JI) = ZSNOW(JI) + PSNOWDZ(JI,JJ)  ! m
   ENDDO
END DO
!
ZWORK(:)=ZSNOW(:)
!
! Calculate new snow albedo at time t
!
ZWORK2(:)=PSNOWALB(:)
!
CALL SNOW3LALB(ZWORK2,ZSPECTRALALBEDO,PSNOWRHO(:,1),PSNOWAGE(:,1),PPERMSNOWFRAC,PPS)
ZWORK3(:) = PSNOWALB(:)/ZWORK2(:)
DO JJ=1,SIZE(ZSPECTRALALBEDO,2)
   DO JI=1,INI
      ZSPECTRALALBEDO(JI,JJ)=ZSPECTRALALBEDO(JI,JJ)*ZWORK3(JI)
   ENDDO
ENDDO
!
!*       2.     Snowfall
!               --------
!
! Calculate uppermost density and thickness changes due to snowfall,
! and add heat content of falling snow
!
 CALL SNOW3LFALL(PTSTEP, PSR, PTAR, PWIND_DRIFT, ZSNOW, PSNOWRHO, PSNOWDZ, &
                PSNOWHEAT, PSNOWHMASS, ZSNOWHMASS1, PSNOWAGE,        &
                PPERMSNOWFRAC, PUNLOAD, ZSNOWMASSNEW)
!
! Calculate new snow albedo at time t if snowfall
!
CALL SNOW3LALB(ZWORK2,ZSPECTRALWORK,PSNOWRHO(:,1),PSNOWAGE(:,1),PPERMSNOWFRAC,PPS)
!
DO JJ=1,SIZE(ZSPECTRALALBEDO,2)
   DO JI=1,INI
      IF(ZWORK(JI)==0.0.AND.ZSNOWMASSNEW(JI)>0.0)THEN
         ZSPECTRALALBEDO(JI,JJ) = ZSPECTRALWORK(JI,JJ)
      ENDIF
   ENDDO
ENDDO
WHERE(ZWORK(:)==0.0.AND.ZSNOWMASSNEW(:)>0.0)
     PSNOWALB(:) = ZWORK2(:)
ENDWHERE
!
!*       3.     Update grid/discretization
!               --------------------------
! Reset grid to conform to model specifications:
!
CALL SNOW3LGRID(ZSNOWDZN,ZSNOW,PSNOWDZ_OLD=PSNOWDZ)
!
! Mass/Heat redistribution:
!
CALL SNOW3LTRANSF(ZSNOW,PSNOWDZ,ZSNOWDZN,PSNOWRHO,PSNOWHEAT,PSNOWAGE)
!
!
!*       4.     Liquid water content and snow temperature
!               -----------------------------------------
!
! First diagnose snow temperatures and liquid
! water portion of the snow from snow heat content:
!
ZSCAP(:,:)     = SNOW3LSCAP(PSNOWRHO)
!
ZSNOWTEMP(:,:) = XTT + ( ((PSNOWHEAT(:,:)/PSNOWDZ(:,:))                   &
                   + XLMTT*PSNOWRHO(:,:))/ZSCAP(:,:) )  
!
PSNOWLIQ(:,:)  = MAX(0.0,ZSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*                  &
                   PSNOWDZ(:,:)/(XLMTT*XRHOLW)  
!
ZSNOWTEMP(:,:) = MIN(XTT,ZSNOWTEMP(:,:))
!
!
!*       5.     Snow Compaction
!               ---------------
!
! Calculate snow density: compaction/aging: density increases
!
CALL SNOW3LCOMPACTN(PTSTEP,XSNOWDZMIN,PSNOWRHO,PSNOWDZ,ZSNOWTEMP,ZSNOW,PSNOWLIQ)
!
! Snow compaction and metamorphism due to drift
!
PSNDRIFT(:) = 0.0
IF (HSNOWDRIFT == 'DFLT') THEN
   CALL SNOW3LDRIFT(PTSTEP,PFORESTFRAC,PWIND_DRIFT,PTAR,PQA,PPS,PRHOA,&
                    PSNOWRHO,PSNOWDZ,ZSNOW,OSNOWDRIFT_SUBLIM,PSNDRIFT)
ENDIF
!
! Update snow heat content (J/m2):
!
ZSCAP(:,:)     = SNOW3LSCAP(PSNOWRHO)
PSNOWHEAT(:,:) = PSNOWDZ(:,:)*( ZSCAP(:,:)*(ZSNOWTEMP(:,:)-XTT)        &
                   - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:)  
!
!
!*       6.     Solar radiation transmission
!               -----------------------------
!
! Heat source (-sink) term due to shortwave
! radiation transmission within the snowpack:
!
CALL SNOW3LRAD(OMEB,XSNOWDZMIN,PSW_RAD,PSNOWALB,      &
               ZSPECTRALALBEDO,PSNOWDZ,PSNOWRHO,PALB, &
               PPERMSNOWFRAC,PZENITH,PSWNETSNOW,      &
               PSWNETSNOWS,ZRADSINK,ZRADXS,PSNOWAGE)  
!
!
!*       7.     Heat transfer and surface energy budget
!               ---------------------------------------
! Snow thermal conductivity:
!
CALL SNOW3LTHRM(PSNOWRHO,ZSCOND,ZSNOWTEMP,PPS)
!
! Precipitation heating term:
! Rainfall renders it's heat to the snow when it enters
! the snowpack:
! NOTE: for now set to zero because we'd need to remove this heat from
!       the atmosphere to conserve energy. This can be done, but should
!       be tested within an atmos model first probably...NOTE, this
!       could be actiavted for use in OFFLINE mode however.
!
!PHPSNOW(:)     = PRR(:)*XCL*(MAX(XTT,PTAR(:))-XTT)    ! (W/m2)
PHPSNOW(:) = 0.0
!
! Surface Energy Budget calculations using ISBA linearized form
! and standard ISBA turbulent transfer formulation
!
IF(OMEB)THEN 

! - surface fluxes prescribed:

   CALL SNOW3LEBUDMEB(PTSTEP,XSNOWDZMIN,                                     &                
                ZSNOWTEMP(:,1),PSNOWDZ(:,1),PSNOWDZ(:,2),                    & 
                ZSCOND(:,1),ZSCOND(:,2),ZSCAP(:,1),                          &
                PSWNETSNOWS,PLWNETSNOW,                                      &
                PHSNOW,PLES3L,PLEL3L,PHPSNOW,                                &   
                ZCT,ZTSTERM1,ZTSTERM2,PGFLUXSNOW)

! - for mass fluxes, snow fraction notion removed

   ZPSN3L(:) = 1.0

ELSE

   ZPSN3L(:) = PPSN3L(:)

! - compute surface energy budget and fluxes:

   CALL SNOW3LEBUD(HSNOWRES, HIMPLICIT_WIND,                                   &
                  PPEW_A_COEF, PPEW_B_COEF,                                    &
                  PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,          &
                  XSNOWDZMIN, PPERMSNOWFRAC,                                   &
                  PZREF,ZSNOWTEMP(:,1),PSNOWRHO(:,1),PSNOWLIQ(:,1),ZSCAP(:,1), &
                  ZSCOND(:,1),ZSCOND(:,2),                                     &
                  PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                           &
                  PLW_RAD,PSW_RAD,PTAC,PQA,PPS,PTSTEP,                         &
                  PSNOWDZ(:,1),PSNOWDZ(:,2),PSNOWALB,PZ0,PZ0EFF,PZ0H,          &
                  ZSFCFRZ,ZRADSINK(:,1),PHPSNOW,                               &
                  ZCT,PEMISNOW,PRHOA,ZTSTERM1,ZTSTERM2,ZRA,PCDSNOW,PCHSNOW,    &
                  ZQSAT, ZDQSAT, ZRSRA, ZUSTAR2_IC, PRI,                       &
                  ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T      )
!
ENDIF
!
! Heat transfer: simple diffusion along the thermal gradient
!
ZSNOWTEMP01(:) = ZSNOWTEMP(:,1) ! save surface snow temperature before update
!
ZGRNDFLUXI(:)  = ZGRNDFLUX(:)
!
CALL SNOW3LSOLVT(OMEB,PTSTEP,XSNOWDZMIN,PSNOWDZ,ZSCOND,ZSCAP,PTG,              &
                   PSOILCOND,PD_G,ZRADSINK,ZCT,ZTSTERM1,ZTSTERM2,              &
                   ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T,    &
                   ZTA_IC,ZQA_IC,ZGRNDFLUX,ZGRNDFLUXO,ZSNOWTEMP,PRESTOREN      )  
!
!
!*       8.     Surface fluxes
!               --------------
! 
! Since surface fluxes already computed under MEB option (and already
! recomputed for the case when T>Tf), Only call if MEB not in use:
!
IF(.NOT.OMEB)THEN 

   CALL SNOW3LFLUX(ZSNOWTEMP(:,1),PEXNS,PEXNA,ZUSTAR2_IC,             &
                  PTSTEP,PSNOWALB,PSW_RAD,                            &
                  PEMISNOW,ZLWUPSNOW,PLW_RAD,PLWNETSNOW,              &
                  ZTA_IC,ZSFCFRZ,ZQA_IC,PHPSNOW,                      &
                  ZSNOWTEMP01,PRESTOREN,ZCT,ZRADSINK(:,1),            &
                  ZQSAT,ZDQSAT,ZRSRA,                                 &
                  PRNSNOW,PHSNOW,PGFLUXSNOW,PLES3L,PLEL3L,PEVAP,      &
                  PUSTAR,GSFCMELT                                     )
!
ENDIF                  
!
!
!*       9.     Snow melt
!               ---------
!
! First Test to see if snowpack vanishes during this time step:
!
CALL SNOW3LGONE(PTSTEP,PLEL3L,PLES3L,PSNOWRHO,                            &
                PSNOWHEAT,ZRADSINK(:,INLVLS),PEVAPCOR,PTHRUFAL,ZGRNDFLUX, &
                PGFLUXSNOW,ZGRNDFLUXO,PSNOWDZ,PSNOWLIQ,ZSNOWTEMP,         &
                PLVTT,PLSTT,ZRADXS,PSNOWMELT,ZSNOWGONE_DELTA              )  
!
! For "normal" melt: transform excess heat content into snow liquid:
!
ZSNOWLIQ0(:,:) = PSNOWLIQ(:,:) ! save liquid water profile before update
!
CALL SNOW3LMELT(PTSTEP,ZSCAP,ZSNOWTEMP,PSNOWDZ,PSNOWRHO,PSNOWLIQ,ZMELTXS)  
!
!
!*      10.     Snow water flow, refreezing, and water budget
!               ---------------------------------------------
!
! Liquid water vertical transfer and possible snowpack runoff
! And refreezing/freezing of meltwater/rainfall (ripening of the snow)
!
CALL SNOW3LREFRZ(PTSTEP,PRR,PSNOWRHO,ZSNOWTEMP,PSNOWDZ,PSNOWLIQ,  &
                  ZSNOWLIQ0,PTHRUFAL,ZDELPHASE_SFC,ZDELPHASE      )
!
ZSCAP(:,:)        = SNOW3LSCAP(PSNOWRHO)
PSNOWHEAT(:,:)    = PSNOWDZ(:,:)*( ZSCAP(:,:)*(ZSNOWTEMP(:,:)-XTT)        &
                      - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:) 
!
! Latent heating due to surface layer snow melt/freeze  (W m-2)
!
PDELPHASEN_SFC(:)=ZSNOWGONE_DELTA(:)*ZDELPHASE_SFC(:)*XLMTT
!
! Snowmelt / water refreezing in the snowpack  [kg/(m2 s)]
!
PSNOWMELT (:)=ZSNOWGONE_DELTA(:)*MAX(0.0, ZDELPHASE(:))
PSNREFREEZ(:)=ZSNOWGONE_DELTA(:)*MAX(0.0,-ZDELPHASE(:))
!
! Latent heating due to snow melt/freeze  (W m-2)
!
PDELPHASEN(:)=ZSNOWGONE_DELTA(:)*ZDELPHASE(:)*XLMTT
!
!*      11.     Snow Evaporation/Sublimation mass updates:
!               ------------------------------------------
!
CALL SNOW3LEVAPN(ZPSN3L,PLES3L,PLEL3L,PTSTEP,PPERMSNOWFRAC,       &
                 ZSNOWTEMP(:,1),PSNOWRHO(:,1),PSNOWDZ,            &
                 PSNOWLIQ(:,1),PTAR,PLVTT,PLSTT,PSNOWHEAT,PSOILCOR)
!
! Update snow temperatures and liquid
! water portion of the snow from snow heat content
! since vertical distribution of enthalpy might have been
! modified in the previous routine:
!
ZSCAP(:,:)     = SNOW3LSCAP(PSNOWRHO)
ZWORK2D(:,:)   = MIN(1.0, PSNOWDZ(:,:)/XSNOWDZMIN)
ZSNOWTEMP(:,:) = XTT + ZWORK2D(:,:)*( ((PSNOWHEAT(:,:)/MAX(XSNOWDZMIN,PSNOWDZ(:,:))) &
                   + XLMTT*PSNOWRHO(:,:))/ZSCAP(:,:) ) 
PSNOWLIQ(:,:)  = MAX(0.0,ZSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*PSNOWDZ(:,:)/(MAX(ZEPSILON,ZWORK2D(:,:))*XLMTT*XRHOLW)
ZSNOWTEMP(:,:) = MIN(XTT,ZSNOWTEMP(:,:))
!
!
! If all snow in uppermost layer evaporates/sublimates, re-distribute
! grid (below could be evoked for vanishingly thin snowpacks):
!
CALL SNOW3LEVAPGONE(PSNOWHEAT,PSNOWDZ,PSNOWRHO,ZSNOWTEMP,PSNOWLIQ)
!
!
!*      12.     Update surface albedo:
!               ----------------------
! Snow clear sky albedo:
!
 CALL SNOW3LALB(PSNOWALB,ZSPECTRALALBEDO,PSNOWRHO(:,1),PSNOWAGE(:,1),PPERMSNOWFRAC,PPS)
!
!
!*      13.     Update snow heat content:
!               -------------------------
! Update the heat content (variable stored each time step)
! using current snow temperature and liquid water content:
!
! First, make check to make sure heat content not too large
! (this can result due to signifigant heating of thin snowpacks):
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZLIQHEATXS(JI) = MAX(0.0,PSNOWLIQ(JI,JJ)*XRHOLW-PSNOWDZ(JI,JJ)*PSNOWRHO(JI,JJ))*XLMTT/PTSTEP  
      PSNOWLIQ(JI,JJ)= PSNOWLIQ(JI,JJ) - ZLIQHEATXS(JI)*PTSTEP/(XRHOLW*XLMTT)
      PSNOWLIQ(JI,JJ)= MAX(0.0, PSNOWLIQ(JI,JJ))
   ENDDO
ENDDO
!
PSNOWTEMP(:,:)    = ZSNOWTEMP(:,:)
!
ZSCAP(:,:)        = SNOW3LSCAP(PSNOWRHO)
!
PSNOWHEAT(:,:)    = PSNOWDZ(:,:)*( ZSCAP(:,:)*(PSNOWTEMP(:,:)-XTT)        &
                  - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:)  
!
!
!*      14.     Snow/Ground heat flux:
!               ----------------------
!
!Add radiation not absorbed by snow to soil/vegetation interface flux
!
PGRNDFLUX(:) = ZGRNDFLUXO(:)+ZRADXS(:)
!
!Comput soil/snow correction heat flux to prevent TG1 numerical oscillations
!Add any excess heat for melting and liquid water :
!
PGFLXCOR(:) = (ZGRNDFLUX(:)-ZGRNDFLUXO(:))+ZMELTXS(:)+ZLIQHEATXS(:)
!
!
!*      15.     Update snow mass and age:
!               -------------------------
!
PSNOWSWE(:,:) = PSNOWDZ(:,:)*PSNOWRHO(:,:)
!
WHERE (PSNOWSWE(:,:)>0)
  PSNOWAGE(:,:)=PSNOWAGE(:,:)+(PTSTEP/XDAY)
ELSEWHERE
  PSNOWAGE(:,:)= XUNDEF
ENDWHERE
!
!
!*      16.     Update surface specific humidity:
!               ---------------------------------
!
IF(OMEB)THEN
   PQS(:) = QSATI(PSNOWTEMP(:,1),PPS) ! purely diagnostic
ELSE
   PQS(:) = ZQSAT(:)
ENDIF
!
!
!*      17.     MEB related checks
!               ------------------
!
IF(OMEB)THEN

! Final adjustment: If using MEB, then all excess heat (above initial guess ground flux)
! is used herein to adjust the ground temperature since for the MEB case, boundary
! fluxes are imposed (the energy budget of the ground has already been computed).
! This correction is needed since snow-soil flux here is implicit, while guess
! used in surface energy soil budget was semi-implicit. This forces the flux
! seen by the ground to be *the same* as that leaving the snow (energy conservation).
! Also, add any excessing cooling from sublimation as snowpack becomes vanishingly thin.
! This is added back to total heat content of the snowpack (and distributed among
! all snow layers while conserving total heat content plus correction)

   ZWORK(:)               = 0.
   DO JJ=1,INLVLS
      DO JI=1,INI
         ZWORK(JI)        = ZWORK(JI) + PSNOWHEAT(JI,JJ)
      ENDDO
   ENDDO
   ZWORK2(:)              = MIN(0.0, ZWORK(:) + PGRNDFLUX(:) - ZGRNDFLUXI(:))
   PGFLXCOR(:)            = PGFLXCOR(:) + MAX(0., ZWORK2(:)) ! add any possible (rare!) excess to soil

   WHERE(ZWORK(:) > -1.E-10)
      ZWORK(:) = 1. ! i.e. no modifs to H profile
   ELSEWHERE
      ZWORK(:) = ZWORK2(:)/ZWORK(:)
   END WHERE

   DO JJ=1,INLVLS
      DO JI=1,INI
         PSNOWHEAT(JI,JJ) = PSNOWHEAT(JI,JJ)*ZWORK(JI)
      ENDDO
   ENDDO

! these are just updated diagnostics at this point:

   ZWORK2D(:,:)   = MIN(1.0, PSNOWDZ(:,:)/XSNOWDZMIN)/MAX(XSNOWDZMIN,PSNOWDZ(:,:))
   PSNOWTEMP(:,:) = XTT + ZWORK2D(:,:)*( (PSNOWHEAT(:,:) + XLMTT*PSNOWSWE(:,:))/ZSCAP(:,:) )
   PSNOWLIQ(:,:)  = MAX(0.0,PSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*PSNOWDZ(:,:)/(XLMTT*XRHOLW)  
   PSNOWTEMP(:,:) = MIN(XTT,PSNOWTEMP(:,:))

ENDIF
!
!
!*      18.     Energy Budget Diagnostics:
!               --------------------------
!
! Snow heat storage change over the time step excluding phase change (W m-2):
!
PDELHEATN(:) = -PDELPHASEN(:)-PSNOWHMASS(:)/PTSTEP
DO JJ=1,INLVLS
   DO JI=1,INI
      PDELHEATN(JI) = PDELHEATN(JI) + (PSNOWHEAT(JI,JJ)-ZSNOWHEAT0(JI,JJ))/PTSTEP
   ENDDO
END DO
!
! where ZSNOWHEAT0 is the value of PSNOWHEAT that enters this routine before any adjustments/processes.
! Errors should be VERY small ( ~0) at each time step.
! The above is "snow relative"...the actual tile or grid box error is multiplied by PPSN3L.
!
!
PDELHEATN_SFC(:) = (PSNOWHEAT(:,1)-ZSNOWHEAT0(:,1)-ZSNOWHMASS1(:))/PTSTEP
!
! The pseudo-heating (W m-2) is derrived as a diagnostic from the surface energy budget
! here...NOTE that the total snow energy budget is explicitly computed, and since
! it balances, then simply compute this term as a residue (which ensures sfc energy balance):
!
PSNOWSFCH(:) = PDELHEATN_SFC(:)+PDELPHASEN_SFC(:)-PSWNETSNOWS(:)-PLWNETSNOW(:) &
                               +PHSNOW(:)+PLES3L(:)+PLEL3L(:)+PRESTOREN(:)
!
! NOTE: To check enthalpy conservation, the error in W m-2 is defined HERE as
!
! total error   = PDELHEATN(:)+PDELPHASEN(:)-PGFLUXSNOW(:)+PGRNDFLUX(:)+PGFLXCOR(:)
!
! surface error = PDELHEATN_SFC(:)-PSNOWSFCH(:)+PSWNETSNOW(:)-PSWNETSNOWS(:)-PGFLUXSNOW(:)+PRESTOREN(:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3L',1,ZHOOK_HANDLE)
!
!
CONTAINS
!
!
!
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LDRIFT(PTSTEP,PFORESTFRAC,PVMOD,PTA,PQA,PPS,PRHOA,&
                       PSNOWRHO,PSNOWDZ,PSNOW,OSNOWDRIFT_SUBLIM,PSNDRIFT)
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_SNOW_PAR, ONLY : XVTIME, XVROMAX, XVROMIN, XVMOB1, &
                          XVDRIFT1, XVDRIFT2, XVDRIFT3,     &
                          XCOEF_FF, XCOEF_EFFECT, XQS_REF
!
USE MODE_THERMOS
!
!
!!    PURPOSE
!!    -------
!     Snow compaction  and metamorphism due to drift
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method inspired from
!     Brun et al. (1997) and Guyomarch
!     
!     - computes a mobility index of each snow layer from its  density
!     - computes a drift index of each layer from its mobility and wind speed
!     - computes a transport index 
!     - increases density in case of transport
!
!     HISTORY:
!     Basic parameterization from Crocus/ARPEGE Coupling (1997) 
!     Adaptation from Crocus SNOWDRIFT subroutine
!
!!    PURPOSE
!!    -------
!     Snow compaction due to overburden and settling.
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. 
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PFORESTFRAC
REAL, DIMENSION(:), INTENT(IN)      :: PVMOD
REAL, DIMENSION(:), INTENT(IN)      :: PTA
REAL, DIMENSION(:), INTENT(IN)      :: PQA
REAL, DIMENSION(:), INTENT(IN)      :: PPS
REAL, DIMENSION(:), INTENT(IN)      :: PRHOA
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW
!
LOGICAL,            INTENT(IN)      :: OSNOWDRIFT_SUBLIM
REAL, DIMENSION(:), INTENT(OUT)     :: PSNDRIFT !blowing snow sublimation (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRMOB
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRDRIFT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDRO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZQS_EFFECT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDRIFT_EFFECT
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZPROFEQU
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZWIND
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZQSATI
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZVT
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZQS
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZW
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZT
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZSNOWDZ1
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZFOREST_EFFECT
!
LOGICAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: GDRIFT
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LDRIFT',0,ZHOOK_HANDLE)
!
! 0. Initialization:
! ------------------
!
INI    = SIZE(PSNOWDZ(:,:),1)
INLVLS = SIZE(PSNOWDZ(:,:),2)
!
ZSNOWRHO2    (:,:) = PSNOWRHO(:,:)
ZRDRIFT      (:,:) = XUNDEF
ZRT          (:,:) = XUNDEF
ZDRO         (:,:) = XUNDEF
ZQS_EFFECT   (:,:) = 0.0
ZDRIFT_EFFECT(:,:) = 0.0
GDRIFT       (:,:) = .FALSE.
!
ZPROFEQU (:) = 0.0
!
ZSNOWDZ1(:) = PSNOWDZ(:,1)
!
! 1. Initialazation of drift and induced settling
! -----------------------------------------------
!
!Limitation of wind speed in Forest : ~ 15% of wind in open area
ZFOREST_EFFECT(:) = 1.0 - 0.85 * PFORESTFRAC(:)
!
ZWIND(:) = XCOEF_FF * ZFOREST_EFFECT(:) * PVMOD(:)
!
DO JJ=1,INLVLS
   DO JI=1,INI
!      
!     mobility index computation of a layer as a function of its density 
      ZRMOB(JI,JJ)= XCOEF_FF-XCOEF_FF*1.0E-3*(MAX(PSNOWRHO(JI,JJ),XVROMIN)-XVROMIN) / XVMOB1
!
!     computation of the drift index inclunding the decay by overburden snow 
      ZRDRIFT(JI,JJ) =  ZRMOB(JI,JJ)-(XVDRIFT1*EXP(-XVDRIFT2*ZWIND(JI))-1.0)
!
   ENDDO
ENDDO
!
!Compaction from the surface to layers beneath until 
!a snow layer having negative drift index
GDRIFT(:,:) = (ZRDRIFT(:,:)>0.0)
DO JI=1,INI
   DO JJ=1,INLVLS
      IF(.NOT.GDRIFT(JI,JJ))THEN
        GDRIFT(JI,JJ:INLVLS)=.FALSE.
        EXIT
      ENDIF
   ENDDO
ENDDO

! 2. Specific case for blowing snow sublimation
! ---------------------------------------------
!
IF(OSNOWDRIFT_SUBLIM)THEN
!
  ZQSATI(:)= QSATI(PTA(:),PPS(:))
!
! computation of wind speed threshold QSATI and RH withe respect to ice
!
  ZW (:)=MAX(-0.99,ZRMOB(:,1))
  ZVT(:)=LOG(XVDRIFT1/(1.0+ZW(:)))/XVDRIFT2
!
  WHERE(ZWIND(:)>0.0)
    ZW(:)=LOG(ZWIND(:)/ZVT(:))
    ZW(:)=EXP(3.6*ZW(:))
  ELSEWHERE
    ZW(:)=0.0
  ENDWHERE
!
! computation of sublimation rate according to Gordon's PhD
!
  ZT(:)=XTT/PTA(:)
  ZT(:)=0.0018*(ZT(:)**4)
!
  ZQS(:)=ZT(:)*ZVT(:)*PRHOA(:)*ZQSATI(:)*(1.-PQA(:)/ZQSATI(:))*ZW(:)
!
! surface depth decrease in case of blowing snow sublimation
!
  PSNOWDZ(:,1)=MAX(0.5*PSNOWDZ(:,1),PSNOWDZ(:,1)-MAX(0.,ZQS(:))*PTSTEP/(XCOEF_FF*PSNOWRHO(:,1)))
!
  PSNDRIFT(:) = (ZSNOWDZ1(:)-PSNOWDZ(:,1))*PSNOWRHO(:,1)/PTSTEP
!
  ZQS_EFFECT(:,1)=MIN(3.,MAX(0.,ZQS(:))/XQS_REF)
!
ENDIF
!
! 3. Computation of drift and induced settling
! --------------------------------------------
!
DO JJ=1,INLVLS
   DO JI=1,INI

!     update the decay coeff by half the current layer
      ZPROFEQU(JI) = ZPROFEQU(JI)  + 0.5 * PSNOWDZ(JI,JJ) * 0.1 * (XVDRIFT3-ZRDRIFT(JI,JJ))
!
      IF(GDRIFT(JI,JJ).AND.PSNOWRHO(JI,JJ)<XVROMAX)THEN
!      
!       computation of the drift index inclunding the decay by overburden snow 
        ZRT(JI,JJ) = MAX(0.0,ZRDRIFT(JI,JJ)*EXP(-ZPROFEQU(JI)*100.0))
!     
        ZDRIFT_EFFECT(JI,JJ) = (ZQS_EFFECT(JI,JJ)+XCOEF_EFFECT)*ZRT(JI,JJ)/(XVTIME*XCOEF_FF)
!
!       settling by wind transport only in case of not too dense snow
        ZDRO(JI,JJ) = (XVROMAX - PSNOWRHO(JI,JJ)) * ZDRIFT_EFFECT(JI,JJ) * PTSTEP
!          
!       Calculate new snow density:
        ZSNOWRHO2(JI,JJ) = MIN(XVROMAX,PSNOWRHO(JI,JJ)+ZDRO(JI,JJ))

!       Conserve mass by decreasing grid thicknesses in response to density increases
        PSNOWDZ(JI,JJ) = PSNOWDZ(JI,JJ)*(PSNOWRHO(JI,JJ)/ZSNOWRHO2(JI,JJ))    
!
      ENDIF
!       
!     update the decay coeff by half the current layer
      ZPROFEQU(JI) = ZPROFEQU(JI)  + 0.5 * PSNOWDZ(JI,JJ) * 0.1 * (XVDRIFT3-ZRDRIFT(JI,JJ))
!      
   ENDDO
ENDDO
!
! 4. Update total snow depth and density profile:
! -----------------------------------------------
!
! Compaction/augmentation of total snowpack depth
!
PSNOW(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOW(JI) = PSNOW(JI) + PSNOWDZ(JI,JJ)
   ENDDO
ENDDO
!
! Update density (kg m-3):
!
PSNOWRHO(:,:)  = ZSNOWRHO2(:,:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LDRIFT',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LDRIFT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LRAD(OMEB, PSNOWDZMIN, PSW_RAD, PSNOWALB,      &
                           PSPECTRALALBEDO, PSNOWDZ, PSNOWRHO, PALB, &
                           PPERMSNOWFRAC, PZENITH,  PSWNETSNOW,      &
                           PSWNETSNOWS, PRADSINK, PRADXS, PSNOWAGE   )  
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave (solar) radiation
!     through the snowpack (using a form of Beer's Law: exponential
!     decay of radiation with increasing snow depth).
!
USE MODD_SNOW_PAR, ONLY : XMINCOSZEN
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
!
USE MODE_SNOW3L,   ONLY : SNOW3LDOPT, SNOW3LRADABS_SFC
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,            INTENT(IN)      :: OMEB ! if=T, then uppermost abs is diagnosed
!                                           !       since fluxes known
!
REAL,               INTENT(IN)      :: PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)      :: PSW_RAD
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWALB
REAL, DIMENSION(:), INTENT(IN)      :: PALB
REAL, DIMENSION(:), INTENT(IN)      :: PPERMSNOWFRAC
REAL, DIMENSION(:), INTENT(IN)      :: PZENITH
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ, PSNOWAGE
REAL, DIMENSION(:,:), INTENT(IN)    :: PSPECTRALALBEDO
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSWNETSNOW, PSWNETSNOWS
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRADXS
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRADSINK
!
!
!*      0.2    declarations of local variables
!
INTEGER                              :: JJ, JI
!
INTEGER                              :: INI
INTEGER                              :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZRADTOT, ZPROJLAT, ZCOSZEN
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDSGRAIN, ZCOEF, ZSNOWDZ, ZAGE
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZBETA1, ZBETA2, ZBETA3, ZWORK
REAL, DIMENSION(SIZE(PSPECTRALALBEDO,1),SIZE(PSPECTRALALBEDO,2)) :: ZSPECTRALALBEDO
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LRAD',0,ZHOOK_HANDLE)
!
INI    = SIZE(PSNOWDZ(:,:),1)
INLVLS = SIZE(PSNOWDZ(:,:),2)
!
ZSPECTRALALBEDO(:,:) = XUNDEF ! Init
!
! 1. Vanishingly thin snowpack check:
! -----------------------------------
!    For vanishingly thin snowpacks, much of the radiation
!    can pass through snowpack into underlying soil, making
!    a large (albeit temporary) thermal gradient: by imposing
!    a minimum thickness, this increases the radiation absorbtion
!    for vanishingly thin snowpacks.
!
ZSNOWDZ(:,:) = MAX(PSNOWDZMIN, PSNOWDZ(:,:))
!
!
! 2. Extinction of net shortwave radiation
! ----------------------------------------
! Fn of snow depth and density (Loth and Graf 1993:
! SNOWCVEXT => from Bohren and Barkstrom 1974
! SNOWAGRAIN and SNOWBGRAIN=> from Jordan 1976)
!
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
!
ZCOSZEN(:)=MAX(XMINCOSZEN,COS(PZENITH(:)))
!
! This formulation is incorrect but it compensate partly the fact that 
! the albedo formulation does not account for zenithal angle.
! Only for polar or glacier regions
!
ZPROJLAT(:)=(1.0-PPERMSNOWFRAC(:))+PPERMSNOWFRAC(:)/ZCOSZEN(:)
!
! Snow optical grain diameter (no age dependency over polar regions):
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZAGE(JI,JJ) = (1.0-PPERMSNOWFRAC(JI))*PSNOWAGE(JI,JJ)
   ENDDO
ENDDO
!
ZDSGRAIN(:,:) = SNOW3LDOPT(PSNOWRHO,ZAGE)
!
! Calculate the transmission of shortwave radiation within the snowpack:
!
IF(OMEB)THEN
!
! Only 2 bands currently considered (for vegetation and soil)
! thus eliminate a band and renormalize (as was done for the surface snow albedo
! for the MEB snow surface energy budget computations)
!
   ZSPECTRALALBEDO(:,1) = PSPECTRALALBEDO(:,1)
   ZSPECTRALALBEDO(:,2) = (PSNOWALB(:) - XSW_WGHT_VIS*ZSPECTRALALBEDO(:,1))/XSW_WGHT_NIR
!
   ZCOEF(:,:)           = SNOW3LRADABS_SFC(PSNOWRHO,ZSNOWDZ,ZSPECTRALALBEDO,PZENITH,PPERMSNOWFRAC,ZDSGRAIN)
!
! Diagnose surface layer coef (should be very close/identical to ZCOEF(:,1) computed above)
!
   ZCOEF(:,1)           = (PSWNETSNOW(:)-PSWNETSNOWS(:))/MAX(1.E-4,PSW_RAD(:))
!
ELSE
!
! Consider 3 bands:
!   
   ZCOEF(:,:)           = SNOW3LRADABS_SFC(PSNOWRHO,ZSNOWDZ,PSPECTRALALBEDO,PZENITH,PPERMSNOWFRAC,ZDSGRAIN)
!
   PSWNETSNOW(:)       = PSW_RAD(:)*(1.-PSNOWALB(:))
   PSWNETSNOWS(:)      = PSWNETSNOW(:) - PSW_RAD(:)*ZCOEF(:,1)
!
ENDIF
!
! 3. Radiation at each level: (W/m2)
! ----------------------------------
!
! NOTE: as this is a "sink" term for heat, it
! as assigned a negative value.
DO JJ=1,INLVLS
   DO JI=1,INI
      PRADSINK(JI,JJ)  = -PSW_RAD(JI)*ZCOEF(JI,JJ)
   ENDDO
ENDDO
!
! For thin snowpacks, radiation might reach base of
! snowpack...so we influence this amount with sfc albedo
! and (outside of this routine) add any excess heat
! to underlying soil:
!
PRADSINK(:,INLVLS) = PRADSINK(:,INLVLS)*(1.0-PALB(:))
!
! 4. Excess radiation not absorbed by snowpack (W/m2):
! ----------------------------------------------------
!
ZRADTOT(:)    = PRADSINK(:,1) + (1.-PSNOWALB(:))*PSW_RAD(:)
DO JJ=2,INLVLS
   DO JI=1,INI
      ZRADTOT(JI) = ZRADTOT(JI) + PRADSINK(JI,JJ)-PRADSINK(JI,JJ-1)
   ENDDO
ENDDO
!
PRADXS(:)     = (1.-PSNOWALB(:))*PSW_RAD(:) - ZRADTOT(:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LRAD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LRAD
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LEBUD(HSNOWRES, HIMPLICIT_WIND,                                     &
                              PPEW_A_COEF, PPEW_B_COEF,                                   &
                              PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,         &
                              PSNOWDZMIN, PPERMSNOWFRAC,                                  &
                              PZREF,PTS,PSNOWRHO,PSNOWLIQ,PSCAP,PSCOND1,PSCOND2,          &
                              PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                          &
                              PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                         &
                              PSNOWDZ1,PSNOWDZ2,PALBT,PZ0,PZ0EFF,PZ0H,                    &
                              PSFCFRZ,PRADSINK,PHPSNOW,                                   &
                              PCT,PEMIST,PRHOA,PTSTERM1,PTSTERM2,PRA,PCDSNOW,PCHSNOW,     &
                              PQSAT,PDQSAT,PRSRA,PUSTAR2_IC,PRI,                          &
                              PPET_A_COEF_T,PPEQ_A_COEF_T,PPET_B_COEF_T,PPEQ_B_COEF_T     )  
!
!!    PURPOSE
!!    -------
!     Calculate surface energy budget linearization (terms) and turbulent
!     exchange coefficients/resistance between surface and atmosphere.
!     (Noilhan and Planton 1989; Giordani 1993; Noilhan and Mahfouf 1996)
!
USE MODD_CSTS,     ONLY : XCPD, XRHOLW, XSTEFAN, XLVTT, XLSTT
USE MODD_SNOW_PAR, ONLY : X_RI_MAX, XEMISSN
!
USE MODE_THERMOS
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
USE MODI_SURFACE_CD
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP, PSNOWDZMIN
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES ! type of sfc resistance
!                                      DEFAULT=Louis (1979), standard ISBA
!                                      method. Option to limit Ri number
!                                      for very stable conditions (RIL, RI1, RI2)
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)      :: PPERMSNOWFRAC
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                         PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF, PTS, PSNOWDZ1, PSNOWDZ2,          &
                                         PRADSINK, PSNOWRHO, PSNOWLIQ, PSCAP,   &
                                         PSCOND1, PSCOND2,                      &
                                         PZ0, PHPSNOW,                          &
                                         PALBT, PZ0EFF, PZ0H  
!
REAL, DIMENSION(:), INTENT(IN)    :: PSW_RAD, PLW_RAD, PTA, PQA, PPS, PRHOA
!
REAL, DIMENSION(:), INTENT(IN)    :: PUREF, PEXNS, PEXNA, PDIRCOSZW, PVMOD
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTSTERM1, PTSTERM2, PEMIST, PRA,         &
                                       PCT, PSFCFRZ, PCDSNOW, PCHSNOW,          &
                                       PQSAT, PDQSAT, PRSRA  
!
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR2_IC,                        &
                                       PPET_A_COEF_T, PPEQ_A_COEF_T,    &
                                       PPET_B_COEF_T, PPEQ_B_COEF_T 
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRI
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTS))        :: ZAC, ZRI, ZCOND1, ZCOND2,          &
                                       ZSCONDA, ZA, ZB, ZC,             &
                                       ZCDN, ZSNOWDZM1, ZSNOWDZM2,      &
                                       ZVMOD, ZUSTAR2, ZTS3, ZLVT,      &
                                       Z_CCOEF  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 1. New saturated specific humidity and derrivative:
! ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEBUD',0,ZHOOK_HANDLE)
!
ZRI   (:) = XUNDEF
!
PQSAT (:) =  QSATI(PTS(:),PPS(:))
PDQSAT(:) = DQSATI(PTS(:),PPS(:),PQSAT(:))
!
!
! 2. Surface properties:
! ----------------------
! It might be of interest to use snow-specific roughness
! or a temperature dependence on emissivity:
! but for now, use ISBA defaults.
!
PEMIST(:) = XEMISSN
!
! 2. Computation of resistance and drag coefficient
! -------------------------------------------------
!
 CALL SURFACE_RI(PTS, PQSAT, PEXNS, PEXNA, PTA, PQA,                  &
                PZREF, PUREF, PDIRCOSZW, PVMOD, ZRI                  )  
!
! Simple adaptation of method by Martin and Lejeune (1998)
! to apply a lower limit to turbulent transfer coef
! by defining a maximum Richarson number for stable
! conditions:
!
IF(HSNOWRES=='RIL') THEN
  DO JJ = 1,SIZE(ZRI)
    ZRI(JJ) = MIN(X_RI_MAX,ZRI(JJ))
  ENDDO
ENDIF
!
PRI(:)=ZRI(:)
!
! Surface aerodynamic resistance for heat transfers
!
CALL SURFACE_AERO_COND(ZRI, PZREF, PUREF, PVMOD, PZ0, PZ0H, PRESA_SV, ZAC, PRA, PCHSNOW, HSNOWRES=HSNOWRES)
!
! For atmospheric model coupling:
!
 CALL SURFACE_CD(ZRI, PZREF, PUREF, PZ0EFF, PZ0H, PCDSNOW, ZCDN)
!
PRSRA(:) = PRHOA(:) / PRA(:)
!
!
! Modify flux-form implicit coupling coefficients:
! - wind components:
!
IF(HIMPLICIT_WIND=='OLD')THEN
! old implicitation
  ZUSTAR2(:) = ( PCDSNOW(:)*PVMOD(:)*PPEW_B_COEF(:)) /        &
               (1.0-PRHOA(:)*PCDSNOW(:)*PVMOD(:)*PPEW_A_COEF(:))  
ELSE
! new implicitation
  ZUSTAR2(:) = (PCDSNOW(:)*PVMOD(:)*(2.*PPEW_B_COEF(:)-PVMOD(:)))  &
               / (1.0-2.0*PRHOA(:)*PCDSNOW(:)*PVMOD(:)*PPEW_A_COEF(:))
ENDIF               
!
ZVMOD(:)       = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
ZVMOD(:)       = MAX(ZVMOD(:),0.)
!
WHERE(PPEW_A_COEF(:)/= 0.)
      ZUSTAR2(:) = MAX( ( ZVMOD(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0.)
ENDWHERE
!
! implicit wind friction
ZUSTAR2(:) = MAX(ZUSTAR2(:),0.)
!
PUSTAR2_IC(:) =  ZUSTAR2(:)
!
! 3. Calculate linearized surface energy budget components:
! ---------------------------------------------------------
! To prevent numerical difficulties for very thin snow
! layers, limit the grid "thinness": this is important as
! layers become vanishing thin:
!
ZSNOWDZM1(:) = MAX(PSNOWDZ1(:), PSNOWDZMIN)
ZSNOWDZM2(:) = MAX(PSNOWDZ2(:), PSNOWDZMIN)
!
! Surface thermal inertia:
!
PCT(:)      = 1.0/(PSCAP(:)*ZSNOWDZM1(:))
!
! Fraction of surface frozen (sublimation) with the remaining
! fraction being liquid (evaporation):
! NOTE: currently, for simplicity, assume no liquid water flux
! OFF: PSFCFRZ(:)  = 1.0 - PSNOWLIQ(:)*XRHOLW/(ZSNOWDZM1(:)*PSNOWRHO(:))
!
PSFCFRZ(:)  = 1.0 
!
! Thermal conductivity between uppermost and lower snow layers:
!
ZCOND1 (:) = ZSNOWDZM1(:)/((ZSNOWDZM1(:)+ZSNOWDZM2(:))*PSCOND1(:))
ZCOND2 (:) = ZSNOWDZM2(:)/((ZSNOWDZM1(:)+ZSNOWDZM2(:))*PSCOND2(:))
!
ZSCONDA(:) = 1.0/(ZCOND1(:)+ZCOND2(:))
!
! Transform implicit coupling coefficients: 
! Note, surface humidity is 100% over snow.
!
! - specific humidity:
!
Z_CCOEF(:)       = 1.0 - PPEQ_A_COEF(:)*PRSRA(:)
!
PPEQ_A_COEF_T(:) = - PPEQ_A_COEF(:)*PRSRA(:)*PDQSAT(:)/Z_CCOEF(:)
!
PPEQ_B_COEF_T(:) = ( PPEQ_B_COEF(:) - PPEQ_A_COEF(:)*PRSRA(:)*(PQSAT(:) - &
                       PDQSAT(:)*PTS(:)) )/Z_CCOEF(:)  
!
! - air temperature:
!   (assumes A and B correspond to potential T):
!
Z_CCOEF(:)       = (1.0 - PPET_A_COEF(:)*PRSRA(:))/PEXNA(:)
!
PPET_A_COEF_T(:) = - PPET_A_COEF(:)*PRSRA(:)/(PEXNS(:)*Z_CCOEF(:))
!
PPET_B_COEF_T(:) = PPET_B_COEF(:)/Z_CCOEF(:)
!
!
! Energy budget solution terms:
!
ZTS3(:) = PEMIST(:) * XSTEFAN * PTS(:)**3
ZLVT(:) = (1.-PSFCFRZ(:))*XLVTT + PSFCFRZ(:)*XLSTT
!
ZA(:)   = 1. / PTSTEP + PCT(:) * (4. * ZTS3(:) +                               &
            PRSRA(:) *  ZLVT(:) * (PDQSAT(:) - PPEQ_A_COEF_T(:))                 &
            + PRSRA(:) * XCPD * ( (1./PEXNS(:))-(PPET_A_COEF_T(:)/PEXNA(:)) )    &
            + (2.*ZSCONDA(:)/(ZSNOWDZM2(:)+ZSNOWDZM1(:))) )  
!
ZB(:)   = 1. / PTSTEP + PCT(:) * (3. * ZTS3(:) +                               &
            PRSRA(:) * PDQSAT(:) *  ZLVT(:) )  
!
ZC(:)   = PCT(:) * (PRSRA(:) * XCPD * PPET_B_COEF_T(:)/PEXNA(:) + PSW_RAD(:) * &
            (1. - PALBT(:)) + PEMIST(:)*PLW_RAD(:) - PRSRA(:) *                  &
            ZLVT(:) *  (PQSAT(:)-PPEQ_B_COEF_T(:))                               &
            + PHPSNOW(:) + PRADSINK(:) )  
!
!
! Coefficients needed for implicit solution
! of linearized surface energy budget:
!
PTSTERM2(:) = 2.*ZSCONDA(:)*PCT(:)/(ZA(:)*(ZSNOWDZM2(:)+ZSNOWDZM1(:)))
!
PTSTERM1(:) = (PTS(:)*ZB(:) + ZC(:))/ZA(:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEBUD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LEBUD
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LSOLVT(OMEB,PTSTEP,PSNOWDZMIN,                  &
                               PSNOWDZ,PSCOND,PSCAP,PTG,              &
                               PSOILCOND,PD_G,                        &
                               PRADSINK,PCT,PTERM1,PTERM2,            &
                               PPET_A_COEF_T,PPEQ_A_COEF_T,           &
                               PPET_B_COEF_T,PPEQ_B_COEF_T,           &
                               PTA_IC, PQA_IC,                        &
                               PGRNDFLUX,PGRNDFLUXO,PSNOWTEMP,        &
                               PRESTOREN                              )  
!
!!    PURPOSE
!!    -------
!     This subroutine solves the 1-d diffusion of 'ZSNOWTEMP' using a
!     layer averaged set of equations which are time differenced
!     using the backward-difference scheme (implicit).
!     The eqs are solved rapidly by taking advantage of the
!     fact that the matrix is tridiagonal. This is a very
!     general routine and can be used for the 1-d diffusion of any
!     quantity as long as the diffusity is not a function of the
!     quantity being diffused. Aaron Boone 8-98. Soln to the eq:
!
!                 c  dQ    d  k dQ    dS
!                    -- =  --   -- -  --
!                    dt    dx   dx    dx
!
!     where k = k(x) (thermal conductivity), c = c(x) (heat capacity)
!     as an eg. for temperature/heat eq. S is a sink (-source) term.
!     Diffusivity is k/c
!
!
USE MODD_CSTS,ONLY : XTT
!
USE MODI_TRIDIAG_GROUND
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,            INTENT(IN)      :: OMEB
!
REAL,               INTENT(IN)      :: PTSTEP, PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)      :: PTG, PSOILCOND, PD_G,         &
                                         PCT, PTERM1, PTERM2  

!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ, PSCOND, PSCAP,       &
                                         PRADSINK
!
REAL, DIMENSION(:), INTENT(IN)      :: PPET_A_COEF_T, PPEQ_A_COEF_T, &
                                         PPET_B_COEF_T, PPEQ_B_COEF_T  
!                                       
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGRNDFLUX, PGRNDFLUXO, PRESTOREN,     &
                                         PTA_IC, PQA_IC   
!
!
!*      0.2    declarations of local variables
!
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PTG))                     :: ZSNOWTEMP_DELTA
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWTEMP, ZDTERM, ZCTERM, &
                                                      ZFRCV, ZAMTRX, ZBMTRX,     &
                                                      ZCMTRX  
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZWORK1, ZWORK2, ZDZDIF,    &
                                                      ZSNOWDZM  
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)-1) :: ZSNOWTEMP_M,             &
                                                      ZFRCV_M, ZAMTRX_M,         &
                                                      ZBMTRX_M, ZCMTRX_M  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LSOLVT',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! ------------------
!
ZSNOWTEMP(:,:)  = PSNOWTEMP(:,:)
INI             = SIZE(PSNOWDZ(:,:),1)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Calculate tri-diagnoal matrix coefficients:
! ----------------------------------------------
! For heat transfer, assume a minimum grid
! thickness (to prevent numerical
! problems for very thin snow cover):
!
ZSNOWDZM(:,:)  = MAX(PSNOWDZ(:,:), PSNOWDZMIN)
!
DO JJ=1,INLVLS-1
   DO JI=1,INI
      ZDZDIF(JI,JJ)  = 0.5*(ZSNOWDZM(JI,JJ)+ZSNOWDZM(JI,JJ+1))
      ZWORK1(JI,JJ)  = ZSNOWDZM(JI,JJ  )/(2.0*ZDZDIF(JI,JJ)*PSCOND(JI,JJ  ))
      ZWORK2(JI,JJ)  = ZSNOWDZM(JI,JJ+1)/(2.0*ZDZDIF(JI,JJ)*PSCOND(JI,JJ+1))
   ENDDO
ENDDO
!
ZDZDIF(:,INLVLS) = 0.5*(ZSNOWDZM(:,INLVLS)+PD_G(:))
ZWORK1(:,INLVLS) = ZSNOWDZM(:,INLVLS)/(2.0*ZDZDIF(:,INLVLS)*PSCOND   (:,INLVLS))
ZWORK2(:,INLVLS) = PD_G    (:       )/(2.0*ZDZDIF(:,INLVLS)*PSOILCOND(:       ))
!
ZDTERM(:,:)      = 1.0/(ZDZDIF(:,:)*(ZWORK1(:,:)+ZWORK2(:,:)))
!
ZCTERM(:,:)      = PSCAP(:,:)*ZSNOWDZM(:,:)/PTSTEP
!
! 2. Set up tri-diagonal matrix
! -----------------------------
!
! Upper BC
!
ZAMTRX(:,1) =  0.0
ZBMTRX(:,1) =  1./(PCT(:)*PTSTEP)
ZCMTRX(:,1) = -PTERM2(:)*ZBMTRX(:,1)
ZFRCV(:,1)  =  PTERM1(:)*ZBMTRX(:,1)
!
!
! Interior Grid
!
DO JJ=2,INLVLS-1
   DO JI=1,INI
      ZAMTRX(JI,JJ) = -ZDTERM(JI,JJ-1)
      ZBMTRX(JI,JJ) =  ZCTERM(JI,JJ) + ZDTERM(JI,JJ-1) + ZDTERM(JI,JJ)
      ZCMTRX(JI,JJ) = -ZDTERM(JI,JJ)
      ZFRCV (JI,JJ) =  ZCTERM(JI,JJ)*PSNOWTEMP(JI,JJ) - (PRADSINK(JI,JJ-1)-PRADSINK(JI,JJ))  
   ENDDO
ENDDO
!
!Lower BC
!
ZAMTRX(:,INLVLS) = -ZDTERM(:,INLVLS-1)
ZBMTRX(:,INLVLS) =  ZCTERM(:,INLVLS) + ZDTERM(:,INLVLS-1) + ZDTERM(:,INLVLS)  
ZCMTRX(:,INLVLS) =  0.0
ZFRCV (:,INLVLS) =  ZCTERM(:,INLVLS)*PSNOWTEMP(:,INLVLS) + ZDTERM(:,INLVLS)*PTG(:) &
                 - (PRADSINK(:,INLVLS-1)-PRADSINK(:,INLVLS))  
!
! - - -------------------------------------------------
!
! 4. Compute solution vector
! --------------------------
!
CALL TRIDIAG_GROUND(ZAMTRX,ZBMTRX,ZCMTRX,ZFRCV,ZSNOWTEMP)
!

! Heat flux between surface and 2nd snow layers: (W/m2)
!
PRESTOREN(:) = ZDTERM(:,1)*(ZSNOWTEMP(:,1) - ZSNOWTEMP(:,2))
!
!
! 5. Snow melt case
! -----------------
! If melting in uppermost layer, assume surface layer
! temperature at freezing point and re-evaluate lower
! snowpack temperatures. This is done as it is most likely
! most signigant melting will occur within a time step in surface layer.
! Surface energy budget (and fluxes) will
! be re-calculated (outside of this routine).
!
! NOTE: if MEB is active, then surface fluxes have been defined outside
! of the snow routine and have been adjusted such that they are evaluated
! at a snow surface temperature no greater than Tf. Thus, the implicit surface temperature
! will likely never greatly exceed Tf (before melt computed and they are adjusted to Tf)
! so we can skip the next block of code when MEB is active.
!
IF(.NOT.OMEB)THEN
!
   ZAMTRX_M(:,1)  =  0.0
   ZBMTRX_M(:,1)  =  ZCTERM(:,2) + ZDTERM(:,1) + ZDTERM(:,2)
   ZCMTRX_M(:,1)  = -ZDTERM(:,2)
   ZFRCV_M(:,1)   =  ZCTERM(:,2)*PSNOWTEMP(:,2) + XTT*ZDTERM(:,1)  -  &
                    (PRADSINK(:,1)-PRADSINK(:,2))  
!
   DO JJ=2,INLVLS-1
      DO JI=1,INI
         ZAMTRX_M   (JI,JJ) = ZAMTRX   (JI,JJ+1)
         ZBMTRX_M   (JI,JJ) = ZBMTRX   (JI,JJ+1)
         ZCMTRX_M   (JI,JJ) = ZCMTRX   (JI,JJ+1)
         ZFRCV_M    (JI,JJ) = ZFRCV    (JI,JJ+1)
         ZSNOWTEMP_M(JI,JJ) = PSNOWTEMP(JI,JJ+1)
      ENDDO
   ENDDO
!
   CALL TRIDIAG_GROUND(ZAMTRX_M,ZBMTRX_M,ZCMTRX_M,ZFRCV_M,ZSNOWTEMP_M)
!
! If melting for 2 consecuative time steps, then replace current T-profile
! with one assuming T=Tf in surface layer:
!
   ZSNOWTEMP_DELTA(:)    = 0.0
!
   WHERE(ZSNOWTEMP(:,1) > XTT .AND. PSNOWTEMP(:,1) == XTT)
      PRESTOREN(:)       = ZDTERM(:,1)*(XTT - ZSNOWTEMP_M(:,1))
      ZSNOWTEMP_DELTA(:) = 1.0
   END WHERE
!
   DO JJ=2,INLVLS
      DO JI=1,INI
         ZSNOWTEMP(JI,JJ) = ZSNOWTEMP_DELTA(JI)*ZSNOWTEMP_M(JI,JJ-1)   &
              + (1.0-ZSNOWTEMP_DELTA(JI))*ZSNOWTEMP(JI,JJ)  
      ENDDO
   ENDDO
!
ENDIF
!
!
! 6. Lower boundary flux:
! -----------------------
! NOTE: evaluate this term assuming the snow layer
! can't exceed the freezing point as this adjustment
! is made in melting routine. Then must adjust temperature
! to conserve energy.
! NOTE: if MEB used, surface fluxes are imposed and a test semi-implicit
! ground-snow flux has already been estimated and snow-free ground fluxes
! are generally weaker (and snow fractions are generally set to go to unity
! faster than for the composite soil-veg case), thus this correction
! is not as essential and is off.
!
PGRNDFLUXO(:)          = ZDTERM(:,INLVLS)*(ZSNOWTEMP(:,INLVLS) - PTG(:))
IF(OMEB)THEN
   PGRNDFLUX(:)        = PGRNDFLUXO(:) 
ELSE
   PGRNDFLUX(:)        = ZDTERM(:,INLVLS)*(MIN(XTT,ZSNOWTEMP(:,INLVLS)) - PTG(:))
ENDIF
!
ZSNOWTEMP(:,INLVLS) = ZSNOWTEMP(:,INLVLS) + (PGRNDFLUXO(:)-PGRNDFLUX(:))/ZCTERM(:,INLVLS)
!
! 7. Update temperatute profile in time:
! --------------------------------------
!
PSNOWTEMP(:,:)      = ZSNOWTEMP(:,:)
!
!
! 8. Compute new (implicit) air T and specific humidity
! -----------------------------------------------------
!
IF(.NOT.OMEB)THEN
!
   PTA_IC(:) = PPET_B_COEF_T(:) + PPET_A_COEF_T(:)* PSNOWTEMP(:,1)
!
   PQA_IC(:) = PPEQ_B_COEF_T(:) + PPEQ_A_COEF_T(:)* PSNOWTEMP(:,1)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOW3LSOLVT',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOW3LSOLVT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LMELT(PTSTEP,PSCAP,PSNOWTEMP,PSNOWDZ, &
                            PSNOWRHO,PSNOWLIQ,PMELTXS       )  
!
!
!!    PURPOSE
!!    -------
!     Calculate snow melt (resulting from surface fluxes, ground fluxes,
!     or internal shortwave radiation absorbtion). It is used to
!     augment liquid water content, maintain temperatures
!     at or below freezing, and possibly reduce the mass
!     or compact the layer(s).
!
!
USE MODD_CSTS,ONLY : XTT, XLMTT, XRHOLW, XRHOLI
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCAP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWTEMP,   &
                                       PSNOWRHO, PSNOWLIQ  
!
REAL, DIMENSION(:), INTENT(OUT)     :: PMELTXS
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE, ZCMPRSFACT,   &
                                                      ZSNOWLWE, ZWHOLDMAX,  &
                                                      ZSNOWMELT, ZSNOWTEMP, &
                                                      ZMELTXS 
!
INTEGER :: INI, INL, JI, JL ! loop counter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ---------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LMELT',0,ZHOOK_HANDLE)
!
INI=SIZE(PSNOWRHO,1)
INL=SIZE(PSNOWRHO,2)
!
ZPHASE(:,:)     = 0.0
ZCMPRSFACT(:,:) = 0.0
ZSNOWLWE(:,:)   = 0.0
ZWHOLDMAX(:,:)  = 0.0
ZSNOWMELT(:,:)  = 0.0
ZSNOWTEMP(:,:)  = 0.0
ZMELTXS(:,:)    = 0.0
!
! 1. Determine amount of melt in each layer:
! ------------------------------------------
!
WHERE(PSNOWDZ(:,:) > 0.0)
!
! Total Liquid equivalent water content of snow (m):
!
   ZSNOWLWE(:,:) = PSNOWRHO(:,:)*PSNOWDZ(:,:)/XRHOLW
!
! Melt snow if excess energy and snow available:
! Phase change (J/m2)
!
   ZPHASE(:,:)  = MIN(PSNOWDZ(:,:)*PSCAP(:,:)*MAX(0.0, PSNOWTEMP(:,:) - XTT), &
                  MAX(0.0,ZSNOWLWE(:,:)-PSNOWLIQ(:,:))*XLMTT*XRHOLW          )  
!
!
! Update snow liq water content and temperature if melting:
! liquid water available for next layer from melting of snow
! which is assumed to be leaving the current layer (m):
!
   ZSNOWMELT(:,:) = ZPHASE(:,:)/(XLMTT*XRHOLW)
!
! Cool off snow layer temperature due to melt:
!
   ZSNOWTEMP(:,:) = PSNOWTEMP(:,:) - ZPHASE(:,:)/(PSCAP(:,:)*PSNOWDZ(:,:))
!
   PSNOWTEMP(:,:) = MIN(XTT, ZSNOWTEMP(:,:))
!
   ZMELTXS(:,:)   = (ZSNOWTEMP(:,:)-PSNOWTEMP(:,:))*PSCAP(:,:)*PSNOWDZ(:,:)
!
END WHERE
!
! Loss of snowpack depth: (m) and liquid equiv (m):
! Compression factor for melt loss: this decreases
! layer thickness and increases density thereby leaving
! total SWE constant. NOTE: All melt water
! in excess of the holding capacity is NOT used
! for compression, rather it decreases the layer
! thickness ONLY, causing a loss of SWE (outside
! of this routine).
!
ZWHOLDMAX(:,:)  = SNOW3LHOLD(PSNOWRHO,PSNOWDZ)
!
WHERE(PSNOWDZ > 0.0)
!
   ZCMPRSFACT(:,:) = (ZSNOWLWE(:,:)-MIN(PSNOWLIQ(:,:)+ZSNOWMELT(:,:),     &
                        ZWHOLDMAX(:,:)))/                                   &
                       (ZSNOWLWE(:,:)-MIN(PSNOWLIQ(:,:),ZWHOLDMAX(:,:)))  
!
   PSNOWDZ(:,:)    = PSNOWDZ(:,:)*ZCMPRSFACT(:,:)
   PSNOWRHO(:,:)   = ZSNOWLWE(:,:)*XRHOLW/PSNOWDZ(:,:)
!
! Make sure maximum density is not surpassed! If it is, lower the density
! and increase the snow thickness accordingly:

   ZCMPRSFACT(:,:) = MAX(XRHOLI, PSNOWRHO(:,:))/XRHOLI
   PSNOWDZ(:,:)    = PSNOWDZ(:,:)*ZCMPRSFACT(:,:)
   PSNOWRHO(:,:)   = ZSNOWLWE(:,:)*XRHOLW/PSNOWDZ(:,:)
!
!
! 2. Add snow melt to current snow liquid water content:
! ------------------------------------------------------
!
   PSNOWLIQ(:,:)   = PSNOWLIQ(:,:) + ZSNOWMELT(:,:)
!
END WHERE
!
! 3. Excess heat from melting and Total snowmelt diag
! ---------------------------------------------------
! use PMELTXS to warm underlying ground/vegetation layer to conserve energy
!
PMELTXS  (:) = 0.
DO JL = 1, INL
   DO JI = 1, INI
      PMELTXS  (JI) = PMELTXS  (JI) + ZMELTXS  (JI,JL) / PTSTEP  ! (W/m2)
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOW3LMELT',1,ZHOOK_HANDLE)
!
!
!
END SUBROUTINE SNOW3LMELT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LREFRZ(PTSTEP,PRR,PSNOWRHO,PSNOWTEMP,PSNOWDZ,PSNOWLIQ, &
                             PSNOWLIQ0,PTHRUFAL,PDELPHASE_SFC,PDELPHASE   )
!
!
!!    PURPOSE
!!    -------
!     Calculate any freezing/refreezing of liquid water in the snowpack.
!     Also, calculate liquid water transmission and snow runoff.
!     Water flow causes layer thickness reduction and can cause
!     rapid densification of a layer.
!
!
USE MODD_CSTS,     ONLY : XTT, XLMTT, XRHOLW
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                      :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)        :: PRR
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ, PSNOWTEMP, PSNOWLIQ, PSNOWRHO
!
REAL, DIMENSION(:,:), INTENT(IN   )   :: PSNOWLIQ0
!
REAL, DIMENSION(:), INTENT(INOUT)     :: PTHRUFAL ! kg/m2/s
!
REAL, DIMENSION(:), INTENT(OUT)       :: PDELPHASE_SFC, PDELPHASE ! kg/m2/s
!
!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ, JI
!
INTEGER                               :: INI
INTEGER                               :: INLVLS
!
REAL, DIMENSION(SIZE(PRR))            :: ZPCPXS, ZTOTWCAP, ZRAINFALL
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) ::   ZFLOWLIQ, ZWORK,     &
                                                        ZSNOWLIQ, ZSNOWRHO,  &
                                                        ZWHOLDMAX, ZSNOWDZ,  &
                                                        ZSNOWTEMP, ZSCAP,    &  
                                                        ZSNOWHEAT
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),0:SIZE(PSNOWRHO,2)):: ZFLOWLIQT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LREFRZ',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! --------------
!
PDELPHASE_SFC(:) = 0.0
PDELPHASE    (:) = 0.0
!
ZSNOWRHO(:,:)  = PSNOWRHO(:,:)
ZSNOWLIQ(:,:)  = PSNOWLIQ(:,:)
ZSNOWTEMP(:,:) = PSNOWTEMP(:,:)
INI            = SIZE(PSNOWDZ(:,:),1)
INLVLS         = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Refreeze due to heat transfer
!    -----------------------------
! Freeze liquid water in any layers which cooled due
! to heat transfer. First, update H and re-diagnose (update) T and Wl:
!
ZSCAP(:,:)     = SNOW3LSCAP(ZSNOWRHO)
!
ZSNOWHEAT(:,:) = PSNOWDZ(:,:)*( ZSCAP(:,:)*(ZSNOWTEMP(:,:)-XTT)        &
                 - XLMTT*ZSNOWRHO(:,:) ) + XLMTT*XRHOLW*ZSNOWLIQ(:,:)  
!
ZSNOWTEMP(:,:) = XTT + ( ((ZSNOWHEAT(:,:)/MAX(PSNOWDZ(:,:),XSNOWDZMIN))    &
                     + XLMTT*ZSNOWRHO(:,:))/ZSCAP(:,:) )
!
ZSNOWLIQ(:,:)  = MAX(0.0,ZSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*PSNOWDZ(:,:)/(XLMTT*XRHOLW)  
!
ZSNOWTEMP(:,:) = MIN(XTT,ZSNOWTEMP(:,:))
!
! 2. Reduce thickness due to snowmelt in excess of holding capacity
!    --------------------------------------------------------------
! Any water in excess of the
! Maximum holding space for liquid water
! amount is drained into next layer down.
! Loss of water due to snowmelt causes a reduction
! in snow layer mass by a reduction in thickness. Owing to a consistent
! decrease in both liq and thickness (and the fact that T=Tf when liquid present), 
! enthalpy is conserved.
!
ZWHOLDMAX(:,:) = SNOW3LHOLD(PSNOWRHO,PSNOWDZ)
!
ZFLOWLIQ(:,:)  = MAX(0.,ZSNOWLIQ(:,:)-ZWHOLDMAX(:,:))
!
ZSNOWLIQ(:,:)  = ZSNOWLIQ(:,:) - ZFLOWLIQ(:,:)
!
ZSNOWDZ(:,:)   = PSNOWDZ(:,:) - ZFLOWLIQ(:,:)*XRHOLW/ZSNOWRHO(:,:)
!
ZSNOWDZ(:,:)   = MAX(0.0, ZSNOWDZ(:,:))  ! to prevent possible very small
!                                          negative values (machine prescision
!                                          as snow vanishes
!
!
! 3. Liquid water flow: liquid precipitation and meltwater
!    -----------------------------------------------------
!
! Rainfall flowing into uppermost snow layer:
! If rainfall is excessive enough (or layers thin enough)
! it is simply routed directly to runoff: First calculate
! the total snowpack  available liquid water holding capacity:
!
ZTOTWCAP(:)   = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      ZTOTWCAP(JI) = ZTOTWCAP(JI) + ZWHOLDMAX(JI,JJ)
   ENDDO
ENDDO
!
! Rain entering snow (m):
!
ZRAINFALL(:)  = PRR(:)*PTSTEP/XRHOLW                ! rainfall (m)
!
ZFLOWLIQT(:,0)= MIN(ZRAINFALL(:),ZTOTWCAP(:))
!
! Rain assumed to directly pass through the pack to runoff (m):
!
ZPCPXS(:)     = ZRAINFALL(:) - ZFLOWLIQT(:,0)
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZFLOWLIQT(JI,JJ) = ZFLOWLIQ(JI,JJ)
   ENDDO
ENDDO
!
! Thickness is maintained during water through-flow,
! so that mass transfer is represented by
! density changes: NOTE a maximum density
! is assumed (XRHOSMAX_ES) so that all flow
! which would result in densities exceeding
! this limit are routed to next layer down.
! First test for saturation, then
! rout excess water down to next layer down
! and repeat calculation. Net gain in liquid (mass) is
! translated into a density increase:
!
ZFLOWLIQ(:,:)  = 0.0                ! clear this array for work
PSNOWLIQ(:,:)  = ZSNOWLIQ(:,:)      ! reset liquid water content
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSNOWLIQ(JI,JJ)  = ZSNOWLIQ(JI,JJ) + ZFLOWLIQT(JI,JJ-1)
      ZFLOWLIQ(JI,JJ)  = MAX(0.0, ZSNOWLIQ(JI,JJ)-ZWHOLDMAX(JI,JJ))
      ZSNOWLIQ(JI,JJ)  = ZSNOWLIQ(JI,JJ) - ZFLOWLIQ(JI,JJ)
      ZFLOWLIQT(JI,JJ) = ZFLOWLIQT(JI,JJ) + ZFLOWLIQ(JI,JJ)
   ENDDO
ENDDO
!
ZWORK    (:,:) = MAX(XSNOWDZMIN,ZSNOWDZ(:,:))
ZSNOWRHO (:,:) = ZSNOWRHO(:,:)+(ZSNOWLIQ(:,:)-PSNOWLIQ(:,:))*XRHOLW/ZWORK(:,:)  
ZSCAP    (:,:) = SNOW3LSCAP(ZSNOWRHO(:,:))
ZSNOWTEMP(:,:) = XTT +(((ZSNOWHEAT(:,:)/ZWORK(:,:))+XLMTT*ZSNOWRHO(:,:))/ZSCAP(:,:))
ZSNOWLIQ (:,:) = MAX(0.0,ZSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*ZSNOWDZ(:,:)/(XLMTT*XRHOLW)  
ZSNOWTEMP(:,:) = MIN(XTT,ZSNOWTEMP(:,:))
!
! Any remaining throughflow after freezing is available to
! the soil for infiltration or surface runoff (m).
! I.E. This is the amount of water leaving the snowpack:
! Rate water leaves the snowpack [kg/(m2 s)]:
!
PTHRUFAL(:)  = PTHRUFAL(:) + ZFLOWLIQT(:,INLVLS)*XRHOLW/PTSTEP
!
! Add excess rain (rain which flowed directly through the snow
! due to saturation):
!
PTHRUFAL(:)  = PTHRUFAL(:) + ZPCPXS(:)*XRHOLW/PTSTEP
!
! 4. Update thickness and density and any freezing:
!    ----------------------------------------------
!
PSNOWTEMP(:,:)= ZSNOWTEMP(:,:)
PSNOWDZ(:,:)  = ZSNOWDZ(:,:)
PSNOWRHO(:,:) = ZSNOWRHO(:,:)
PSNOWLIQ(:,:) = ZSNOWLIQ(:,:)
!
! 5. Compute phase change according to water liquid budget for diag :
!    ----------------------------------------------------------------
!
PDELPHASE_SFC(:) = (PSNOWLIQ(:,1)-PSNOWLIQ0(:,1)-ZFLOWLIQT(:,0)+ZFLOWLIQT(:,1))*XRHOLW/PTSTEP
!
PDELPHASE(:) = PTHRUFAL(:) - PRR(:)
DO JJ=1,INLVLS
   DO JI=1,INI
      PDELPHASE(JI) = PDELPHASE(JI) + (PSNOWLIQ(JI,JJ)-PSNOWLIQ0(JI,JJ))*XRHOLW/PTSTEP
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOW3LREFRZ',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LREFRZ
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LFLUX(PSNOWTEMP, PEXNS,PEXNA, PUSTAR2_IC,       &
                              PTSTEP,PALBT,PSW_RAD,PEMIST,PLWUPSNOW,  &
                              PLW_RAD,PLWNETSNOW,                     &
                              PTA,PSFCFRZ,PQA,PHPSNOW,                &
                              PSNOWTEMPO1,PRESTOREN,PCT,PRADSINK,     &
                              PQSAT,PDQSAT,PRSRA,                     &
                              PRN,PH,PGFLUX,PLES3L,PLEL3L,PEVAP,      &
                              PUSTAR,OSFCMELT                         )  
!
!
!!    PURPOSE
!!    -------
!     Calculate the surface fluxes (atmospheric/surface).
!     (Noilhan and Planton 1989; Noilhan and Mahfouf 1996)
!
!
USE MODD_CSTS,ONLY : XSTEFAN, XCPD, XLSTT, XLVTT, XTT
!
USE MODE_THERMOS
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWTEMPO1, PRESTOREN, PCT, &
                                         PRADSINK, PEXNS, PEXNA  
!
REAL, DIMENSION(:), INTENT(IN)      :: PALBT, PSW_RAD, PEMIST, PLW_RAD,        &
                                         PTA, PSFCFRZ, PQA,                    &
                                         PHPSNOW, PQSAT, PDQSAT, PRSRA,        &
                                         PUSTAR2_IC  
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRN, PH, PGFLUX, PLES3L, PLEL3L,        &
                                         PEVAP, PLWUPSNOW, PUSTAR,             &  
                                         PLWNETSNOW
!
LOGICAL, DIMENSION(:), INTENT(OUT)  :: OSFCMELT
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PCT))      :: ZEVAPC, ZLE, ZSNOWTEMP, ZSMSNOW, ZGFLUX,  &
                                   ZDELTAT, ZSNOWTO3  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LFLUX',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! --------------
!
ZSNOWTEMP(:)  = PSNOWTEMP(:)
ZLE(:)        = 0.0
ZSMSNOW(:)    = 0.0
ZGFLUX(:)     = 0.0
!
OSFCMELT(:)   = .FALSE.
!
ZSNOWTO3(:)   = PSNOWTEMPO1(:) ** 3  ! to save some CPU time, store this
!
! 1. Flux calculations when melt not occuring at surface (W/m2):
! --------------------------------------------------------------
!
!
ZDELTAT(:)   = PSNOWTEMP(:) - PSNOWTEMPO1(:)   ! surface T time change:
!
PLWUPSNOW(:) = PEMIST(:) * XSTEFAN * ZSNOWTO3(:)*( PSNOWTEMPO1(:) + 4.* ZDELTAT(:) )
!
PLWNETSNOW(:)= PEMIST(:) * PLW_RAD(:) -  PLWUPSNOW(:)
!
PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PLWNETSNOW(:)
!
PH(:)        = PRSRA(:) * XCPD * (PSNOWTEMP(:)/PEXNS(:) - PTA(:)/PEXNA(:))
!
ZEVAPC(:)    = PRSRA(:) * ( (PQSAT(:) - PQA(:)) + PDQSAT(:)*ZDELTAT(:) )
!
PLES3L(:)    = PSFCFRZ(:)     * XLSTT * ZEVAPC(:)
!
PLEL3L(:)    = (1.-PSFCFRZ(:))* XLVTT * ZEVAPC(:)
!
ZLE(:)       = PLES3L(:) + PLEL3L(:)
!
PGFLUX(:)    = PRN(:) - PH(:) - ZLE(:) + PHPSNOW(:)
!
!
! 2. Initial melt adjustment
! --------------------------
! If energy avalabile to melt snow, then recalculate fluxes
! at the freezing point and add residual heat to layer
! average heat.
!
! A) If temperature change is > 0 and passes freezing point this timestep,
!    then recalculate fluxes at freezing point and excess energy
!    will be used outside of this routine to change snow heat content:
!
WHERE (PSNOWTEMP(:) > XTT .AND. PSNOWTEMPO1(:) < XTT)
!
   OSFCMELT(:)= .TRUE.
!
   ZDELTAT(:)    = XTT - PSNOWTEMPO1(:)
!
   PLWUPSNOW(:) = PEMIST(:) * XSTEFAN * ZSNOWTO3(:)*( PSNOWTEMPO1(:) + 4.* ZDELTAT(:) ) 
!
   PLWNETSNOW(:)= PEMIST(:) * PLW_RAD(:) -  PLWUPSNOW(:)
!
   PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PLWNETSNOW(:)
!
   PH(:)        = PRSRA(:) * XCPD * (XTT/PEXNS(:) - PTA(:)/PEXNA(:))   
!
   ZEVAPC(:)    = PRSRA(:) * ( (PQSAT(:) - PQA(:)) + PDQSAT(:)*ZDELTAT(:) )
!
   PLES3L(:)    = PSFCFRZ(:)     * XLSTT * ZEVAPC(:)
!
   PLEL3L(:)    = (1.-PSFCFRZ(:))* XLVTT * ZEVAPC(:)
!
   ZLE(:)       = PLES3L(:) + PLEL3L(:)
!
   ZGFLUX(:)    = PRN(:) - PH(:) - ZLE(:) + PHPSNOW(:)
!
   ZSMSNOW(:)   = PGFLUX(:) - ZGFLUX(:)
!
   PGFLUX(:)  = ZGFLUX(:)
!
! This will be used to change heat content of snow:
!
   ZSNOWTEMP(:) = PSNOWTEMP(:) - ZSMSNOW(:)*PTSTEP*PCT(:)
!
END WHERE
!
! 3. Ongoing melt adjustment: explicit solution
! ---------------------------------------------
!    If temperature change is 0 and at freezing point this timestep,
!    then recalculate fluxes and surface temperature *explicitly*
!    as this is *exact* for snow at freezing point (Brun, Martin)
!
WHERE(PSNOWTEMP(:) > XTT .AND. PSNOWTEMPO1(:) >= XTT)
!
   OSFCMELT(:)  = .TRUE.   
!
   PLWUPSNOW(:) = PEMIST(:) * XSTEFAN * (XTT ** 4) 
!
   PLWNETSNOW(:)= PEMIST(:) * PLW_RAD(:) -  PLWUPSNOW(:)
!
   PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PLWNETSNOW(:)
!
   PH(:)        = PRSRA(:) * XCPD * (XTT/PEXNS(:) - PTA(:)/PEXNA(:))
!
   ZEVAPC(:)    = PRSRA(:) * (PQSAT(:) - PQA(:))
!
   PLES3L(:)    = PSFCFRZ(:)     * XLSTT * ZEVAPC(:)
!
   PLEL3L(:)    = (1.-PSFCFRZ(:))* XLVTT * ZEVAPC(:)
!
   ZLE(:)       = PLES3L(:) + PLEL3L(:)
!
   PGFLUX(:)    = PRN(:) - PH(:) - ZLE(:) + PHPSNOW(:)
!
   ZSNOWTEMP(:) = XTT + PTSTEP*PCT(:)*(PGFLUX(:) + PRADSINK(:) - PRESTOREN(:))
!
END WHERE
!
! 4. Update surface temperature:
! ------------------------------
!
PSNOWTEMP(:) = ZSNOWTEMP(:)
!
! 5. Final evaporative flux (kg/m2/s)
!
PEVAP(:)     = ZEVAPC(:)
!
! 6. Friction velocity
! --------------------
!
PUSTAR(:) = SQRT(PUSTAR2_IC(:)) 
!
IF (LHOOK) CALL DR_HOOK('SNOW3LFLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LFLUX
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LEVAPN(PPSN3L,PLES3L,PLEL3L,PTSTEP,PPERMSNOWFRAC,&
                             PSNOWTEMP,PSNOWRHO,PSNOWDZ,PSNOWLIQ,PTA,  &
                             PLVTT,PLSTT,PSNOWHEAT,PSOILCOR            )
!
!
!!    PURPOSE
!!    -------
!     Remove mass from uppermost snow layer in response to
!     evaporation (liquid) and sublimation.
!
!
USE MODD_CSTS,     ONLY : XRHOLW, XLMTT, XCI, XTT
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PPERMSNOWFRAC
!
REAL, DIMENSION(:), INTENT(IN)      :: PPSN3L 
!
REAL, DIMENSION(:), INTENT(IN)      :: PLES3L, PLEL3L   ! (W/m2)
!
REAL, DIMENSION(:), INTENT(IN)      :: PTA, PLVTT, PLSTT
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWRHO, PSNOWLIQ, &
                                       PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSOILCOR
!
!*      0.2    declarations of local variables
!
INTEGER                             :: INI, INLVLS, JJ, JI
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWEVAPS, ZSNOWEVAP, ZSNOWEVAPX,  &
                                       ZSNOWDZ, ZSNOWHEAT, ZSCAP, ZSNOWTEMP
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZXSE, ZISNOWD, ZTDIF

!*      0.3    declarations of local parameters

REAL, PARAMETER                     :: ZSNOWDEMIN = 1.E-4 ! m
REAL, PARAMETER                     :: ZTDIF1      = 15.  ! K To prevent a possible
                                                          !   decoupling of sfc snow T
                                                          !   when vanishingly thin, impose
                                                          !   this max T-diff based on obs...
                                                          !   between Ta-Ts
REAL, PARAMETER                     :: ZTDIF2      = 50.  ! Same but for Greenland/Antractic
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPN',0,ZHOOK_HANDLE)
!
PSOILCOR(:)    = 0.0
!
ZSNOWEVAPS(:)  = 0.0
ZSNOWEVAP(:)   = 0.0
ZSNOWEVAPX(:)  = 0.0
ZSNOWDZ(:)     = 0.0
ZSCAP(:)       = 0.0
ZSNOWHEAT(:)   = 0.0
ZSNOWTEMP(:)   = 0.0
ZXSE(:)        = 0.0
!
INI            = SIZE(PSNOWDZ(:,:),1)
INLVLS         = SIZE(PSNOWDZ(:,:),2)
!
ZTDIF(:)=(1.0-PPERMSNOWFRAC(:))*ZTDIF1+PPERMSNOWFRAC(:)*ZTDIF2
!
! 1. Evaporation of snow liquid water
! -----------------------------------
! Evaporation reduces liq water equivalent content
! of snowpack either by reducing the snow density
! (evaporation) or the layer thickness (sublimation).
! Condensation does the reverse.
!
! Multiply evaporation components by snow fraction
! to be consistent with fluxes from the snow covered
! fraction of grid box
!
WHERE(PSNOWDZ(:,1) > 0.0)
!
! Evaporation:
! Reduce density and liquid water content:
!
   ZSNOWEVAP(:)   = PPSN3L(:)*PLEL3L(:)*PTSTEP/(PLVTT(:)*XRHOLW)
   ZSNOWEVAPX(:)  = MIN(ZSNOWEVAP(:),PSNOWLIQ(:))
!
!  This should not change enthalpy (since already accounted for 
!  in sfc e budget), so make density change insuring constant enthalpy:
!
   PSNOWLIQ(:)    = PSNOWLIQ(:) - ZSNOWEVAPX(:)
   PSNOWRHO(:) = (PSNOWHEAT(:,1)-XLMTT*XRHOLW*PSNOWLIQ(:))/ &
                 (PSNOWDZ(:,1)*(XCI*(PSNOWTEMP(:)-XTT)-XLMTT))
!
! Budget check: If density drops below minimum, then extract the
! corresponding water mass from soil (for vanishingly thin snow covers):
!
   PSOILCOR(:)    = MAX(0.0,XRHOSMIN_ES-PSNOWRHO(:))*PSNOWDZ(:,1)/PTSTEP
   PSNOWRHO(:)    = MAX(XRHOSMIN_ES,PSNOWRHO(:))
!
END WHERE
!
! 2. Update heat capacity:
! ---------------------------
!
ZSCAP(:) = SNOW3LSCAP(PSNOWRHO(:))
!
!
! 3. Sublimation of snow ice
! ---------------------------
!
WHERE(PSNOWDZ(:,1) > 0.0)

! Budget check: as last traces of liquid in snow evaporates, it is possible
! evaporation could exceed liquid present in snow. If this is the case,
! remove residual mass from solid snow in order to maintain high
! accuracy water balance:

   ZSNOWEVAPX(:)  = MAX(0.0, ZSNOWEVAP(:) - ZSNOWEVAPX(:))
   ZSNOWDZ(:)     = PSNOWDZ(:,1) - ZSNOWEVAPX(:)*XRHOLW/PSNOWRHO(:)
   PSNOWDZ(:,1)   = MAX(0.0, ZSNOWDZ(:))
   PSOILCOR(:)    = PSOILCOR(:) + MAX(0.0,-ZSNOWDZ(:))*PSNOWRHO(:)/PTSTEP
   
! Sublimation: Reduce layer thickness and total snow depth
! if sublimation: add to correction term if potential
! sublimation exceeds available snow cover.
!
   ZSNOWEVAPS(:)  = PPSN3L(:)*PLES3L(:)*PTSTEP/(PLSTT(:)*PSNOWRHO(:))
   ZSNOWDZ(:)     = PSNOWDZ(:,1) - ZSNOWEVAPS(:)
   PSNOWDZ(:,1)   = MAX(0.0, ZSNOWDZ(:))
   PSOILCOR(:)    = PSOILCOR(:) + MAX(0.0,-ZSNOWDZ(:))*PSNOWRHO(:)/PTSTEP

! Note that the effect on enthalpy is already accounted for via the surface energy
! budget, so the change in enthalpy
! here owing to a mass loss must be offset by possible adjustments to T and Wl
! to conserve Enthalpy.

   PSNOWTEMP(:)   = XTT + ( ((PSNOWHEAT(:,1)/MAX(XSNOWDZMIN,PSNOWDZ(:,1)))   &
                    + XLMTT*PSNOWRHO(:))/ZSCAP(:) )

   PSNOWLIQ(:)    = MAX(0.0,PSNOWTEMP(:)-XTT)*ZSCAP(:)*                  &
                    PSNOWDZ(:,1)/(XLMTT*XRHOLW)  

   PSNOWTEMP(:)   = MIN(XTT,PSNOWTEMP(:)) 

! For vanishingly thin snow, a decoupling can occur between forcing level T
! and snow sfc T (as snowpack vanishes owing to sublimation), so a simple limit
! imposed on this gradient:

! Bertrand : I am not sure that it is physical... to be tested without...

   PSNOWTEMP(:)   = MAX(MIN(XTT,PTA(:)-ZTDIF(:)), PSNOWTEMP(:))

! Update surface enthalpy:

   ZSNOWHEAT(:)   = PSNOWHEAT(:,1)
   PSNOWHEAT(:,1) = PSNOWDZ(:,1)*( ZSCAP(:)*(PSNOWTEMP(:)-XTT)             &
                    - XLMTT*PSNOWRHO(:) ) + XLMTT*XRHOLW*PSNOWLIQ(:) 

   ZXSE(:)        = PSNOWHEAT(:,1) - ZSNOWHEAT(:) ! excess cooling (removed from sfc)
 
END WHERE
!
! If the sfc T-Gradient was limited (and thus the enthalpy in the uppermost layer changed),
! distribute this energy correction (cooling) in the layers below (thickness-based weighting)
! to conserve total snowpack enthalpy. Normally this only occurs for this snowpacks,
! but it can rarely occur otherwise also.
!
ZISNOWD(:) = 0.
DO JJ=2,INLVLS
   DO JI=1,INI
      ZISNOWD(JI) = ZISNOWD(JI) + PSNOWDZ(JI,JJ)        ! m
   ENDDO
END DO
ZISNOWD(:)        = ZXSE(:)/MAX(ZISNOWD(:),XSNOWDMIN)  ! J kg-1 m-1 
!
DO JJ=2,INLVLS
   DO JI=1,INI
      PSNOWHEAT(JI,JJ) = PSNOWHEAT(JI,JJ) - PSNOWDZ(JI,JJ)*ZISNOWD(JI)
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LEVAPN
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LGONE(PTSTEP,PLEL3L,PLES3L,PSNOWRHO,                 &
                   PSNOWHEAT,PRADSINK,PEVAPCOR,PTHRUFAL,PGRNDFLUX,   &
                   PGFLUXSNOW,PGRNDFLUXO,PSNOWDZ,PSNOWLIQ,PSNOWTEMP, &
                   PLVTT,PLSTT,PRADXS,PSNOWMELT,PSNOWGONE_DELTA      )  
!
!
!
!!    PURPOSE
!!    -------
!     Account for the case when the last trace of snow melts
!     during a time step: ensure mass and heat balance of
!     snow AND underlying surface.
!
!
USE MODD_CSTS,        ONLY : XTT, XLSTT, XLVTT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PLEL3L, PLES3L, PGFLUXSNOW, &
                                       PRADSINK, PGRNDFLUXO,       &
                                       PLVTT, PLSTT
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWHEAT
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PGRNDFLUX, PRADXS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWLIQ, PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL   ! melt water [kg/(m2 s)]
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWMELT  ! total melt water [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)     :: PEVAPCOR   ! [kg/(m2 s)]
!                                      PEVAPCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance.
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWGONE_DELTA
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWHEATC, ZSNOW
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LGONE',0,ZHOOK_HANDLE)
!
INI                   = SIZE(PSNOWDZ(:,:),1)
INLVLS                = SIZE(PSNOWDZ(:,:),2)
!
PEVAPCOR(:)           = 0.0
PSNOWMELT(:)          = 0.0
!
ZSNOWHEATC(:)         = 0.
ZSNOW(:)              = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSNOWHEATC(JI) = ZSNOWHEATC(JI) + PSNOWHEAT(JI,JJ) ! total heat content (J m-2)
      ZSNOW(JI)      = ZSNOW(JI)      + PSNOWDZ(JI,JJ)   ! total snow depth (m)
   ENDDO
ENDDO
PSNOWGONE_DELTA(:)    = 1.0
!
! 1. Simple test to see if snow vanishes:
! ---------------------------------------
! If so, set thicknesses (and therefore mass and heat) and liquid content
! to zero, and adjust fluxes of water, evaporation and heat into underlying
! surface. Note, test with flux computed *before* correction since this represents
! actual inflow of heat from below (as heat content correction owing to a corrected
! flux has not yet been done: here we compare to pre-corrected heat content).
!
WHERE(PGFLUXSNOW(:) + PRADSINK(:) >= (-ZSNOWHEATC(:)/PTSTEP) )
   PGRNDFLUX(:)       = PGFLUXSNOW(:) + (ZSNOWHEATC(:)/PTSTEP)
   PEVAPCOR(:)        = (PLEL3L(:)/PLVTT(:)) + (PLES3L(:)/PLSTT(:))
   PRADXS(:)          = 0.0
   PSNOWGONE_DELTA(:) = 0.0          ! FLAG...if=0 then snow vanishes, else=1
END WHERE
!
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOWMELT(JI) = PSNOWMELT(JI) + (1.0-PSNOWGONE_DELTA(JI))*PSNOWRHO(JI,JJ)*PSNOWDZ(JI,JJ)/PTSTEP
   ENDDO
END DO
!
PTHRUFAL(:) = PSNOWMELT(:)
!
! 2. Final update of snow state
! -----------------------------
! (either still present or not):
!
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOWDZ  (JI,JJ) =                                 PSNOWDZ  (JI,JJ)*PSNOWGONE_DELTA(JI)
      PSNOWLIQ (JI,JJ) =                                 PSNOWLIQ (JI,JJ)*PSNOWGONE_DELTA(JI)
      PSNOWTEMP(JI,JJ) = (1.0-PSNOWGONE_DELTA(JI))*XTT + PSNOWTEMP(JI,JJ)*PSNOWGONE_DELTA(JI)
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOW3LGONE',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOW3LGONE
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LEVAPGONE(PSNOWHEAT,PSNOWDZ,PSNOWRHO,PSNOWTEMP,PSNOWLIQ)
!
!!    PURPOSE
!!    -------
!
!     If all snow in uppermost layer evaporates/sublimates, re-distribute
!     grid (below assumes very thin snowpacks so layer-thicknesses are
!     constant).
!
!
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XRHOSMAX_ES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWRHO   ! snow density profile                (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ    ! snow layer thickness profile        (m)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWHEAT  ! snow heat content/enthalpy          (J/m2)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWTEMP  ! snow temperature profile            (K)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWLIQ   ! snow liquid water profile           (m)
!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ, JI
!
INTEGER                               :: INI
INTEGER                               :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOWHEAT_1D ! total heat content                (J/m2)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOW        ! total snow depth                  (m)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZMASS        ! total mass
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSCAP  ! Snow layer heat capacity          (J/K/m3)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Initialize:
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPGONE',0,ZHOOK_HANDLE)
INI             = SIZE(PSNOWDZ,1)
INLVLS          = SIZE(PSNOWDZ,2)
!
! First, determine where uppermost snow layer has completely
! evaporated/sublimated (as it becomes thin):
!
ZSNOWHEAT_1D(:) = 0.
ZSNOW(:)        = 0.
ZMASS(:)        = 0.
!
ZSCAP(:,:) = SNOW3LSCAP(PSNOWRHO(:,:))
!
DO JJ=2,INLVLS
   DO JI=1,INI
      IF(PSNOWDZ(JI,1) == 0.0)THEN
         ZSNOWHEAT_1D(JI) = ZSNOWHEAT_1D(JI) + XLMTT*XRHOLW*PSNOWLIQ(JI,JJ)  &
                          + PSNOWDZ(JI,JJ)*(ZSCAP(JI,JJ)*(PSNOWTEMP(JI,JJ)-XTT) &
                          - XLMTT*PSNOWRHO(JI,JJ) )           
         ZSNOW       (JI) = ZSNOW(JI) + PSNOWDZ(JI,JJ)
         ZMASS       (JI) = ZMASS(JI) + PSNOWDZ(JI,JJ)*PSNOWRHO(JI,JJ)    
       ENDIF
    ENDDO
ENDDO
!
! Where uppermost snow layer has vanished, redistribute vertical
! snow mass and heat profiles (and associated quantities):
!
DO JJ=1,INLVLS
   DO JI=1,INI
      IF(ZSNOW(JI)/= 0.0)THEN
        ZSNOW    (JI)    = MAX(0.5*XSNOWDMIN,ZSNOW(JI))
        PSNOWDZ  (JI,JJ) = ZSNOW(JI)/REAL(INLVLS)
        PSNOWHEAT(JI,JJ) = ZSNOWHEAT_1D(JI)/REAL(INLVLS)
        PSNOWRHO (JI,JJ) = ZMASS (JI)/ZSNOW(JI)
      ENDIF
    ENDDO
ENDDO        
!
ZSCAP(:,:) = SNOW3LSCAP(PSNOWRHO(:,:))
!
DO JJ=1,INLVLS
   DO JI=1,INI
      IF(ZSNOW(JI)/= 0.0)THEN
        PSNOWTEMP(JI,JJ) = XTT + ( ((PSNOWHEAT(JI,JJ)/PSNOWDZ(JI,JJ))   &
                               + XLMTT*PSNOWRHO(JI,JJ))/ZSCAP(JI,JJ) )  
        PSNOWLIQ (JI,JJ) = MAX(0.0,PSNOWTEMP(JI,JJ)-XTT)*ZSCAP(JI,JJ)*  &
                                   PSNOWDZ(JI,JJ)/(XLMTT*XRHOLW)  
        PSNOWTEMP(JI,JJ) = MIN(XTT,PSNOWTEMP(JI,JJ))
      ENDIF
    ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPGONE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW3LEVAPGONE
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LEBUDMEB(PTSTEP, PSNOWDZMIN,                        &
           PTS, PSNOWDZ1, PSNOWDZ2, PSCOND1, PSCOND2, PSCAP,        &
           PSWNETSNOWS, PLWNETSNOW,                                 &
           PHSNOW, PLES3L, PLEL3L, PHPSNOW,                         &
           PCT, PTSTERM1, PTSTERM2, PGFLUXSNOW                      )                          
!
!!    PURPOSE
!!    -------
!     Calculate surface energy budget with surface fluxes imposed.
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,    INTENT(IN)               :: PTSTEP, PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)    :: PTS, PSNOWDZ1, PSNOWDZ2, PSCOND1, PSCOND2, PSCAP,  &
                                     PHSNOW, PLES3L, PLEL3L, PHPSNOW,                   &
                                     PSWNETSNOWS, PLWNETSNOW
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCT, PTSTERM1, PTSTERM2, PGFLUXSNOW
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTS))        :: ZSCONDA, ZA, ZB, ZC,             &
                                     ZSNOWDZM1, ZSNOWDZM2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEBUDMEB',0,ZHOOK_HANDLE)
!
! Calculate surface energy budget components:
! ---------------------------------------------------------
! To prevent numerical difficulties for very thin snow
! layers, limit the grid "thinness": this is important as
! layers become vanishing thin:
!
ZSNOWDZM1(:)  = MAX(PSNOWDZ1(:), PSNOWDZMIN)
ZSNOWDZM2(:)  = MAX(PSNOWDZ2(:), PSNOWDZMIN)
!
! Surface thermal inertia:
!
PCT(:)        = 1.0/(PSCAP(:)*ZSNOWDZM1(:))
!
! Surface fluxes entering the snowpack (radiative and turbulent):
!
PGFLUXSNOW(:) = PSWNETSNOWS(:) + PLWNETSNOW(:) - PHSNOW(:) - PLES3L(:) - PLEL3L(:)
!
! Thermal conductivity between uppermost and lower snow layers:
!
ZSCONDA(:)    = (ZSNOWDZM1(:)+ZSNOWDZM2(:))/                           &
               ((ZSNOWDZM1(:)/PSCOND1(:)) + (ZSNOWDZM2(:)/PSCOND2(:)))
!
!
! Energy budget solution terms (with surface flux imposed):
!
ZB(:)         = 1./PTSTEP 
!
ZA(:)         = ZB(:) + PCT(:)*(2*ZSCONDA(:)/(ZSNOWDZM2(:)+ZSNOWDZM1(:)))
!
ZC(:)         = PCT(:)*( PGFLUXSNOW(:) + PHPSNOW(:) )
!
! Coefficients needed for implicit solution
! of linearized surface energy budget:
!
PTSTERM2(:)   = 2*ZSCONDA(:)*PCT(:)/(ZA(:)*(ZSNOWDZM2(:)+ZSNOWDZM1(:)))
!
PTSTERM1(:)   = (PTS(:)*ZB(:) + ZC(:))/ZA(:)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOW3LEBUDMEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW3LEBUDMEB
!####################################################################
!####################################################################
!####################################################################
!
!
!
END SUBROUTINE SNOW3L
