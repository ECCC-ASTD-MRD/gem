!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
module snow_svs2_mod
  implicit none
  public
contains
SUBROUTINE SNOW_SVS2(  PSNOWSWE,PSNOWTEMP, PSNOWLIQ,PSNOWRHO,PSNOWALB,  &
                      PSNOWAGE, PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST,     & 
                      PTSTEP,PTG, PCT, PSOILHCAPZ, & 
                      PSOILCONDZ, PPS, PTA,PSW_RAD, PQA, PVMOD,  PWIND_DRIFT,    &
                      PLW_RAD,PRR, PSR, PUNLOAD, PRESA_SN, PRHOA, PUREF,         &
                      PZREF, PALB, PD_G, PDZG,         &
                      PTHRUFAL, PGRNDFLUX,PRNSNOW,PHSNOW, PGFLUXSNOW,        &
                      PSWNETSNOW, PLWNETSNOW,PSUBLDRIFT,                     &
                      PHPSNOW, PPSN,PZ0NAT, PZ0EFF, PZ0HNAT, &
                      PLES3L, PLEL3L, PEVAP, &
                      PZENITH, PLAT, PLON, PFOREST,      &
                      PSNOWTYPE, PHVEGAPOL, PAGINGCOEF,  &
                      N ,  NSOIL)
!     ######################################################################################
!
!!****  *SNOW_SVS2*  
!!
!!    PURPOSE
!!    -------
!
!     This routine acts as an interface between the snowpack schames SNOWES & Crocus
!     and the main code of SVS2
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
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!	V. Vionnet           * ECCC *
!!      based on interface routine snow3L_isba.F90 from SURFEX
!!
!!    MODIFICATIONS
!!    -------------

!-------------------------------------------------------------------------------
!

!
USE MODD_CSTS,       ONLY : XTT, XPI, XDAY, XLMTT, XLSTT, XLVTT,XCI,XP00, XRD, XCPD,  &
                                   XRHOLW,XSTEFAN, XG
USE MODD_SNOW_PAR,   ONLY : XRHOSMAX_ES, XSNOWDMIN, XRHOSMIN_ES, XEMISSN 
USE MODD_SURF_PAR,   ONLY : XUNDEF 
USE MODD_PREP_SNOW,   ONLY : NIMPUR
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MU_JDATE_MOD,    ONLY: MU_JS2YMDHMS
!
USE MODI_SNOW3L
USE MODI_SNOWCRO
USE MODI_SNOWCRO_DIAG
!
! Namelist containing the options of the snowpack schemes. 
use sfc_options
!
IMPLICIT NONE
!
!
!
!*      0.1    declarations of arguments
!
!

!
!VV TYPE(DATE_TIME), INTENT(IN)         :: TPTIME     ! current date and time
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
!
!INTEGER N, NSL, NSOIL
INTEGER N,  NSOIL
!                                      N         = Size of "TRNCH" i.e.,of row passed to SVS
!                                      NSL       = Number of snow layers
!                                      NSOIL     = Number of soil layers
!
REAL, DIMENSION(N), INTENT(IN)    :: PFOREST ! = 1 for High vegeation or 0 for low vegeation 
!                                              ! Use to modify the effect of the snowdrift routine forested areas. 
!
REAL, DIMENSION(N), INTENT(INOUT)   :: PSNOWALB
!                                      PSNOWALB = Prognostic surface snow albedo
!                                                 (does not include anything but
!                                                 the actual snow cover)
!
REAL, DIMENSION(N,NSL), INTENT(INOUT) :: PSNOWTEMP, PSNOWRHO, PSNOWSWE,PSNOWLIQ
!                                      PSNOWTEMP = Snow layer(s) temperature (K)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)

REAL, DIMENSION(N,NSL), INTENT(INOUT) :: PSNOWAGE    ! Snow grain age
!
REAL, DIMENSION(N,NSL), INTENT(INOUT) :: PSNOWDIAMOPT  ! Snow grain characteristic
                                                     !  PSNOWDIAMOPT is the optical diameter
REAL, DIMENSION(N,NSL), INTENT(INOUT) :: PSNOWSPHERI  ! Snow grain characteristic
                                                     !  PSNOWSPHERI is the sphericity
REAL, DIMENSION(N,NSL), INTENT(INOUT) :: PSNOWHIST   ! Historical variables in Crocus metamorphism scheme

REAL, DIMENSION(N,NSL), INTENT(OUT) :: PSNOWTYPE   ! Diagnostic snow grain type in Crocus 

!
REAL, DIMENSION(N), INTENT(INOUT)   :: PRNSNOW, PHSNOW,PHPSNOW 
!                                      PHPSNOW       = heat release from rainfall (W/m2)
!                                      PRNSNOW       = net radiative flux from snow (W/m2)
!                                      PHSNOW        = sensible heat flux from snow (W/m2)
!
REAL, DIMENSION(N), INTENT(OUT)     :: PGFLUXSNOW
!                                      PGFLUXSNOW    = net heat flux from snow (W/m2)
!
REAL, DIMENSION(N), INTENT(IN)      :: PPSN,PLAT, PLON
!                                      PPSN     = Snow cover fraction (total) 
!                                      PLAT     = Latitude (in rad) 
!                                      PLON     = Longitude (in rad)
!
REAL, DIMENSION(N), INTENT(IN)      :: PZ0NAT, PZ0EFF, PZ0HNAT
!                                      PZ0EFF    = roughness length for momentum 
!                                      PZ0NAT    = grid box average roughness length
!                                      PZ0HNAT   = grid box average roughness length
!
!
REAL, DIMENSION(N,NSOIL), INTENT(INOUT) :: PTG
!                                      PTG       = Soil temperature profile (K)
!
REAL, DIMENSION(N,NSOIL), INTENT(IN)    :: PSOILHCAPZ, PD_G, PDZG
REAL, DIMENSION(N),   INTENT(IN)    :: PCT, PSOILCONDZ  
!                                      PD_G      = Depth to bottom of each soil layer (m)
!                                      PDZG      = Soil layer thicknesses (m)
!                                      PCT       = area-averaged surface heat capacity [(K m2)/J]
!                                      PSOILCONDZ= soil thermal conductivity (W m-1 K-1)
!                                      PSOILHCAPZ= soil heat capacity (J m-3 K-1)
!
REAL, DIMENSION(N), INTENT(IN)      :: PPS, PTA, PSW_RAD, PQA,                       &
                                       PVMOD, PWIND_DRIFT, PLW_RAD, PSR, PRR, PUNLOAD  
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = snow rate (SWE) [kg/(m2 s)]
!                                      PUNLOAD = unload rate (for snow below canopy) [kg/(m2 s)]
!                                      PTA     = atmospheric temperature at level za (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PWIND_DRIFT   = modulus of the wind under canopy if PVMOD is above canopy, otherwise PWIND_DRIFT = PVMOD (m/s)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
!
REAL, DIMENSION(N), INTENT(INOUT)      :: PRESA_SN    
!                                         PRESA_SN = Aerodynamic  resistance for snow under canopy (cf. Gouttevin et al. 2015), = 0 if open
!
REAL, DIMENSION(N), INTENT(IN)      :: PZREF, PUREF,PRHOA, PALB
!                                      PZREF     = reference height of the first
!                                                  atmospheric level
!                                      PUREF     = reference height of the wind
!                                      PRHOA     = air density
!                                      PALB      = soil/vegetation albedo
!
REAL, DIMENSION(N), INTENT(IN)      :: PHVEGAPOL   ! Mean polar vegetation height (m)
!
REAL, DIMENSION(N), INTENT(IN)      :: PAGINGCOEF   ! Snow aging coefficient (-)
!
REAL, DIMENSION(N), INTENT(INOUT)   :: PLES3L, PLEL3L, PEVAP, PGRNDFLUX
!                                      PLEL3L        = evaporation heat flux from snow (W/m2)
!                                      PLES3L        = sublimation (W/m2)
!                                      PEVAP         = total evaporative flux from snow (kg/m2/s)
!                                      PGRNDFLUX     = soil/snow interface heat flux (W/m2)
!
REAL, DIMENSION(N), INTENT(OUT)      :: PSWNETSNOW, PLWNETSNOW
!                                      PSWNETSNOW = net shortwave radiation entering top of snowpack 
!                                                  (W m-2) 
!                                      PLWNETSNOW = net longwave radiation entering top of snowpack 
!                                                  (W m-2) 
!
REAL, DIMENSION(N), INTENT(OUT)      :: PSUBLDRIFT
!                                      PSUBLDRIFT    = rate of mass loss due to blowing snow sublimation (kg/m2/s)
!                                                   (negative sign since it corresponds to a mass loss for the snow cover)
!
!
REAL, DIMENSION(N), INTENT(OUT)     :: PTHRUFAL
!                                      PTHRUFAL  = rate that liquid water leaves snow pack: 
!                                                  paritioned into soil infiltration/runoff 
!                                                  by ISBA [kg/(m2 s)]
!
!
! ajout_EB pour prendre en compte angle zenithal du soleil dans LRAD
! puis plus tard dans LALB
REAL, DIMENSION(N), INTENT(IN)      :: PZENITH    ! solar zenith angle




!  
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                     :: ZCHECK_TEMP = 150.0 
!                                      Limit to check suspicious low temperature (K)
!
INTEGER                             :: JWRK, JJ,JIMP ! Loop control
!
INTEGER                             :: INLVLS   ! maximum number of snow layers
INTEGER                             :: INLVLG   ! number of ground layers
INTEGER                             :: IBLOWSNW     ! number of blowing snow variables
!
REAL, DIMENSION(SIZE(PTA))          :: ZRRSNOW, ZSOILCOND, ZSNOW, ZSNOWFALL,  &
                                       ZSNOWABLAT_DELTA, ZSNOWSWE_1D, ZSNOWD, & 
                                       ZSNOWH, ZSNOWH1, ZGRNDFLUXN, ZPSN,     &
                                       ZSOILCOR, ZSNOWSWE_OUT, ZTHRUFAL,      &
                                       ZSNOW_MASS_BUDGET, ZWGHT, ZWORK, ZC2     
!                                      ZSOILCOND    = soil thermal conductivity [W/(m K)]
!                                      ZRRSNOW      = rain rate over snow [kg/(m2 s)]
!                                      ZSNOW        = snow depth (m) 
!                                      ZSNOWFALL    = minimum equivalent snow depth
!                                                     for snow falling during the
!                                                     current time step (m)
!                                      ZSNOWABLAT_DELTA = FLAG =1 if snow ablates completely
!                                                     during current time step, else=0
!                                      ZSNOWSWE_1D  = TOTAL snowpack SWE (kg m-2)
!                                      ZSNOWD       = snow depth
!                                      ZSNOWH       = snow total heat content (J m-2)
!                                      ZSNOWH1      = snow surface layer heat content (J m-2)
!                                      ZGRNDFLUXN   = corrected snow-ground flux (if snow fully ablated during timestep)
!                                      ZPSN         = snow fraction working array
!                                      ZSOILCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!                                      ZSNOW_MASS_BUDGET = snow water equivalent budget (kg/m2/s)
!                                      ZWGHT        = MEB surface layer weight for distributing energy
!                                                     between litter and ground layers for the case
!                                                     of total ablation during a timestep (-).
!                                      ZWORK        = local working variable (*)
!                                      ZC2          = sub-surface heat capacity [(K m2)/J]
!
!###########################################################################################
!########       Local variables for phasing with SVS
!###########################################################################################

!
TYPE(DATE_TIME)  :: TPTIME      ! current date and time

REAL, DIMENSION(SIZE(PTA))          :: ZEXNS,ZEXNA,ZDIRCOSZW,ZSLOPEDIR, ZAZIM
!                                      ZEXNS     = Exner function at surface
!                                      ZEXNA     = Exner function at lowest atmos level
!                                      ZDIRCOSZW = Cosinus of the angle between the 
!                                                  normal to the surface and the vertical
!                                      ZSLOPEDIR   = Slope direction
!                                      ZAZIM       =  azimuthal angle      (radian from North, clockwise)

LOGICAL     :: LMEB                ! Activate the Multi-Energy Budget Option present in SURFEX (set to FALSE since snow-vegeation
                                   ! interactions is handle differently in SVS)


CHARACTER(LEN=3)   :: CIMPLICIT_WIND   ! wind implicitation option
!                                                  ! 'OLD' = direct
!                                                  ! 'NEW' = Taylor serie, order 1 (default in SVS)
!
CHARACTER(LEN=3)    :: CSNOWZREF 
                                         ! reference height is constant or variable from the snow surface
                                         ! CSNOWZREF='CST' constant reference height from the snow surface  (Default in SVS-Cro)
                                         ! CSNOWZREF='VAR' variable reference height from the snow surface (i.e. constant from the ground)
CHARACTER(LEN=4)    ::  CSNOWFPAPPUS
                                       ! Option to force fresh snow characteristics when snowpappus is used (Not activated in SVS2)
                                       ! 'GM98' => forces fresh snow characteristics as if HSNOWDRIFT = 'NONE'
                                       ! 'VI13' => """" 'VI13'
                                       ! 'NONE' => no influence ( it is forced when snowpappus is not activated )
!
LOGICAL     :: LSNOW_ABS_ZENITH  ! Activate parametrization of solar absorption for polar regions (default FALSE in SVS)
                                 ! Need to be tested in the Arctic
LOGICAL     :: LATMORAD          ! Activate atmotartes scheme in the TARTES snowpack radiative scheme (default FALSE in SVS-Cro)
!
LOGICAL     :: LFOREST           ! Logical to determine when Crocus in SVS2 is called in a forested environement 
!

LOGICAL     :: LSNOWCOMPACT_BOOL,LSNOWMAK_BOOL,LSNOWTILLER,LSELF_PROD,& ! Logical for artificial snow making (default FALSE in SVS-Cro)
                LSNOWMAK_PROP

REAL, DIMENSION(SIZE(PTA))     ::      ZPEW_A_COEF, ZPEW_B_COEF,                   &
                                       ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF,      &
                                       ZPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient
!                                      PPEW_B_COEF = wind coefficient
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient   

REAL, DIMENSION(SIZE(PTA),SIZE(PSNOWSWE,2))     :: ZSNOWHEAT, ZSNOWDZ, ZSCAP
!                                      ZSNOWHEAT = Snow layer(s) heat content (J/m2)
!                                      ZSNOWDZ   = Snow layer(s) thickness (m)
!                                      ZSCAP      = Snow layer(s) heat capacity [J/(K m3)
!
REAL, DIMENSION(SIZE(PTA))     ::  ZLAT, ZLON
!					ZLAT = Latitude (deg)
!					ZLON = Longitude (deg)
!
REAL, DIMENSION(SIZE(PTA))       :: ZLVTT, ZLSTT ! = latent heats for hydrology
!
REAL, DIMENSION(SIZE(PTA))  :: ZFLSN_COR, ZEVAPCOR, ZSNOWHMASS, ZRESTOREN,               &
                                       ZDELHEATG, ZDELHEATG_SFC
!                                      ZFLSN_COR = soil/snow correction heat flux (W/m2) (not MEB)

!                                      ZEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!                                      ZSNOWHMASS = heat content change due to mass
!                                                   changes in snowpack (J/m2): for budget
!                                                   calculations only.
!                                      ZRESTOREN  = heat flux between the surface and sub-surface 
!                                                   snow layers (for energy budget diagnostics) (W/m2)
!                                      ZDELHEATG     = ground heat content change (diagnostic) (W/m2)
!                                                      note, modified if ground-snow flux adjusted
!                                      ZDELHEATG_SFC = ground heat content change in sfc only (diagnostic) (W/m2)
!                                                      note, modified if ground-snow flux adjusted

REAL, DIMENSION(SIZE(PTA))   :: ZSRSFC, ZRRSFC, ZSNOWSFCH, ZDELHEATN, ZDELHEATN_SFC,ZDELPHASEN, &
                                       ZDELPHASEN_SFC, ZMELTSTOT, ZSNREFREEZ
!                                      ZSRSFC = snow rate on soil/veg surface when SNOWES in use
!                                      ZRRSFC = rain rate on soil/veg surface when SNOWES in use
!                                      ZSNOWSFCH = snow surface layer pseudo-heating term owing to
!                                                  changes in grid thickness            (W m-2)
!                                      ZDELHEATN = total snow heat content change in the surface layer (W m-2)
!                                      ZDELHEATN_SFC = total snow heat content change during the timestep (W m-2)
!                                      ZDELPHASEN = latent heating due to snow melt/freeze  (W m-2)
!                                      ZDELPHASEN_SFC = latent heating due to surface layer snow melt/freeze  (W m-2)
!                                      ZMELTSTOT = snowmelt in the snowpack  [kg/(m2 s)]
!                                      ZSNREFREEZ = refreezing of water in the snowpack  [kg/(m2 s)]
!
REAL, DIMENSION(SIZE(PTA))  :: ZUSTARSNOW, ZCDSNOW, ZCHSNOW, ZRI,ZEMISNOW
!                                      ZCDSNOW    = drag coefficient for momentum over snow (-)
!                                      ZUSTARSNOW = friction velocity over snow (m/s)
!                                      ZCHSNOW    = drag coefficient for heat over snow (-)
!                                      ZRI        = Richardson number (-)
!                                      ZEMISNOW      = snow surface emissivity
!
REAL, DIMENSION(SIZE(PTA))        :: ZQS
!                                      ZQS = surface humidity (kg/kg)

REAL, DIMENSION(SIZE(PTA))        :: ZPA
!                                      ZQS = air pressure at the atmospheric forcing level (Pa)
!

INTEGER                           :: ZYY, ZMO, ZDD, ZHH, ZMN, ZSEC
!
!###########################################################################################
!########      END Definition of local variables for phasing with SVS
!###########################################################################################

!*      0.3    declarations of packed  variables
!
INTEGER                            :: ISIZE_SNOW ! number of points where computations are done
INTEGER, DIMENSION(SIZE(PTA))      :: NMASK      ! indices correspondance between arrays
!
LOGICAL, DIMENSION(SIZE(PTA))      :: LREMOVE_SNOW
!
REAL, DIMENSION(SIZE(PTA))         :: ZSWNET_NS
!
REAL, DIMENSION(SIZE(PPS))         :: ZSNOWMAK
!
!
! - - ---------------------------------------------------
!
!
!*       0.     Initialize variables:
!               ---------------------
!
ZSNOWHMASS(:)  = 0.0

ZSRSFC(:)      = PSR(:)         ! these are snow and rain rates passed to ISBA,
ZRRSFC(:)      = PRR(:)         ! so initialize here if SNOW3L not used:
!
ZFLSN_COR(:)   = 0.0
PTHRUFAL(:)    = 0.0
PSWNETSNOW(:)  = 0.0 
PLWNETSNOW(:)  = 0.0 
PSUBLDRIFT(:)    = 0.
ZEVAPCOR(:)    = 0.0
ZSOILCOR(:)    = 0.0
ZSNOW_MASS_BUDGET(:) = 0.0
ZQS(:)         = XUNDEF
!
ZMELTSTOT (:) = 0.0
ZSNREFREEZ(:) = 0.0
!
ZSNOW(:)       = 0.0
ZSNOWD(:)      = 0.0
ZGRNDFLUXN(:)  = 0.0
ZSNOWH(:)      = 0.0
ZSNOWH1(:)     = 0.0
ZSNOWSWE_1D(:) = 0.0
ZSNOWSWE_OUT(:)= 0.0
ZSOILCOND(:)   = 0.0
ZRRSNOW(:)     = 0.0
ZSNOWFALL(:)   = 0.0
ZSNOWABLAT_DELTA(:) = 0.0
ZWGHT(:)       = 0.0
ZWORK(:)       = 0.0
ZC2(:)         = PCT(:)
ZSNOWMAK(:) = 0.0
!
!
ZWGHT(:)       = 0.0
ZWORK(:)       = 0.0
ZC2(:)         = PCT(:)
!
!
INLVLS          = SIZE(PSNOWSWE(:,:),2)    
INLVLG          = MIN(SIZE(PD_G(:,:),2),SIZE(PTG(:,:),2)) 


!VV  IBLOWSNW       = SIZE(ZBLOWSNW(:,:),2)

!
!###########################################################################################
!########       Initialization of variables for phasing with SVS
!###########################################################################################
LMEB=.FALSE.

CIMPLICIT_WIND = 'NEW'
CSNOWZREF='CST'
CSNOWFPAPPUS='NONE'

LSNOW_ABS_ZENITH = .FALSE.
LATMORAD = .FALSE.

LSNOWCOMPACT_BOOL = .FALSE.
LSNOWMAK_BOOL = .FALSE.
LSNOWTILLER = .FALSE.
LSELF_PROD = .FALSE.
LSNOWMAK_PROP = .FALSE.


ZPEW_A_COEF(:) = 0.
ZPEW_B_COEF(:) = PVMOD(:)
ZPET_A_COEF(:) = 0.
! We shoudl use PPA in this formulation. Currently we use PPS
ZPET_B_COEF(:) =  PTA / (PPS/XP00)**(XRD/XCPD) 
ZPEQ_A_COEF(:) = 0.
ZPEQ_B_COEF(:) = PQA(:)  

! Thickness of snow layers (m)
ZSNOWDZ(:,:) = PSNOWSWE(:,:)/PSNOWRHO(:,:)

! WARNING  : calcul a reprendre en fonction des valeurs par defauts dans SVS
ZSCAP(:,:)     = PSNOWRHO(:,:)*XCI

ZSNOWHEAT(:,:)  = 0.
! Compute snow layers heat content (J/m2)
WHERE(PSNOWSWE(:,:)>0.) 
      ZSNOWHEAT(:,:) = ZSNOWDZ(:,:)*( ZSCAP(:,:)*(PSNOWTEMP(:,:)-XTT)        &
                   - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:)  
END WHERE

! Compute Exner function at the surface 
ZEXNS(:)   = (PPS(:)/XP00)**(XRD/XCPD)

! Compute pressure at the height of the atmospheric focring and get Exner function at this height
ZPA(:) = PPS(:)- PRHOA(:)*PZREF(:)*XG
ZEXNA(:)   = (ZPA(:)/XP00)**(XRD/XCPD)

! We assume so far a flat terrain.
ZDIRCOSZW(:) = 1.
ZSLOPEDIR(:) = 0.
ZAZIM(:) = 0. ! Azimuthal angle set to zero since it is not use internally in SVS2/Crocus so far
           ! To be modified when TARTES will be used

! Latent heats are initialized as in init_veg_pgdn.F90 from SURFEX
ZLSTT(:)   = XLSTT
ZLVTT(:)   = XLVTT

! No blowing snow variables are considerd in SVS at the moment
! Set the default value to 4 to avoid issues when running with Crocus
IBLOWSNW = 4
! Variables to keep track of heat exchanges when snow is remove are set to 0.
ZDELHEATG(:) = 0.
ZDELHEATG_SFC(:) = 0.

! Initialisation of date
CALL MU_JS2YMDHMS(JDATEO, ZYY, ZMO, ZDD, ZHH, ZMN, ZSEC)

! Initialize type date for Crocus
TPTIME%TDATE%YEAR = ZYY
TPTIME%TDATE%MONTH = ZMO
TPTIME%TDATE%DAY = ZDD
TPTIME%TIME = ZHH*3600.+ZMN*60+ZSEC

! Convert latitude and longitude from rad to deg
ZLAT(:) = PLAT(:)* 180./ACOS(-1.)
ZLON(:) = PLON(:)* 180./ACOS(-1.)

! Determine if snow scheme is called in a snow environement
LFOREST= ANY(PFOREST==1.)

!###########################################################################################
!########      End Initialization of variables for phasing with SVS
!###########################################################################################

!
  !
! If MEB activated, these values are input, else initialize here:
!
PGRNDFLUX(:)      = 0.0
PLES3L(:)         = 0.0
PLEL3L(:)         = 0.0    
PEVAP(:)          = 0.0
ZRI(:)            = XUNDEF  
ZEMISNOW(:) = XEMISSN    
PRNSNOW(:)    = 0.0
PHSNOW(:)     = 0.0
PGFLUXSNOW(:) = 0.0
PHPSNOW(:)    = 0.0
ZUSTARSNOW(:) = 0.0
ZCDSNOW(:)    = 0.0
ZCHSNOW(:)    = 0.0
!
! Initialize diagnostic to default value
PSNOWTYPE(:,:)= XUNDEF
!
!
ZSWNET_NS(:)      = 0.0


! - Snow and rain falling onto the ES grid space:
!
  ZSRSFC(:) = 0.0
!
  DO JJ=1,SIZE(PSR)
!VV    ZRRSNOW(JJ)        = PEK%XPSN(JJ)*PRR(JJ)
    ZRRSNOW(JJ)        = PRR(JJ)   ! VV: Modification compare to ISBA to close water budget in SVS
    ZRRSFC(JJ)    = PRR(JJ) - ZRRSNOW(JJ)
    ZSNOWFALL(JJ)      = PSR(JJ)*PTSTEP/XRHOSMAX_ES    ! maximum possible snowfall depth (m)
  ENDDO
!

!
! - Soil thermal conductivity
ZSOILCOND(:) = PSOILCONDZ(:)
!


!
! Calculate preliminary snow depth (m)


  ZSNOW(:)=0.
  ZSNOWH(:)=0.
  ZSNOWSWE_1D(:)=0.
  ZSNOWH1(:)    = ZSNOWHEAT(:,1)*PSNOWSWE(:,1)/PSNOWRHO(:,1) ! sfc layer only
!
  
  DO JWRK=1,SIZE(PSNOWSWE(:,:),2)
    DO JJ=1,SIZE(PSNOWSWE(:,:),1)
        ZSNOWSWE_1D(JJ)     = ZSNOWSWE_1D(JJ) + PSNOWSWE(JJ,JWRK)
        ZSNOW(JJ)           = ZSNOW(JJ)       + PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
        ZSNOWH(JJ)          = ZSNOWH(JJ)      + ZSNOWHEAT(JJ,JWRK)*PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
    END DO
  ENDDO
!
!
! ===============================================================
! === Packing: Only call snow model when there is snow on the surface
!              exceeding a minimum threshold OR if the equivalent
!              snow depth falling during the current time step exceeds 
!              this limit.
!
! counts the number of points where the computations will be made
!
!
  ISIZE_SNOW = 0
  NMASK(:) = 0

!
  DO JJ=1,SIZE(ZSNOW)
    IF (ZSNOW(JJ) >= XSNOWDMIN .OR. ZSNOWFALL(JJ) >= XSNOWDMIN) THEN
      ISIZE_SNOW = ISIZE_SNOW + 1
      NMASK(ISIZE_SNOW) = JJ
    ENDIF
  ENDDO
!  
  IF (ISIZE_SNOW>0) CALL CALL_MODEL(ISIZE_SNOW,INLVLS,INLVLG,IBLOWSNW,NMASK)
!
! ===============================================================
!
! Remove trace amounts of snow and reinitialize snow prognostic variables
! if snow cover is ablated.
! If MEB used, soil T already computed, therefore correct heating/cooling
! effect of updated snow-soil flux
!
  ZSNOWD(:) = 0.
  ZSNOWSWE_OUT(:) = 0.
  DO JWRK=1,SIZE(PSNOWSWE(:,:),2)
    DO JJ=1,SIZE(PSNOWSWE(:,:),1)
      ZSNOWD      (JJ) = ZSNOWD      (JJ) + PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
      ZSNOWSWE_OUT(JJ) = ZSNOWSWE_OUT(JJ) + PSNOWSWE(JJ,JWRK)
    ENDDO
  END DO
!
! Remove snow is activated when preexisting snow cover or new snow fall is ablated
!
  LREMOVE_SNOW(:)=((ZSNOWD(:)<XSNOWDMIN*1.1) .AND. (ZSNOW(:) >= XSNOWDMIN .OR. ZSNOWFALL(:) >= XSNOWDMIN)) 
!
!
!    To Conserve mass in ISBA without MEB, 
!    EVAP must be weighted by the snow fraction
!    in the calulation of THRUFAL
  ZPSN(:) = PPSN(:)
!
  ZSNOWABLAT_DELTA(:) = 0.0
  ZTHRUFAL        (:) = PTHRUFAL(:)
!
  WHERE(LREMOVE_SNOW(:))
    !
    ZSNOWSWE_OUT(:)     = 0.0
    PLES3L(:)           = MIN(PLES3L(:), XLSTT*(ZSNOWSWE_1D(:)/PTSTEP + PSR(:)))
    PLEL3L(:)           = 0.0
    PEVAP(:)            = PLES3L(:)/ZLSTT(:)
    PTHRUFAL(:)         = MAX(0.0, ZSNOWSWE_1D(:)/PTSTEP + PSR(:) - PEVAP(:)*ZPSN(:) + ZRRSNOW(:)) ! kg m-2 s-1
    ZTHRUFAL(:)         = MAX(0.0, ZSNOWSWE_1D(:)/PTSTEP + PSR(:) - PEVAP(:)         + ZRRSNOW(:)) ! kg m-2 s-1
    ZMELTSTOT(:)    = PTHRUFAL(:)
    ZSNREFREEZ(:)   = 0.0    
    !
    ZSRSFC(:)       = 0.0
    ZRRSFC(:)       = ZRRSFC(:)
    !
    ZSNOWABLAT_DELTA(:) = 1.0
    !
    !VV Change  for SVS
    !PSNOWALB(:)    = XUNDEF
    PSNOWALB(:)         = 0.1
    !
    ZEVAPCOR(:)         = 0.0
    ZSOILCOR(:)         = 0.0
    !
    PGFLUXSNOW(:)   = PRNSNOW(:) - PHSNOW(:) - PLES3L(:) - PLEL3L(:)
    ZSNOWHMASS(:)   = -PSR(:)*(XLMTT*PTSTEP)
    !
    ZRESTOREN(:)        = 0.0
    ZDELHEATN(:)        = -ZSNOWH(:) /PTSTEP
    ZDELHEATN_SFC(:)    = -ZSNOWH1(:)/PTSTEP
    ZDELPHASEN(:)    = 0.0
    ZDELPHASEN_SFC(:)= 0.0     

    ZSNOWSFCH(:)        = ZDELHEATN_SFC(:) - (ZSWNET_NS(:) + PLWNETSNOW(:)    &
                           - PHSNOW(:) - PLES3L(:) - PLEL3L(:)) + ZRESTOREN(:)      &
                           - ZSNOWHMASS(:)/PTSTEP 
    ZGRNDFLUXN(:)       = (ZSNOWH(:)+ZSNOWHMASS(:))/PTSTEP + PGFLUXSNOW(:)
    ZWORK(:)            = PTSTEP * ZPSN(:) * (ZGRNDFLUXN(:) - PGRNDFLUX(:) - ZFLSN_COR(:))
    PTG(:,1)            = PTG(:,1) + ZWORK(:)*(1.-ZWGHT(:))*PCT(:)
    PTG(:,2)            = PTG(:,2) + ZWORK(:)*    ZWGHT(:) *ZC2(:)
    ZWORK(:)            = ZWORK(:) / PTSTEP
    ZDELHEATG(:)        = ZDELHEATG(:)     + ZWORK(:)  
    ZDELHEATG_SFC(:)    = ZDELHEATG_SFC(:) + ZWORK(:)  
    PGRNDFLUX(:)        = ZGRNDFLUXN(:)
    ZFLSN_COR(:)        = 0.0
  !
  END WHERE
  !
!
!
  DO JWRK=1,INLVLS
    DO JJ=1,SIZE(PSNOWSWE(:,:),1)
      PSNOWSWE(JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWSWE(JJ,JWRK)
!VV    PEK%TSNOW%HEAT (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%HEAT (JJ,JWRK)
      PSNOWRHO  (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWRHO  (JJ,JWRK)  + &
                                    ZSNOWABLAT_DELTA(JJ)*XRHOSMIN_ES  
      PSNOWAGE(JJ,JWRK)    = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWAGE (JJ,JWRK)
      PSNOWTEMP(JJ,JWRK)    = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWTEMP(JJ,JWRK) + &
                                    ZSNOWABLAT_DELTA(JJ)*XTT  
      PSNOWLIQ (JJ,JWRK)    = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWLIQ(JJ,JWRK)        
!VV      DMK%XSNOWDZ  (JJ,JWRK)    = (1.0-ZSNOWABLAT_DELTA(JJ))*DMK%XSNOWDZ (JJ,JWRK)
    ENDDO
  ENDDO

!  
    DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PSNOWDIAMOPT(:,:),1)
        PSNOWDIAMOPT(JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWDIAMOPT(JJ,JWRK) 
        PSNOWSPHERI(JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWSPHERI(JJ,JWRK)
        PSNOWHIST (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWHIST (JJ,JWRK)
      ENDDO
    ENDDO

   DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PSNOWSWE,1)
         IF(PSNOWSWE (JJ,JWRK)==0.) THEN
             PSNOWRHO (JJ,JWRK)  = XRHOSMIN_ES  
             PSNOWTEMP(JJ,JWRK)  = XTT  
             PSNOWLIQ (JJ,JWRK)  = 0.    
             PSNOWAGE (JJ,JWRK)  = 0.
             IF(JWRK==1.) THEN
               PSNOWALB(JJ) = 0.1
             ENDIF
         ENDIF
      ENDDO
    ENDDO

!VV  ENDIF
!
!  ===============================================================
!
!  Compute snow mass budget 
!
  ZSNOW_MASS_BUDGET(:) = (ZSNOWSWE_1D(:)-ZSNOWSWE_OUT(:))/PTSTEP + PSR     (:)+ZRRSNOW (:) &
                                                                 - PEVAP   (:)-ZTHRUFAL(:) &
                                                                 + ZEVAPCOR(:)+ZSOILCOR(:)
!
!
!  ===============================================================
!
!  To Conserve mass in ISBA, the latent heat flux part of 
!  the EVAPCOR term must be weignted by the snow fraction 
!
  ZEVAPCOR (:) = ZEVAPCOR(:)*ZPSN(:) + ZSOILCOR(:)
!
! ===============================================================
!
! check suspicious low temperature
!
  DO JWRK=1,INLVLS
    !
    DO JJ=1,SIZE(PSNOWSWE,1)
      !
      IF (PSNOWSWE(JJ,JWRK)>0.0) THEN
        !
        IF (PSNOWTEMP(JJ,JWRK)<ZCHECK_TEMP) THEN
          WRITE(*,*) 'Suspicious low temperature :',PSNOWTEMP(JJ,JWRK)
          WRITE(*,*) 'At point and location      :',JJ,'LAT=',ZLAT(JJ),'LON=',ZLON(JJ)
          WRITE(*,*) 'At snow level / total layer:',JWRK,'/',INLVLS
          WRITE(*,*) 'SNOW MASS BUDGET (kg/m2/s) :',ZSNOW_MASS_BUDGET(JJ)
          WRITE(*,*) 'SWE BY LAYER      (kg/m2)  :',PSNOWSWE (JJ,1:INLVLS)
          WRITE(*,*) 'DEPTH BY LAYER      (m)    :',ZSNOWDZ  (JJ,1:INLVLS)
          WRITE(*,*) 'DENSITY BY LAYER   (kg/m3) :',PSNOWRHO(JJ,1:INLVLS)
          WRITE(*,*) 'TEMPERATURE BY LAYER (K)   :',PSNOWTEMP(JJ,1:INLVLS)
!VV          CALL ABOR1_SFX('SNOW3L_ISBA: Suspicious low temperature')  
            CALL ABORT               
        ENDIF
        !
      ELSE
        !
        !Prognostic variables forced to XUNDEF for correct outputs
        ZSNOWDZ(JJ,JWRK)=XUNDEF
        ! Careful : to compute average surface temperature in ISBA_SNOW_AGR
        ! PSNOWTEMP(JJ,1) is required when PPSN(JJ)>0 even if PSNOWSWE(JJ,1)==0
        ! (vanishing snowpack)
       ! IF (.NOT.(PEK%XPSN(JJ)>0.0.AND.JWRK==1))

 !VV   Modification to avoid crash in SVS    
 ! 
 !VV       PSNOWTEMP(JJ,JWRK) = XUNDEF
 !VV       PSNOWLIQ  (JJ,JWRK) = XUNDEF
 !VV       ZSNOWHEAT(JJ,JWRK) = XUNDEF
 !VV       PSNOWRHO (JJ,JWRK) = XUNDEF
 !VV       PSNOWAGE (JJ,JWRK) = XUNDEF
 !VV       PSNOWDIAMOPT(JJ,JWRK) = XUNDEF
 !VV       PSNOWSPHERI(JJ,JWRK) = XUNDEF
 !VV       PSNOWHIST (JJ,JWRK) = XUNDEF
      ENDIF               
    ENDDO
  ENDDO
!

! ===============================================================
!
!
!
 CONTAINS
!
!================================================================
SUBROUTINE CALL_MODEL(KSIZE1,KSIZE2,KSIZE3,KSIZE4,KMASK)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KSIZE1
INTEGER, INTENT(IN) :: KSIZE2
INTEGER, INTENT(IN) :: KSIZE3
INTEGER, INTENT(IN) :: KSIZE4
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSWE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDZ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRHO
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHEAT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTEMP
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWLIQ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDIAMOPT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSPHERI
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHIST
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWAGE
REAL, DIMENSION(KSIZE1,KSIZE2,NIMPUR) :: ZP_SNOWIMPUR 
REAL, DIMENSION(KSIZE1,NIMPUR) :: ZP_IMPWET
REAL, DIMENSION(KSIZE1,NIMPUR) :: ZP_IMPDRY 
REAL, DIMENSION(KSIZE1,KSIZE4) :: ZP_BLOWSNW
REAL, DIMENSION(KSIZE1)        :: ZP_VFRIC_T
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWALB
REAL, DIMENSION(KSIZE1)        :: ZP_SWNETSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SWNETSNOWS
REAL, DIMENSION(KSIZE1)        :: ZP_LWNETSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_PS
REAL, DIMENSION(KSIZE1)        :: ZP_SRSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_RRSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_UNLOAD
REAL, DIMENSION(KSIZE1)        :: ZP_PSN3L
REAL, DIMENSION(KSIZE1)        :: ZP_TAR
REAL, DIMENSION(KSIZE1)        :: ZP_TAC
REAL, DIMENSION(KSIZE1)        :: ZP_CT
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_TG
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_D_G
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_DZG
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_SOILHCAPZ
REAL, DIMENSION(KSIZE1)        :: ZP_SOILD
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATG
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATG_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_DELPHASEN
REAL, DIMENSION(KSIZE1)        :: ZP_DELPHASEN_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_SW_RAD
REAL, DIMENSION(KSIZE1)        :: ZP_QA
REAL, DIMENSION(KSIZE1)        :: ZP_LVTT
REAL, DIMENSION(KSIZE1)        :: ZP_LSTT
REAL, DIMENSION(KSIZE1)        :: ZP_VMOD
REAL, DIMENSION(KSIZE1)        :: ZP_VMOD_DRIFT
REAL, DIMENSION(KSIZE1)        :: ZP_LW_RAD
REAL, DIMENSION(KSIZE1)        :: ZP_RHOA
REAL, DIMENSION(KSIZE1)        :: ZP_UREF
REAL, DIMENSION(KSIZE1)        :: ZP_PRESA_SN
REAL, DIMENSION(KSIZE1)        :: ZP_EXNS
REAL, DIMENSION(KSIZE1)        :: ZP_EXNA
REAL, DIMENSION(KSIZE1)        :: ZP_DIRCOSZW
REAL, DIMENSION(KSIZE1)        :: ZP_SLOPEDIR
REAL, DIMENSION(KSIZE1)        :: ZP_AZIM
REAL, DIMENSION(KSIZE1)        :: ZP_ZREF
REAL, DIMENSION(KSIZE1)        :: ZP_Z0NAT
REAL, DIMENSION(KSIZE1)        :: ZP_Z0HNAT
REAL, DIMENSION(KSIZE1)        :: ZP_Z0EFF
REAL, DIMENSION(KSIZE1)        :: ZP_ALB
REAL, DIMENSION(KSIZE1)        :: ZP_SOILCOND
REAL, DIMENSION(KSIZE1)        :: ZP_THRUFAL
REAL, DIMENSION(KSIZE1)        :: ZP_MELTSTOT
REAL, DIMENSION(KSIZE1)        :: ZP_SNREFREEZ
REAL, DIMENSION(KSIZE1)        :: ZP_GRNDFLUX
REAL, DIMENSION(KSIZE1)        :: ZP_FLSN_COR
REAL, DIMENSION(KSIZE1)        :: ZP_RESTOREN
REAL, DIMENSION(KSIZE1)        :: ZP_EVAPCOR
REAL, DIMENSION(KSIZE1)        :: ZP_SOILCOR
REAL, DIMENSION(KSIZE1)        :: ZP_GFLXCOR
REAL, DIMENSION(KSIZE1)        :: ZP_RNSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_HSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_GFLUXSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATN
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATN_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWSFCH
REAL, DIMENSION(KSIZE1)        :: ZP_HPSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_LES3L
REAL, DIMENSION(KSIZE1)        :: ZP_LEL3L
REAL, DIMENSION(KSIZE1)        :: ZP_EVAP
REAL, DIMENSION(KSIZE1)        :: ZP_SNDRIFT
REAL, DIMENSION(KSIZE1)        :: ZP_RI
REAL, DIMENSION(KSIZE1)        :: ZP_QS
REAL, DIMENSION(KSIZE1)        :: ZP_EMISNOW
REAL, DIMENSION(KSIZE1)        :: ZP_CDSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_USTARSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_CHSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWHMASS
REAL, DIMENSION(KSIZE1)        :: ZP_VEGTYPE
REAL, DIMENSION(KSIZE1)        :: ZP_FOREST
REAL, DIMENSION(KSIZE1)        :: ZP_HVEGAPOL
REAL, DIMENSION(KSIZE1)        :: ZP_AGINGCOEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEW_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEW_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PET_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PET_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEQ_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEQ_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_ZENITH
REAL, DIMENSION(KSIZE1)        :: ZP_ANGL_ILLUM    ! Angle entre le soleil et la normal au sol et le soleil (=zenith sans pente au sol) utilisé dans TARTES
REAL, DIMENSION(KSIZE1)        :: ZP_LAT,ZP_LON
REAL, DIMENSION(KSIZE1)        :: ZP_PSN_INV
REAL, DIMENSION(KSIZE1)        :: ZP_PSN
REAL, DIMENSION(KSIZE1)        :: ZP_PSN_GFLXCOR
REAL, DIMENSION(KSIZE1)        :: ZP_WORK
!
! Declare output diagnostics
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDEND
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSPHER
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSIZE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSSA
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTYPEMEPRA
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRAM
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSHEAR
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_ACC_RAT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_NAT_RAT
!
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_12H
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_1DY
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_3DY
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_5DY
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_7DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_1DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_3DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_5DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_7DY
REAL, DIMENSION(KSIZE1) :: ZP_SNRAM_SONDE
REAL, DIMENSION(KSIZE1) :: ZP_SN_WETTHCKN
REAL, DIMENSION(KSIZE1) :: ZP_SN_REFRZNTHCKN
REAL, DIMENSION(KSIZE1) :: ZP_DEP_HIG
REAL, DIMENSION(KSIZE1) :: ZP_DEP_MOD
REAL, DIMENSION(KSIZE1) :: ZP_DEP_SUP
REAL, DIMENSION(KSIZE1) :: ZP_DEP_TOT
REAL, DIMENSION(KSIZE1) :: ZP_DEP_HUM
REAL, DIMENSION(KSIZE1) :: ZP_ACC_LEV
REAL, DIMENSION(KSIZE1) :: ZP_NAT_LEV
REAL, DIMENSION(KSIZE1) :: ZP_PRO_SUP_TYP
REAL, DIMENSION(KSIZE1) :: ZP_PRO_INF_TYP
REAL, DIMENSION(KSIZE1) :: ZP_AVA_TYP
!
REAL, DIMENSION(KSIZE1,KSIZE2,NIMPUR) :: ZP_SNOWIMP_CONC

REAL, DIMENSION(KSIZE1) :: ZP_SNOWMAK
REAL, DIMENSION(KSIZE1,1) :: ZP_DIR_SW !VV: Use of dimension of 1 for the number of spectral bands
REAL, DIMENSION(KSIZE1,1) :: ZP_SCA_SW !VV: Need to be adjusted in the future version of SVS
REAL, DIMENSION(KSIZE1,1) :: ZP_SPEC_ALB !VV
REAL, DIMENSION(KSIZE1,1) :: ZP_DIFF_RATIO !VV
REAL, DIMENSION(KSIZE1,1) :: ZP_SPEC_TOT !VV

!
REAL, PARAMETER :: ZDEPTHABS = 0.60 ! m
!
INTEGER :: JWRK, JJ, JI
!
LOGICAL :: GCOMPUTECRODIAG ! flag to compute Crocus-MEPRA diagnostics
!
!
!
! Initialize:
!
GCOMPUTECRODIAG = .FALSE.
ZP_PSN_GFLXCOR(:)  = 0.
ZP_WORK(:)         = 0.
ZP_SOILD(:)        = 0.
ZP_UNLOAD(:)       = 0.
!
!
! pack the variables
!
DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWSWE (JJ,JWRK) = PSNOWSWE(JI,JWRK)
      ZP_SNOWRHO (JJ,JWRK) = PSNOWRHO(JI,JWRK)
      ZP_SNOWHEAT(JJ,JWRK) = ZSNOWHEAT(JI,JWRK)
      ZP_SNOWAGE (JJ,JWRK) = PSNOWAGE (JI,JWRK)      
      ZP_SNOWTEMP(JJ,JWRK) = PSNOWTEMP(JI,JWRK)
      ZP_SNOWLIQ (JJ,JWRK) = PSNOWLIQ (JI,JWRK)
      ZP_SNOWDZ  (JJ,JWRK) = ZSNOWDZ  (JI,JWRK)
   ENDDO
ENDDO
!
IF (HSNOWSCHEME=='CRO') THEN

DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWDIAMOPT(JJ,JWRK) = PSNOWDIAMOPT (JI,JWRK)
      ZP_SNOWSPHERI(JJ,JWRK) = PSNOWSPHERI (JI,JWRK)
      ZP_SNOWHIST (JJ,JWRK) = PSNOWHIST  (JI,JWRK)
    ENDDO
ENDDO

    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_VFRIC_T    (JJ) = 0. !VV Value not used since threshold friction velocity is computed in Crocus. 
    ENDDO

     DO JIMP=1,NIMPUR 
      DO JWRK=1,KSIZE2
         DO JJ=1,KSIZE1
            JI = KMASK(JJ)
!VV            ZP_SNOWIMPUR(JJ,JWRK,JIMP) =PEK%TSNOW%IMPUR(JI,JWRK,JIMP)  
            ZP_SNOWIMPUR(JJ,JWRK,JIMP) = 0. !VV Set impurity content to 0 (No treatment of impurities in SVS so far)
         ENDDO
     ENDDO
    ENDDO

   DO JIMP=1,NIMPUR
    DO JJ=1,KSIZE1
       JI = KMASK(JJ)
!VV       ZP_IMPWET(JJ,JIMP)=PIMPWET(JI,JIMP)
!VV       ZP_IMPDRY(JJ,JIMP)=PIMPDRY(JI,JIMP)
       ZP_IMPWET(JJ,JIMP)=0. ! No incoming fluxes of impurities so far in SVS
       ZP_IMPDRY(JJ,JIMP)=0. ! No incoming fluxes of impurities so far in SVS
     ENDDO
   ENDDO 

   DO JWRK=1,KSIZE4
     DO JJ=1,KSIZE1
       JI = KMASK(JJ)
!VV       ZP_BLOWSNW(JJ,JWRK) = ZBLOWSNW(JI,JWRK)
          ZP_BLOWSNW(JJ,JWRK) = 0.! No incoming blowing snow fluxes in SVS
     ENDDO
   ENDDO

ELSE

   
    DO JIMP=1,NIMPUR
      DO JJ=1,KSIZE1
        ZP_IMPWET(JJ,JIMP)=XUNDEF
        ZP_IMPDRY(JJ,JIMP)=XUNDEF
        DO JWRK=1,KSIZE2
            ZP_SNOWIMPUR(JJ,JWRK,JIMP) = XUNDEF       
         ENDDO
      ENDDO
    ENDDO

   DO JWRK=1,KSIZE4
      DO JJ=1,KSIZE1
         ZP_BLOWSNW(JJ,JWRK) = XUNDEF
      ENDDO
   ENDDO
ENDIF
!  
DO JWRK=1,KSIZE3
   DO JJ=1,KSIZE1
      JI                    = KMASK           (JJ)
      ZP_TG       (JJ,JWRK) = PTG        (JI,JWRK)
      ZP_D_G      (JJ,JWRK) = PD_G       (JI,JWRK)
      ZP_SOILHCAPZ(JJ,JWRK) = PSOILHCAPZ (JI,JWRK)
   ENDDO
ENDDO
!
IF (LMEB) THEN
  DO JWRK=1,KSIZE3
    DO JJ=1,KSIZE1
      JI                    = KMASK           (JJ)
      ZP_DZG      (JJ,JWRK) = PDZG       (JI,JWRK)
    ENDDO
  ENDDO
ENDIF
!
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
   ZP_LVTT    (JJ) = ZLVTT (JI)
   ZP_LSTT    (JJ) = ZLSTT (JI)   
   ZP_EMISNOW (JJ) = ZEMISNOW(JI)   
   ZP_SNOWALB (JJ) = PSNOWALB (JI)   
   ZP_PSN3L   (JJ) = PPSN      (JI)   
   ZP_Z0NAT   (JJ) = PZ0NAT   (JI)
   ZP_Z0HNAT  (JJ) = PZ0HNAT  (JI)
   ZP_Z0EFF   (JJ) = PZ0EFF(JI)   
   ZP_RNSNOW  (JJ) = PRNSNOW (JI)
   ZP_HSNOW   (JJ) = PHSNOW  (JI)   
   ZP_HPSNOW  (JJ) = PHPSNOW (JI)   

   ZP_PS      (JJ) = PPS      (JI)
   ZP_SRSNOW  (JJ) = PSR      (JI)
   ZP_UNLOAD  (JJ) = PUNLOAD  (JI)
   ZP_CT      (JJ) = PCT      (JI)
   ZP_TAR     (JJ) = PTA      (JI)
   ZP_TAC     (JJ) = PTA      (JI)  ! Use same value for Canopy air space 
   ZP_DELHEATG(JJ) = ZDELHEATG(JI)
   ZP_DELHEATG_SFC(JJ) = ZDELHEATG_SFC(JI)
   ZP_SW_RAD  (JJ) = PSW_RAD  (JI)
   ZP_QA      (JJ) = PQA      (JI)
   ZP_VMOD    (JJ) = PVMOD    (JI)
   ZP_VMOD_DRIFT(JJ) = PWIND_DRIFT (JI)
   ZP_LW_RAD  (JJ) = PLW_RAD  (JI)
   ZP_RHOA    (JJ) = PRHOA    (JI)
   ZP_UREF    (JJ) = PUREF    (JI)
   ZP_PRESA_SN   (JJ) = PRESA_SN   (JI)
   ZP_EXNS    (JJ) = ZEXNS    (JI)
   ZP_EXNA    (JJ) = ZEXNA    (JI)
   ZP_DIRCOSZW(JJ) = ZDIRCOSZW(JI)
   ZP_SLOPEDIR(JJ) = ZSLOPEDIR(JI)
   ZP_AZIM(JJ)     = ZAZIM(JI)
   ZP_ZREF    (JJ) = PZREF    (JI)
   ZP_ALB     (JJ) = PALB     (JI)

   ZP_RRSNOW  (JJ) = ZRRSNOW  (JI)   
   ZP_SOILCOND(JJ) = ZSOILCOND(JI)

   ZP_SNDRIFT(JJ)  = PSUBLDRIFT(JI)
 
   !  
   ZP_PEW_A_COEF(JJ) = ZPEW_A_COEF(JI)
   ZP_PEW_B_COEF(JJ) = ZPEW_B_COEF(JI)
   ZP_PET_A_COEF(JJ) = ZPET_A_COEF(JI)
   ZP_PEQ_A_COEF(JJ) = ZPEQ_A_COEF(JI)      
   ZP_PET_B_COEF(JJ) = ZPET_B_COEF(JI)
   ZP_PEQ_B_COEF(JJ) = ZPEQ_B_COEF(JI)
   !
   ZP_LAT  (JJ)      = ZLAT(JI)
   ZP_LON  (JJ)      = ZLON(JI)

   ZP_ZENITH(JJ)     = PZENITH  (JI)

! VV   ZP_ANGL_ILLUM(JJ)     = PANGL_ILLUM  (JI)
   ZP_ANGL_ILLUM(JJ) = PZENITH  (JI) ! Since SVS assume a flat terrain
!
   ZP_GRNDFLUX    (JJ) = PGRNDFLUX    (JI)
   ZP_DELHEATN    (JJ) = ZDELHEATN    (JI)
   ZP_DELHEATN_SFC(JJ) = ZDELHEATN_SFC(JI)
   ZP_SNOWSFCH    (JJ) = ZSNOWSFCH    (JI)
   ZP_RESTOREN    (JJ) = ZRESTOREN    (JI)
   ZP_LES3L       (JJ) = PLES3L       (JI) 
   ZP_LEL3L       (JJ) = PLEL3L       (JI)  
   ZP_EVAP        (JJ) = PEVAP        (JI)
   !
   ZP_SWNETSNOW   (JJ) = PSWNETSNOW (JI) 
   ZP_SWNETSNOWS  (JJ) = ZSWNET_NS  (JI) 
   ZP_LWNETSNOW   (JJ) = PLWNETSNOW (JI) 
!
   ZP_SNOWMAK (JJ) = ZSNOWMAK  (JI)

   ZP_AGINGCOEF (JJ) = PAGINGCOEF  (JI)
!  
ENDDO
!
!VV DO JWRK=1,SIZE(P_DIR_SW,2)
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
!VV    ZP_DIR_SW(JJ,JWRK)=P_DIR_SW(JI,JWRK)
!VV    ZP_SCA_SW(JJ,JWRK)=P_SCA_SW(JI,JWRK)
     ZP_DIR_SW(JJ,:)= 0. !VV Need to be modified for future use of T17 snow radiation scheme in Crocus
     ZP_SCA_SW(JJ,:)= 0. !VV Need to be modified for future use of T17 snow radiation scheme in Crocus
  ENDDO
!VV!
!VVENDDO
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
!VV   ZP_VEGTYPE (JJ) = PVEGTYPE (JI,NVT_SNOW)

   ZP_VEGTYPE (JJ) = 0. ! VEGTYPE here represents the fraction of permanent snow (glacier/ice cap) on a grid point
                        ! This value is set to zero for SVS since at the moment the snowpack schemes in SVS 
                        ! only represent snow on land and do not represent snow on glaciers. 

!VV   ZP_FOREST  (JJ) = PVEGTYPE(JI,NVT_TEBD) + PVEGTYPE(JI,NVT_TRBE) + PVEGTYPE(JI,NVT_BONE)   &
!                   + PVEGTYPE(JI,NVT_TRBD) + PVEGTYPE(JI,NVT_TEBE) + PVEGTYPE(JI,NVT_TENE)   & 
!                   + PVEGTYPE(JI,NVT_BOBD) + PVEGTYPE(JI,NVT_BOND) + PVEGTYPE(JI,NVT_SHRB)   
     ZP_FOREST  (JJ) =  PFOREST(JI)

     ZP_HVEGAPOL (JJ) = PHVEGAPOL(JI)
ENDDO
!
!
! ===============================================================
! conversion of snow heat from J/m3 into J/m2
!WHERE(ZP_SNOWSWE(:,:)>0.) &
!  ZP_SNOWHEAT(:,:) = ZP_SNOWHEAT(:,:) / ZP_SNOWRHO (:,:) * ZP_SNOWSWE (:,:)  
! ===============================================================
!
ZP_PSN_INV(:)       = 0.0
ZP_PSN(:)           = ZP_PSN3L(:)
!
IF(LMEB)THEN
!
!   MEB (case of imposed surface fluxes)
!   - Prepare inputs for explicit snow scheme(s):
!     If using MEB, these are INPUTs ONLY:
!     divide fluxes by snow fraction to make "snow-relative"
!
   ZP_PSN(:)         = MAX(1.E-4, ZP_PSN3L(:))
   ZP_PSN_INV(:)     = 1.0/ZP_PSN(:)
!
   ZP_RNSNOW(:)      = ZP_RNSNOW(:)      *ZP_PSN_INV(:)
   ZP_SWNETSNOW(:)   = ZP_SWNETSNOW(:)   *ZP_PSN_INV(:)
   ZP_SWNETSNOWS(:)  = ZP_SWNETSNOWS(:)  *ZP_PSN_INV(:)
   ZP_LWNETSNOW(:)   = ZP_LWNETSNOW(:)   *ZP_PSN_INV(:)
   ZP_HSNOW(:)       = ZP_HSNOW(:)       *ZP_PSN_INV(:)
   ZP_GFLUXSNOW(:)   = ZP_GFLUXSNOW(:)   *ZP_PSN_INV(:) 
   ZP_RESTOREN(:)    = ZP_RESTOREN(:)    *ZP_PSN_INV(:) 
   ZP_SNOWHMASS(:)   = ZP_SNOWHMASS(:)   *ZP_PSN_INV(:)   
   ZP_LES3L(:)       = ZP_LES3L(:)       *ZP_PSN_INV(:)
   ZP_LEL3L(:)       = ZP_LEL3L(:)       *ZP_PSN_INV(:)
   ZP_GRNDFLUX(:)    = ZP_GRNDFLUX(:)    *ZP_PSN_INV(:)
   ZP_EVAP(:)        = ZP_EVAP(:)        *ZP_PSN_INV(:)
   ZP_HPSNOW(:)      = ZP_HPSNOW(:)      *ZP_PSN_INV(:)
   ZP_DELHEATN(:)    = ZP_DELHEATN(:)    *ZP_PSN_INV(:)
   ZP_DELHEATN_SFC(:)= ZP_DELHEATN_SFC(:)*ZP_PSN_INV(:)
   ZP_DELPHASEN(:)   = ZP_DELPHASEN(:)   *ZP_PSN_INV(:)
   ZP_DELPHASEN_SFC(:)= ZP_DELPHASEN_SFC(:)*ZP_PSN_INV(:)   
   ZP_SNOWSFCH(:)    = ZP_SNOWSFCH(:)    *ZP_PSN_INV(:)

   ZP_SRSNOW(:)      = ZP_SRSNOW(:)      *ZP_PSN_INV(:)
   ZP_RRSNOW(:)      = ZP_RRSNOW(:)      *ZP_PSN_INV(:)

   DO JJ=1,KSIZE2
      DO JI=1,KSIZE1
         ZP_SNOWSWE(JI,JJ)  = ZP_SNOWSWE(JI,JJ) *ZP_PSN_INV(JI)
         ZP_SNOWHEAT(JI,JJ) = ZP_SNOWHEAT(JI,JJ)*ZP_PSN_INV(JI)
         ZP_SNOWDZ(JI,JJ)   = ZP_SNOWDZ(JI,JJ)  *ZP_PSN_INV(JI)
      ENDDO
   ENDDO
   !
ENDIF
!
! Call ISBA-SNOW3L model:  
!  
IF (HSNOWSCHEME=='CRO') THEN 
      CALL SNOWCRO(HSNOWRES, TPTIME, LMEB, CIMPLICIT_WIND,    &
                ZP_PEW_A_COEF, ZP_PEW_B_COEF, ZP_PET_A_COEF, ZP_PEQ_A_COEF,   &
                ZP_PET_B_COEF, ZP_PEQ_B_COEF, ZP_SNOWSWE, ZP_SNOWRHO,         &
                ZP_SNOWHEAT, ZP_SNOWALB, ZP_SNOWDIAMOPT, ZP_SNOWSPHERI,          &
                ZP_SNOWHIST, ZP_SNOWAGE,ZP_SNOWIMPUR, PTSTEP, ZP_PS,          &
                ZP_SRSNOW,ZP_UNLOAD,ZP_RRSNOW, ZP_PSN3L, ZP_PRESA_SN,         &
                ZP_TAR, ZP_TG(:,1),ZP_SW_RAD,                                 &
                ZP_QA,ZP_VMOD,ZP_VMOD_DRIFT, ZP_LW_RAD, ZP_RHOA, ZP_UREF,     &
                ZP_EXNS, ZP_EXNA,                                             &
                ZP_DIRCOSZW,ZP_SLOPEDIR, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF, ZP_Z0HNAT,          &
                ZP_ALB, ZP_SOILCOND, ZP_D_G(:,1), ZP_SNOWLIQ, ZP_SNOWTEMP,    &
                ZP_SNOWDZ, ZP_THRUFAL, ZP_GRNDFLUX, ZP_EVAPCOR,ZP_GFLXCOR,    &
                ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW,ZP_RNSNOW,ZP_HSNOW, & 
                ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L, ZP_LEL3L,ZP_EVAP,          &
                ZP_SNDRIFT, ZP_RI,ZP_EMISNOW, ZP_CDSNOW,ZP_USTARSNOW,         &
                ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_VEGTYPE, ZP_ZENITH,        &
                ZP_AZIM, ZP_LAT, ZP_LON, ZP_BLOWSNW,                    &
                HSNOWDRIFT_CRO, CSNOWFPAPPUS, LSNOWDRIFT_SUBLIM,LSNOW_ABS_ZENITH,    &
                HSNOWMETAMO,HSNOWRAD,LATMORAD,ZP_DIR_SW,                      &
                ZP_SCA_SW,ZP_SPEC_ALB, ZP_DIFF_RATIO,ZP_SPEC_TOT,ZP_RESTOREN, &
                ZP_IMPWET,ZP_IMPDRY,                                          &
                HSNOWFALL, HSNOWCOND, HSNOWHOLD, HSNOWCOMP,                   &
                CSNOWZREF,ZP_SNOWMAK, LSNOWCOMPACT_BOOL,                      &
                LSNOWMAK_BOOL,LSNOWTILLER,LSELF_PROD,                         &
                LSNOWMAK_PROP,ZP_VFRIC_T, ZP_HVEGAPOL, LSNOWAGING_VAR,        &
                ZP_AGINGCOEF, LFOREST )

      ! Purely diagnostic, this routine can be called only at output time steps
!      CALL SNOWCRO_DIAG(HSNOWHOLD, HSNOWMETAMO, ZP_SNOWDZ, ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWDIAMOPT, ZP_SNOWSPHERI, ZP_SNOWAGE,  &
!              ZP_SNOWHIST, ZP_SNOWTEMP, ZP_SNOWLIQ, ZP_DIRCOSZW, ZP_SNOWIMPUR, ZP_SNOWDEND, &
!              ZP_SNOWSPHER,ZP_SNOWSIZE, ZP_SNOWSSA, ZP_SNOWTYPEMEPRA, ZP_SNOWRAM, ZP_SNOWSHEAR, &
!              ZP_ACC_RAT, ZP_NAT_RAT,   &
!              ZP_SNDPT_12H, ZP_SNDPT_1DY, ZP_SNDPT_3DY, ZP_SNDPT_5DY, ZP_SNDPT_7DY, &
!              ZP_SNSWE_1DY, ZP_SNSWE_3DY, ZP_SNSWE_5DY, ZP_SNSWE_7DY, &
!              ZP_SNRAM_SONDE, ZP_SN_WETTHCKN, ZP_SN_REFRZNTHCKN,ZP_SNOWIMP_CONC,  &
!              ZP_DEP_HIG, ZP_DEP_MOD, ZP_DEP_SUP, ZP_DEP_TOT, ZP_DEP_HUM, &
!              ZP_ACC_LEV, ZP_NAT_LEV, ZP_PRO_SUP_TYP, ZP_PRO_INF_TYP, ZP_AVA_TYP)        

!VV!
  ZP_GFLXCOR (:) = 0.0
  ZP_FLSN_COR(:) = 0.0
  ZP_SOILCOR (:) = 0.0
  ZP_MELTSTOT (:) = 0.0
  ZP_SNREFREEZ(:) = 0.0  

! Compute net shortwave and longwave fluxes
   
  DO JI=1,KSIZE1
      ZP_SWNETSNOW(JI)  = ZP_SW_RAD(JI) *( 1 - ZP_SNOWALB(JI) )
      ZP_LWNETSNOW(JI)  = ZP_RNSNOW(JI) -  ZP_SWNETSNOW(JI)
  ENDDO
!

 ELSE 
!
  CALL SNOW3L(HSNOWRES, LMEB, CIMPLICIT_WIND,                   &
             ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                 &
             ZP_PET_A_COEF, ZP_PEQ_A_COEF,ZP_PET_B_COEF, ZP_PEQ_B_COEF,    &
             ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWHEAT, ZP_SNOWALB,              &
             ZP_SNOWAGE, PTSTEP,ZP_PS, ZP_SRSNOW, ZP_UNLOAD,               &
             ZP_RRSNOW, ZP_PSN3L, ZP_PRESA_SN, ZP_TAR, ZP_TAC, ZP_TG(:,1), &
             ZP_SW_RAD, ZP_QA, ZP_VMOD, ZP_VMOD_DRIFT, ZP_LW_RAD, ZP_RHOA, ZP_UREF,       &
             ZP_EXNS, ZP_EXNA, ZP_DIRCOSZW, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,   &
             ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G(:,1),                  &
             ZP_LVTT, ZP_LSTT, ZP_SNOWLIQ,                                 &
             ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL,ZP_MELTSTOT, ZP_SNREFREEZ, &
             ZP_GRNDFLUX,ZP_EVAPCOR, ZP_SOILCOR, ZP_GFLXCOR, ZP_SNOWSFCH,  &
             ZP_DELHEATN, ZP_DELHEATN_SFC, ZP_DELPHASEN, ZP_DELPHASEN_SFC, &
             ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW, ZP_RESTOREN,       &
             ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
             ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI,                         &
             ZP_EMISNOW, ZP_CDSNOW, ZP_USTARSNOW,                          &
             ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_VEGTYPE,  ZP_FOREST,       &
              ZP_ZENITH, HSNOWDRIFT_ES, LSNOWDRIFT_SUBLIM )
 ENDIF
!



  IF(LMEB)THEN
!
! - reverse transform: back to surface-relative
!
     ZP_RNSNOW(:)      = ZP_RNSNOW(:)      /ZP_PSN_INV(:)
     ZP_SWNETSNOW(:)   = ZP_SWNETSNOW(:)   /ZP_PSN_INV(:)
     ZP_SWNETSNOWS(:)  = ZP_SWNETSNOWS(:)  /ZP_PSN_INV(:)
     ZP_LWNETSNOW(:)   = ZP_LWNETSNOW(:)   /ZP_PSN_INV(:)
     ZP_HSNOW(:)       = ZP_HSNOW(:)       /ZP_PSN_INV(:)
     ZP_LES3L(:)       = ZP_LES3L(:)       /ZP_PSN_INV(:)
     ZP_LEL3L(:)       = ZP_LEL3L(:)       /ZP_PSN_INV(:)
     ZP_GRNDFLUX(:)    = ZP_GRNDFLUX(:)    /ZP_PSN_INV(:)
     ZP_EVAP(:)        = ZP_EVAP(:)        /ZP_PSN_INV(:)
     ZP_HPSNOW(:)      = ZP_HPSNOW(:)      /ZP_PSN_INV(:)
     ZP_GFLUXSNOW(:)   = ZP_GFLUXSNOW(:)   /ZP_PSN_INV(:) 
     ZP_DELHEATN(:)    = ZP_DELHEATN(:)    /ZP_PSN_INV(:) 
     ZP_DELHEATN_SFC(:)= ZP_DELHEATN_SFC(:)/ZP_PSN_INV(:) 
     ZP_DELPHASEN(:)    = ZP_DELPHASEN(:)    /ZP_PSN_INV(:)
     ZP_DELPHASEN_SFC(:)= ZP_DELPHASEN_SFC(:)/ZP_PSN_INV(:)     
     ZP_SNOWSFCH(:)    = ZP_SNOWSFCH(:)    /ZP_PSN_INV(:) 
     ZP_RESTOREN(:)    = ZP_RESTOREN(:)    /ZP_PSN_INV(:) 

     ZP_SRSNOW(:)      = ZP_SRSNOW(:)      /ZP_PSN_INV(:)
     ZP_RRSNOW(:)      = ZP_RRSNOW(:)      /ZP_PSN_INV(:)
     DO JJ=1,KSIZE2
        DO JI=1,KSIZE1
           ZP_SNOWSWE(JI,JJ)  = ZP_SNOWSWE(JI,JJ) /ZP_PSN_INV(JI)
           ZP_SNOWHEAT(JI,JJ) = ZP_SNOWHEAT(JI,JJ)/ZP_PSN_INV(JI)
           ZP_SNOWDZ(JI,JJ)   = ZP_SNOWDZ(JI,JJ)  /ZP_PSN_INV(JI)
        ENDDO
     ENDDO
     
     ZP_SNOWHMASS(:)  = ZP_SNOWHMASS(:)/ZP_PSN_INV(:)
     ZP_THRUFAL(:)    = ZP_THRUFAL(:)  /ZP_PSN_INV(:)
     ZP_MELTSTOT(:)   = ZP_MELTSTOT(:) /ZP_PSN_INV(:)
     ZP_SNREFREEZ(:)  = ZP_SNREFREEZ(:)/ZP_PSN_INV(:)     
!
!    Final Adjustments:
!    ------------------
!    Add cooling/heating flux correction to underlying soil.
!    This term is usually active for vanishingly thin snowpacks..
!    it is put outside of the snow scheme owing to it's dependence on
!    snow fraction. It is related to a possible correction to the ground-snow
!    heat flux when it is imposed (using MEB).
!    Also, it is added as a heat sink/source here since
!    fluxes have already be computed and should not be adjusted at this point:
!    applying it to the soil has the same impact as soil freeze-thaw, in the
!    sense it is computed after the fluxes have been updated.
!    (and update heat storage diagnostic in a consistent manner)
!
!    Energy is thickness weighted, thus thicker layers receive more energy and energy
!    is evenly distributed to depth ZDEPTHABS. An
!    alternate method is to weight near surface layers more and diminish weights
!    (thus eenrgy received by each layer) with depth. Both methods conserve energy as 
!    long as vertical weights are normalized.

!    i) Determine soil depth for energy absorption:

     ZP_SOILD(:) = ZP_DZG(:,1)
     DO JJ=2,KSIZE3
        DO JI=1,KSIZE1
           IF(ZP_DZG(JI,JJ) <= ZDEPTHABS)THEN
              ZP_SOILD(JI) = ZP_DZG(JI,JJ)
           ENDIF
        ENDDO
     ENDDO

!    ii) Distribute (possible) energy to absorb vertically over some layer (defined above):

     ZP_PSN_GFLXCOR(:)  = ZP_PSN(:)*ZP_GFLXCOR(:)                                ! (W/m2)
     ZP_WORK(:)         = ZP_PSN_GFLXCOR(:)*PTSTEP/ZP_SOILD(:)

     ZP_TG(:,1)         = ZP_TG(:,1)         + ZP_WORK(:)*ZP_CT(:)*ZP_D_G(:,1)   ! (K)
     DO JJ=2,KSIZE3
        DO JI=1,KSIZE1
           IF (ZP_SOILD(JI) <= ZDEPTHABS) THEN
              ZP_TG(JI,JJ) = ZP_TG(JI,JJ)    + ZP_WORK(JI)/ZP_SOILHCAPZ(JI,JJ)   ! K
           ENDIF
        ENDDO
     ENDDO

     ZP_DELHEATG(:)     = ZP_DELHEATG(:)     + ZP_PSN_GFLXCOR(:)                 ! (W/m2)
     ZP_DELHEATG_SFC(:) = ZP_DELHEATG_SFC(:) + ZP_PSN_GFLXCOR(:)                 ! (W/m2)
!
     ZP_FLSN_COR(:)     = 0.0
!
  ELSE
!
!    To conserve energy in ISBA, the correction flux must be distributed at least
!    over the first 60cm depth. This method prevent numerical oscillations
!    especially when explicit snow vanishes. Final Adjustments are done in ISBA_CEB
!
     ZP_FLSN_COR(:) = ZP_GFLXCOR(:) ! (W/m2)
!
  ENDIF
!
!
!===============================================================
!conversion of snow heat from J/m2 into J/m3
!WHERE(ZP_SNOWSWE (:,:)>0.)
!      ZP_SNOWHEAT(:,:)=ZP_SNOWHEAT(:,:)*ZP_SNOWRHO(:,:)/ZP_SNOWSWE(:,:)  
!ENDWHERE
!===============================================================
!
! === Packing:
!
! unpack variables
!
DO JWRK=1,KSIZE2
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    PSNOWSWE(JI,JWRK) = ZP_SNOWSWE  (JJ,JWRK)
    PSNOWRHO  (JI,JWRK) = ZP_SNOWRHO  (JJ,JWRK)
    ZSNOWHEAT (JI,JWRK) = ZP_SNOWHEAT (JJ,JWRK)
    PSNOWAGE  (JI,JWRK) = ZP_SNOWAGE  (JJ,JWRK)
    PSNOWTEMP(JI,JWRK)   = ZP_SNOWTEMP (JJ,JWRK)
    PSNOWLIQ (JI,JWRK)   = ZP_SNOWLIQ  (JJ,JWRK)
    ZSNOWDZ  (JI,JWRK)   = ZP_SNOWDZ   (JJ,JWRK)
  ENDDO
ENDDO
!
!
DO JWRK=1,KSIZE2
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    PSNOWDIAMOPT(JI,JWRK) = ZP_SNOWDIAMOPT(JJ,JWRK)
    PSNOWSPHERI(JI,JWRK) = ZP_SNOWSPHERI(JJ,JWRK)
    PSNOWHIST (JI,JWRK) = ZP_SNOWHIST (JJ,JWRK)
    !M.A.: WATCH-OUT !!
    ! Commenting the PSNOWTYPE line to be able to print stats
    ! for the whole physics bus in SPS
    ! Can be un-commented if call to SNOWCRO_DIAG is un-commented
    ! Otherwise  ZP_SNOWTYPEMEPRA not initialized
    !PSNOWTYPE (JI,JWRK) = ZP_SNOWTYPEMEPRA (JJ,JWRK)
  ENDDO
ENDDO

!
DO JWRK=1,KSIZE3
   DO JJ=1,KSIZE1
      JI              = KMASK          (JJ)
      PTG    (JI,JWRK)= ZP_TG        (JJ,JWRK)
   ENDDO
ENDDO
!
DO JJ=1,KSIZE1
  JI                  = KMASK          (JJ)
  PSNOWALB (JI) = ZP_SNOWALB     (JJ)
  ZEMISNOW(JI) = ZP_EMISNOW     (JJ)
  ZCDSNOW   (JI) = ZP_CDSNOW      (JJ)
  ZUSTARSNOW(JI) = ZP_USTARSNOW   (JJ)
  ZCHSNOW   (JI) = ZP_CHSNOW      (JJ)
  ZSNOWHMASS(JI) = ZP_SNOWHMASS   (JJ)
  PRNSNOW   (JI) = ZP_RNSNOW      (JJ)
  PHSNOW    (JI) = ZP_HSNOW       (JJ)
  PHPSNOW   (JI) = ZP_HPSNOW      (JJ)
  PGFLUXSNOW(JI) = ZP_GFLUXSNOW   (JJ)
  !
  ZDELHEATG    (JI)   = ZP_DELHEATG    (JJ)
  ZDELHEATG_SFC(JI)   = ZP_DELHEATG_SFC(JJ)
  PTHRUFAL     (JI)   = ZP_THRUFAL     (JJ)
  ZEVAPCOR     (JI)   = ZP_EVAPCOR     (JJ)
  ZRI          (JI)   = ZP_RI          (JJ)
  ZQS          (JI)   = ZP_QS          (JJ)
  PGRNDFLUX    (JI)   = ZP_GRNDFLUX    (JJ)
  ZFLSN_COR    (JI)   = ZP_FLSN_COR    (JJ)
  ZDELHEATN    (JI)   = ZP_DELHEATN    (JJ)
  ZDELHEATN_SFC(JI)   = ZP_DELHEATN_SFC(JJ)
  ZDELPHASEN    (JI)  = ZP_DELPHASEN    (JJ)
  ZDELPHASEN_SFC(JI)  = ZP_DELPHASEN_SFC(JJ)  
  ZRESTOREN    (JI)   = ZP_RESTOREN    (JJ)
  ZSNOWSFCH    (JI)   = ZP_SNOWSFCH    (JJ)
  PLES3L       (JI)   = ZP_LES3L       (JJ)
  PLEL3L       (JI)   = ZP_LEL3L       (JJ)
  PEVAP        (JI)   = ZP_EVAP        (JJ)
  ZSOILCOR      (JI) = ZP_SOILCOR     (JJ)
  ZMELTSTOT    (JI)   = ZP_MELTSTOT    (JJ)
  ZSNREFREEZ   (JI)   = ZP_SNREFREEZ   (JJ)  
  PRESA_SN     (JI)   = ZP_PRESA_SN   (JJ) 

  !
  PSWNETSNOW   (JI) = ZP_SWNETSNOW   (JJ)
  ZSWNET_NS    (JI) = ZP_SWNETSNOWS  (JJ)
  PLWNETSNOW   (JI) = ZP_LWNETSNOW   (JJ)
  PSUBLDRIFT     (JI) = -1. * ZP_SNDRIFT(JJ) ! Make sure it has a negative sign (mass loss for the snow cover)
ENDDO
!
!
END SUBROUTINE CALL_MODEL
!
END SUBROUTINE SNOW_SVS2
end module snow_svs2_mod
