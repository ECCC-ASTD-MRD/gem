!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########

SUBROUTINE SNOWCRO_DIAG(HSNOWHOLD, HSNOWMETAMO, PSNOWDZ, PSNOWSWE, PSNOWRHO, PSNOWDIAMOPT, PSNOWSPHERI, PSNOWAGE, &
                        PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PDIRCOSZW, PSNOWIMP, PSNOWDEND, PSNOWSPHER, &
                        PSNOWSIZE, PSNOWSSA, PSNOWTYPEMEPRA, PSNOWRAM, PSNOWSHEAR, &
                        PACC_RAT, PNAT_RAT, &
                        PSNOWDEPTH_12H, PSNOWDEPTH_1DAYS, PSNOWDEPTH_3DAYS, PSNOWDEPTH_5DAYS, PSNOWDEPTH_7DAYS,&
                        PSNOWSWE_1DAYS, PSNOWSWE_3DAYS, PSNOWSWE_5DAYS,PSNOWSWE_7DAYS,&
                        PSNOWRAM_SONDE, PSNOW_WETTHICKNESS, PSNOW_REFTHICKNESS,PSNOWIMP_CONC,&
                        PDEP_HIG, PDEP_MOD, PDEP_SUP, PDEP_TOT, PDEP_HUM,&
                        PACC_LEV, PNAT_LEV, PPRO_SUP_TYP, PPRO_INF_TYP, PAVA_TYP)
! Diagnostics of Crocus snowpack model
! Author: M. Lafaysse, Meteo-France, October 2015
! Authors: P. Hagenmuller, Meteo-France, July 2016
! Modified Summer 2017 (P. Hagenmuller)
! Modified 10/02/20 (L.Viallon-Galinier): if PPRO_SUP_TYP is NAN -> set VLO risk
!                                         (it was sometimes the case sometimes not)
! Modified 14/04/20 (L.Viallon-Galinier): Compute maximum liquid water retention in snow layers
!                                         coherently compared to snowcro routines, accounting
!                                         to physical water retention option (SNOWHOLD)
! Modified 04/2021 (Matthieu Baron) use of function GETGRAINSIZE to compute old grain size
!
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
! Note that the Mepra diagnosis is the exact copy of the original Mepra (version in snowtools)
! and that this version explicitely contains incoherences (see comments in code and list below).
! In consequence, Mepra results should be considered for what they are worth.
!
!
!########################Mepra overall organization################################################!
!   0) Initialization of working variables
!   1) Loop on points and layers, to compute layer properties:
!       a) grain morphology: size (PSNOWSIZE), dendricity (PSNOWDEND), sphericity (PSNOWSPHER) and
!          snow type (PSNOWTYPEMEPRA)
!       b) mechanical properties: ram strength and shear strength
!
!TO DO
!1)Reduce memory print by choising appropriate types (e.g. INTEGER*1 for classes instead of REAL)
!2)Where are the projections on slope done?
!3)The input PSNWOLIQ is not coorect (too much projection)
!4)Strange definition of sphericity classes
!5)Strange definition of grain size classes
!6)Update of risk nat depending on time step
!7)Type of avalanche is weird
!8)+0.05 pour nat_rat dans le cas pro sup new
!9)add relevant threshold for ZSW (not strictly 0)

!########################Import of parameters#####################################################!
!
USE MODD_SURF_PAR,      ONLY : XUNDEF, XSURF_EPSILON
USE MODD_CSTS,          ONLY : XRHOLI, XRHOLW
USE MODD_SNOW_METAMO, ONLY : XUEPSI, XVDIAM1, XVDIAM6
USE MODD_PREP_SNOW ,ONLY : NIMPUR
USE MODE_SNOW3L
!
USE MODD_SNOW_PAR,ONLY :&
XX,XD1,XD2,XD3,&
XPERCENTAGEPORE, &
!
JPTAB_DEND,JPTAB_NODEND,&
JP_PP_PP,JP_PP_DF,JP_DF_DF,JP_DF_RG,JP_DF_FC,JP_RG_RG,JP_RG_MF,&
JP_RG_FC,JP_FC_FC,JP_FC_DH,JP_DH_DH,JP_MF_MF,JP_MF_DH,JP_MF_FC,&
!
JPACC_HIG,JPACC_MOD,JPACC_LOW,JPACC_NUL,JPACC_NAN,&
XACC_RAT_HIG,XACC_RAT_MOD,XACC_SLA_STR,&
!
JPNAT_VLO,JPNAT_LOW,JPNAT_MOA,JPNAT_MOD,JPNAT_HIG,JPNAT_VHI,JPNAT_NAN,&
JPNAT_TAB,JPNAT_ACT,&
XNAT_RAT_HIG,XNAT_RAT_MOD,&
XNAT_HEI_HIG,XNAT_HEI_MOD,XNAT_HEI_LOW, XNAT_HEI_MIN,&
!
JPPRO_SUP_NAN,JPPRO_SUP_NEW,JPPRO_SUP_WET,JPPRO_SUP_FRO,&
XPRO_SUP_CRU,XPRO_SUP_DEP,&
!
JPPRO_INF_NAN,JPPRO_INF_HAR,JPPRO_INF_SOF,&
XPRO_INF_RAM,XPRO_INF_COE,&
!
JPAVA_NEW_DRY,JPAVA_NEW_WET,JPAVA_NEW_MIX,JPAVA_SLA_SUR,&
JPAVA_MEL_SUR,JPAVA_MEL_GRO,JPAVA_NAN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
!########################Declaration of input arguments############################################!
!
IMPLICIT NONE
! Options
CHARACTER(3)        , INTENT(IN)    :: HSNOWHOLD           ! liquid water retention option
CHARACTER(3)        , INTENT(IN)    :: HSNOWMETAMO         ! metamorphism option
! Prognostic of snowcro (layer variables)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ             ! slope-perpendicular thickness (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSWE            ! mass (kg/m2)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO            ! density (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDIAMOPT        ! grain morphology variable 1 (optical diameter)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSPHERI         ! grain morphology variable 2 (sphericity)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE            ! age since snowfall (day)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWHIST           ! historical parameter (-) in {0-5}
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP           ! temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWLIQ            ! vertical liquid water content (m)
!
! Characteristics of slope
REAL, DIMENSION(:),   INTENT(IN)    :: PDIRCOSZW           ! cosine of slope angle (-)
!
! Diagnostic variables of Mepra and snowpro
! Layer variables
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDEND           ! dendricity (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSPHER          ! sphericity (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSIZE           ! grain size (m)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSSA            ! specific surface area (m2/kg)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWTYPEMEPRA      ! snow type (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWRAM            ! ram penetration strength (kgf = 9.81 N)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSHEAR          ! shear strength (kgf/dm2 = 0.981 kPa)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PACC_RAT            ! accidental ratio shear strength/stress
REAL, DIMENSION(:,:), INTENT(OUT)   :: PNAT_RAT            ! natural ratio shear strength/stress
! Point variables
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_12H    ! slope-perpendicular thickness of snow with age <= 12h  (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_1DAYS    ! slope-perpendicular thickness of snow with age <= 1 day  (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_3DAYS    ! slope-perpendicular thickness of snow with age <= 3 days (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_5DAYS    ! slope-perpendicular thickness of snow with age <= 5 days (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWDEPTH_7DAYS    ! slope-perpendicular thickness of snow with age <= 7 days (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_1DAYS      ! swe with age <= 1 day  (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_3DAYS      ! swe with age <= 3 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_5DAYS      ! swe with age <= 5 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSWE_7DAYS      ! swe with age <= 7 days (kg m-2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWRAM_SONDE      ! ramsonde top penetration slope-perpendicular thickness (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOW_WETTHICKNESS  ! top continous wet snow slope-perpendicular thickness (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOW_REFTHICKNESS  ! top continous refrozen snow slope-perpendicular thickness (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDEP_HIG            ! vertical depth of high instability (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDEP_MOD            ! vertical depth of moderate instability (m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PDEP_SUP            ! vertical depth of superior profile bottom (m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PDEP_TOT            ! total vertical snow depth (m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PDEP_HUM            ! vertical thickness of the uppest continuous block of humid snow in the sup profile
REAL, DIMENSION(:),   INTENT(OUT)   :: PACC_LEV            ! accidental risk level (0-4)
REAL, DIMENSION(:),   INTENT(INOUT) :: PNAT_LEV            ! natural risk level (0-6)
REAL, DIMENSION(:),   INTENT(INOUT) :: PPRO_SUP_TYP        ! type of superior profile (0, 4, 5, 6)
REAL, DIMENSION(:),   INTENT(OUT)   :: PPRO_INF_TYP        ! type of inferior profile (0, 1, 6)
REAL, DIMENSION(:),   INTENT(INOUT) :: PAVA_TYP            ! type of avalanche (0-6)
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PSNOWIMP_CONC
REAL, DIMENSION(:,:,:), INTENT(IN) :: PSNOWIMP
!
!########################Declaration of local variables############################################!
! scalar (no initialization required)
REAL    :: ZDIAM,&                                         ! optical diameter for SSA diagnostic (m)
           ZRAM_FIN, ZRAM_DEN, ZRAM_ANG,&                  ! for ram strength calculations
           ZSHE_SPH, ZSHE_DEN, ZSHE_MTS,&                  ! for shear strength calculations
           ZSHE_FE , ZSHE_FRE, ZSHE_FRE_SEC,&              ! for shear strength calculations
           ZEC,&                                           ! ratio between the liquid water content (kg/m3) and the maximum retention capacity
           Z_TGM,&                                         ! kind of simplified optical radius
           ZSCW,&                                          ! liquid water content (kg/m3)
           ZLWCMAX, &                                      ! maximum liquid water retention (kg/m3)
           ZSNOW_DEPTH_TOP,&                               ! vertical depth of the top of the layer (m)
           ZSNOW_HEIGHT_TOP,&                              ! vertical height (from ground) of the top of the layer (m)
           ZNAT_LEV_TMP                                    ! natural risk level before time actualization
INTEGER :: JJ,JIMP, &                                      ! loop control over point
           JST,&                                           ! loop control over layer
           ICLASS_DEND,&                                   ! dendricity class
           ICLASS_SPHER,&                                  ! sphericity class
           ICLASS_SIZE,&                                   ! grain size class
           IPRO_CLASS,&                                    ! pro sup classes for natural risk calculations
           INAT_LEV_HIG,&                                  ! naturel risk level derived only from the high instability
           INAT_LEV_MOD,&                                  ! naturel risk level derived only from the moderate instability
           INAT_LEV_LOW,&                                  ! naturel risk level derived only from no instability
           IACC_FROMNAT                                    ! accidental risk level derived from natural one
LOGICAL :: GMF,&                                           ! to be MF-like (MF/RG, MF, MF/DH, MF/FC)
           GPP,&                                           ! to be PP-like(PP, PP/DF, DF, DF/RG, DF/FC) not strictly equivalent of being dendritic
           GTHERMSTATE,GDENDRITIC                                     ! thermal_state <= 2
!
!1d (point dimension). 1D and not scalar because of loop first on JST and after on JJ.
REAL   , DIMENSION(SIZE(PSNOWSWE,1)) :: ZSNOW_DEPTH,&      ! vertical depth of current layer (m)
                                        ZWEIGH_STRESS,&    ! shear stress due to snow weight (kgf/dm2)
                                        ZSKIER_STRESS,&    ! shear stress due to a skier (kgf/dm2)
                                        ZBETA,&            ! bridging factor for skier stress calculations
                                        ZCRUST_THICKNESS,& ! vertical thickness of some crust used for weak layer detection (m)
                                        ZPRO_CRUST,&       ! vertical thickness of some crust used for sup porfile classification(m)
                                        ZDEP_INF,&         ! vertical height of the inferior porfile (m), (no init required)
                                        ZHUMTHICK,&        ! some weird height (m) used for avalanche type calculations. Non-sense.
                                        ZTHICK_BOT,&       ! vertical thickness of the bottom layer used for avalanche type calculations. Non-sense
                                        ZACC_HIG_DEP,&     ! vertical depth of high accidental instability (m)
                                        ZACC_MOD_DEP,&     ! vertical depth of mod. accidental instability (m)
                                        ZNAT_HIG_DEP,&     ! vertical depth of high natural instability (m)
                                        ZNAT_MOD_DEP,&     ! vertical depth of mod. natural instability (m)
                                        ZDEP_SUP_PRE,&     ! previous vertical depth of superior profile bottom (m)
                                        ZDEP_TOT_PRE,&     ! previous total vertical snow depth (m)
                                        ZDEP_HUM_PRE,&     ! previous vertical thickness of the uppest continuous block of humid snow in the sup profile (m)
                                        ZAVA_TYP_PRE,&     ! previous type of avalanche
                                        ZNAT_LEV_PRE,&     ! previous natural risk level
                                        ZPRO_SUP_TYP_PRE   ! previous type of superior profile
INTEGER, DIMENSION(SIZE(PSNOWSWE,1)) :: IPRO_SUP_LIM       ! index of the bottom layer of the superior profile (included)
LOGICAL, DIMENSION(SIZE(PSNOWSWE,1)) :: GRAM,&             ! no ram>2 found yet
                                        GWET,&             ! used twice with different meanings. no liq=0 found yet. no thermal sate>2 found yet.
                                        GREFROZEN,&        ! no not refrozen snow found yet
                                        GBELOW_SLAB,&      ! to be below a potnetial slab
                                        GPREV_ISNOT_MF,&   ! previous layer is not MF-like
                                        GACC_MOD_TEMP,&    ! has a real accidental moderate instabilty been found (used to detect PP weak layer)
                                        GCOULD_BE_NEW,&    ! the sup. profile could be of type new
                                        GIKNOTMF,&         ! used to determine sup. profile type
                                        GFOUND,&           ! has a layer of type MF or being at a depth > 3 cm been found
                                        GHUM,&             ! does all layers in sup pro and heigth > thr have thermal state > 2?
                                        GDRY,&             ! does all layers in sup pro and height > thr have thermal state <=2?
                                        GMEL_GRO,&         ! is there one layer in sup pro between depth > thr and height >thr with thermal state <= 2
                                        GTHERMSTATE_BOT    ! thermal_state <= 2 of bottom layer
!
REAL :: ZDRYDENSITY
INTEGER :: IMAX_USE
INTEGER, DIMENSION(SIZE(PSNOWSWE,1)) :: INLVLS_USE
INTEGER, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ISNOWHIST
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO_DIAG',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!########################Initialization of output and local variables##############################!
! Type conversion to only deal with integers
ISNOWHIST(:,:) = NINT(PSNOWHIST(:,:))
!
! Two dimensional variables (intitalization not absolutely necessary)
PSNOWDEND       = XUNDEF
PSNOWSPHER      = XUNDEF
PSNOWSIZE       = XUNDEF
PSNOWSSA        = XUNDEF
PSNOWTYPEMEPRA  = XUNDEF
PSNOWRAM        = XUNDEF
PSNOWSHEAR      = XUNDEF
PSNOWIMP_CONC=XUNDEF
PACC_RAT        = XUNDEF
PNAT_RAT        = XUNDEF
!
!Saving previous time step variables
DO JJ=1,SIZE(PSNOWSWE,1)
  ZDEP_SUP_PRE    (JJ) = PDEP_SUP    (JJ)
  ZDEP_TOT_PRE    (JJ) = PDEP_TOT    (JJ)
  ZDEP_HUM_PRE    (JJ) = PDEP_HUM    (JJ)
  ZAVA_TYP_PRE    (JJ) = PAVA_TYP    (JJ)
  ZNAT_LEV_PRE    (JJ) = PNAT_LEV    (JJ)
  ZPRO_SUP_TYP_PRE(JJ) = PPRO_SUP_TYP(JJ)
ENDDO
!
! 1d outputs
PSNOWDEPTH_12H          = 0.
PSNOWDEPTH_1DAYS        = 0.
PSNOWDEPTH_3DAYS        = 0.
PSNOWDEPTH_5DAYS        = 0.
PSNOWDEPTH_7DAYS        = 0.
PSNOWSWE_1DAYS          = 0.
PSNOWSWE_3DAYS          = 0.
PSNOWSWE_5DAYS          = 0.
PSNOWSWE_7DAYS          = 0.
PSNOWRAM_SONDE          = 0.
PSNOW_WETTHICKNESS      = 0.
PSNOW_REFTHICKNESS      = 0.
PDEP_HIG                = 0.
PDEP_MOD                = 0.
PDEP_SUP                = 0.
PDEP_TOT                = 0.
PDEP_HUM                = 0.
PACC_LEV                = JPACC_LOW
PNAT_LEV                = JPNAT_NAN
PPRO_SUP_TYP            = JPPRO_SUP_NAN
PPRO_INF_TYP            = JPPRO_INF_NAN
PAVA_TYP                = JPAVA_NAN
PAVA_TYP                = JPAVA_NAN
!
! 1d local
ZSNOW_DEPTH             = 0.
ZWEIGH_STRESS           = 0.
ZSKIER_STRESS           = 0.
ZBETA                   = 0.
ZCRUST_THICKNESS        = 0.
ZPRO_CRUST              = 0.
ZHUMTHICK               = 0.
ZTHICK_BOT              = 0.
ZACC_HIG_DEP            = XUNDEF
ZACC_MOD_DEP            = XUNDEF
ZNAT_HIG_DEP            = 0
ZNAT_MOD_DEP            = 0
IPRO_SUP_LIM            = 0
GRAM                    = .TRUE.
GWET                    = .TRUE.
GREFROZEN               = .TRUE.
GBELOW_SLAB             = .FALSE.
GPREV_ISNOT_MF          = .TRUE.
GACC_MOD_TEMP           = .FALSE.
GCOULD_BE_NEW           = .FALSE.
GIKNOTMF                = .FALSE.
GFOUND                  = .FALSE.
GHUM                    = .TRUE.
GDRY                    = .TRUE.
GMEL_GRO                = .TRUE.
GTHERMSTATE_BOT         = .FALSE.
!
IMAX_USE = 0
INLVLS_USE = 0
!
DO JJ=1,SIZE(PSNOWSWE,1)
  DO JST = 1,SIZE(PSNOWSWE,2)
    IF ( PSNOWSWE(JJ,JST)>0. ) THEN
      INLVLS_USE(JJ) = JST
      IF (IMAX_USE < JST) IMAX_USE= JST 
    ELSE 
      EXIT
    ENDIF
  ENDDO  !  end loop snow layers
  
ENDDO    ! end loop grid points
!IMAX_USE = MAXVAL(INLVLS_USE)
!
DO JST=1,IMAX_USE
  DO JJ=1,SIZE(PSNOWSWE,1)
!
    IF (JST <= INLVLS_USE(JJ)) THEN
!   Do something only in case of non-empty layer
!
!     WARNING: PSNOWLIQ is not well exported before (too much projections) RIGHT PLACE TO DO THAT?
      PSNOWLIQ(JJ,JST) = PSNOWLIQ(JJ,JST) * PDIRCOSZW(JJ)
      ZSCW = XRHOLW * PSNOWLIQ(JJ,JST) / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) !LWC (kg/m3)
!
!     ########################Grain morphology#####################################################!
!     Computes dendricity, sphericity, grain size, optical dimater and snow type from variables
!     grain1 and grain2.
!
      ! Carmagnola and Morin had the terrible idea to use SNOWGRAN1 and SNOWGRAN2 with different meanings depending on the physical options.
      ! This complicates everything and is about to change.
      ! For now :
      !
      GDENDRITIC = ( PSNOWDIAMOPT(JJ,JST)<XVDIAM6*(4.-PSNOWSPHERI(JJ,JST))-XUEPSI )   
!
      IF (GDENDRITIC) THEN
!     Dendritic case
!
!  
        PSNOWSIZE(JJ,JST)  = XUNDEF
        PSNOWSPHER(JJ,JST) = PSNOWSPHERI(JJ,JST)
        PSNOWDEND(JJ,JST)  = MAX(0.,MIN(1.,((1/XVDIAM6) * PSNOWDIAMOPT(JJ,JST)-4. + PSNOWSPHER(JJ,JST)) / &
                                                 (PSNOWSPHER(JJ,JST) - 3.))) 
!
!       10 classes of dendricity 0:[0,0.1[, ..., 9:[0.9,1.0[ (value 1.0 does not exist)
        ICLASS_DEND = MIN(INT(10 * PSNOWDEND(JJ,JST)), 9)

!       10 classes of sphericity 0:[0,0.05[, 1:[0.05,0.15[, ..., 9:[0.85,1.0]. Strange
        ICLASS_SPHER = MIN(INT(10 * PSNOWSPHER(JJ,JST) + 0.05),9)
!
!       Overall 10x10 classes from 1 to 100 (included)
        PSNOWTYPEMEPRA(JJ,JST) = JPTAB_DEND(1 + ICLASS_DEND + ICLASS_SPHER * 10)
!
!
      ELSE
!     Non dendritic case
!       
        PSNOWDEND(JJ,JST)=0.
        PSNOWSPHER(JJ,JST) = PSNOWSPHERI(JJ,JST)
        
        IF ((HSNOWMETAMO == 'B21') .OR. (HSNOWMETAMO == 'F06').OR. &
           (HSNOWMETAMO == 'S-F').OR.(HSNOWMETAMO == 'T07'))  THEN
            CALL GETGRAINSIZE_B21(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),PSNOWSIZE(JJ,JST))
        ELSE
            ! take care strange results with that
            CALL GETGRAINSIZE(PSNOWDIAMOPT(JJ,JST),PSNOWSPHERI(JJ,JST),PSNOWSIZE(JJ,JST),HSNOWMETAMO)
        ENDIF
!        
!       10 classes of sphericity 0:[0,0.05[, 1:[0.05,0.15[, ..., 9:[0.85,1.0]. Strange
        ICLASS_SPHER = MIN(INT(10 * PSNOWSPHER(JJ,JST) + 0.05),9)
!
!       3 classes of grain size in mm 0:[0,0.55[, 1:[0.55,1.05[, 2:[1.05, +inf[. Strange +0.05
        IF    (INT(10000*PSNOWSIZE(JJ,JST)+0.5).LE.5) THEN
          ICLASS_SIZE = 0
        ELSEIF(INT(10000*PSNOWSIZE(JJ,JST)+0.5).LE.10) THEN
          ICLASS_SIZE = 1
        ELSE
          ICLASS_SIZE = 2
        ENDIF
!       Overall 10x3x6 classes from 1 to 180 (included)
!       Historical variable {0,1,...,5} already defines 6 classes
        PSNOWTYPEMEPRA(JJ,JST) = JPTAB_NODEND(1 + ICLASS_SPHER + ICLASS_SIZE * 10 +&
                                              ISNOWHIST(JJ,JST) * 30)
      ENDIF
!     For all snow types (dendritic and non-dendritic)
!
!     Additional condition to define MF/RG
!     WARNING: the rounding on ZSCW here and snowttols causes very slight differences an the condition
!     WARNING: could be changed to 1e-6
      IF ((ZSCW > XSURF_EPSILON) .AND. (ISNOWHIST(JJ,JST) < 2)) THEN
        PSNOWTYPEMEPRA(JJ,JST) = JP_RG_MF
      ENDIF
!
!     ########################Thermal state######################################################!
!     Table of thermal state definition
!
!     # T(deg C) tel(kg/m3)# tel < 5      # 5 <= tel < 50     # 50 <= tel#
!     ####################################################################
!     #        T < -2      # 0            # 0                 # 0        #
!     # -2  <= T < -0.2    # 1            # 1                 # 1        #
!     #-0.2 <= T           # 2            # 3                 # 4        #
!
!     In practice, only the distinction between thermal state >2 or <= 2 is considered here, i.e.
      GTHERMSTATE = (PSNOWTEMP(JJ,JST) < 272.96) .OR. (ZSCW < 5)           !Thermstate <= 2
!
!     ########################Ram strength#########################################################!
!     Computes penetration resistance force (PSNWORAM, kgf) of ramsonde as a function of snow type
!     (PSNOWTYPEMEPRA), dendricity (PSNWODEND), sphericity (PSNOWSPHER), grain diameter (PSNWOSIZE,m)
!     and density (PSNOWRHO,kg/m3)
!
!     WARNING: ram_strength expressed in kgf (i.e. weight of 1kg = 9.81 N) corresponding to a
!     standard ramsonde. Roughly, the higher the ram strength, the higher the penetration resistance.
!     But this is not an intrinsic material property, it is dependent on the measurement instrument.
!     All threshold on ram strength in Mepra are thus also expressed in kgf.
!
!     Variable ZRAM_FIN used several times
      ZRAM_FIN = 0.17 * PSNOWRHO(JJ,JST) - 31
!
      IF (PSNOWDEND(JJ,JST) > 0) THEN
!     Dendritic case
        ZRAM_DEN = MAX(1., 0.018 * PSNOWRHO(JJ,JST) - 1.363)
        ZRAM_ANG = MAX(2., ZRAM_FIN * PSNOWSPHER(JJ,JST) + (1 - PSNOWSPHER(JJ,JST)) * (0.5 * ZRAM_FIN + 0.6))
        PSNOWRAM(JJ,JST) = PSNOWDEND(JJ,JST) * ZRAM_DEN + (1 - PSNOWDEND(JJ,JST)) * ZRAM_ANG
      ELSE
!     Non-dendritic cases
!
        SELECT CASE (INT(PSNOWTYPEMEPRA(JJ,JST)))
!       RG type
        CASE (JP_RG_RG)
          IF (PSNOWRHO(JJ,JST) < 200) THEN
            PSNOWRAM(JJ,JST) = 3
          ELSE
            PSNOWRAM(JJ,JST) = ZRAM_FIN
          ENDIF
!       RG/FC type
        CASE (JP_RG_FC)
          IF (PSNOWRHO(JJ,JST) < 200) THEN
            PSNOWRAM(JJ,JST) = 2
          ELSE
            PSNOWRAM(JJ,JST) = PSNOWSPHER(JJ,JST)  * ZRAM_FIN +&
            (1- PSNOWSPHER(JJ,JST)) * (ZRAM_FIN * (0.8-1000*PSNOWSIZE(JJ,JST)) + 2000*PSNOWSIZE(JJ,JST))
          ENDIF
!       FC and FC/DH types
        CASE (JP_FC_FC,JP_FC_DH)
          IF(PSNOWSIZE(JJ,JST) > 0.0008) THEN
            PSNOWRAM(JJ,JST) = 2
          ELSEIF (PSNOWRHO(JJ,JST) < 200) THEN
            PSNOWRAM(JJ,JST) = 3        * (0.8 - 1000 * PSNOWSIZE(JJ,JST)) + 2000 * PSNOWSIZE(JJ,JST)
          ELSE
            PSNOWRAM(JJ,JST) = ZRAM_FIN * (0.8 - 1000 * PSNOWSIZE(JJ,JST)) + 2000 * PSNOWSIZE(JJ,JST)
          ENDIF
!       MF/RG, MF, MF/DH, MF/FC types
        CASE (JP_RG_MF,JP_MF_MF,JP_MF_DH,JP_MF_FC)
          IF (GTHERMSTATE) THEN
            PSNOWRAM(JJ,JST) = MAX(10., 0.103 * PSNOWRHO(JJ,JST) - 19.666)
          ELSEIF (PSNOWRHO(JJ,JST) < 250) THEN
            PSNOWRAM(JJ,JST) = 1
          ELSEIF (PSNOWRHO(JJ,JST) < 350) THEN
            PSNOWRAM(JJ,JST) = 2
          ELSE
            PSNOWRAM(JJ,JST) = 0.16 * PSNOWRHO(JJ,JST) - 54
          ENDIF
!       All other cases
        CASE DEFAULT
          PSNOWRAM(JJ,JST) = 2
        END SELECT
      ENDIF
!
!
!     ########################Shear strength#######################################################!
!     Return the shear strength (kgf/dm2= 0.981 kPa) as a function
!     of snow type, dendricity, sphericity, grain diameter (m), historical variable, density (kg/m3),
!     liquid water content (kg/m3) and thermal state.
!
!     WARNING: One of the most time-consuming calculation, here. The evaluation of the paramaterization
!     by J.P. Navarre is elusive. The parameterization is threshold sensitive. Moreover, the complexity
!     of the parameterization appears, in some cases (e.g. C_den, C_spher), to be useless.
!
      IF ((.NOT.GTHERMSTATE).AND.((PSNOWTYPEMEPRA(JJ,JST)==JP_RG_FC).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_RG).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_RG).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC).OR.&
                                  (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF))) THEN
!     Case of wet snow of type RG, DF+RG, MF+RG, MF+DH, MF+FC, MF and RG+FC
!
        IF     (PSNOWRHO(JJ,JST) < 200) THEN
          PSNOWSHEAR(JJ,JST) = 0.1
        ELSEIF (PSNOWRHO(JJ,JST) < 320) THEN
          PSNOWSHEAR(JJ,JST) = 0.02  * PSNOWRHO(JJ,JST) - 3.9
        ELSE
          PSNOWSHEAR(JJ,JST) = 0.068 * PSNOWRHO(JJ,JST) - 18.64
        ENDIF
!
      ELSE
!     Other cases
!
!       General formula for not wet snow
!       shear_strength = ZSHE_SPH * ZSHE_DEN * ZSHE_MTS * ZSHE_FE * ZSHE_FRE) * (rho^2/10^4 - 0.6) + 0.12 with
!       ZSHE_SPH coefficient of sphericity/faceting effect,
!       ZSHE_DEN coefficient of dendrites effect,
!       ZSHE_MTS coefficient of size effect
!       ZSHE_FE  coefficient of water content
!       ZSHE_FRE coefficient of melt/freeze cycles
!
!       ZEC corresponds to the ratio between the liquid water content (kg/m3) and
!       the theoretical maximum liquid water content retaind by the snow layer,
!       defined as 5% of pore volume in B92 option.
        ! Cluzet et al 2016 : different lwc options
        ! To be optimized
        ZDRYDENSITY = PSNOWRHO(JJ,JST) - ZSCW
        !
        IF ( HSNOWHOLD == 'B92' ) THEN
          ZLWCMAX = XPERCENTAGEPORE * XRHOLW * (1 - ZDRYDENSITY/XRHOLI)
!         ZLWCMAX =  SNOWCROHOLD_0D(PSNOWRHO(JJ,JST), PSNOWLIQ(JJ,JST), PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) *  &
!                   XRHOLW / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ))
        ELSE IF ( HSNOWHOLD == 'B02' ) THEN 
          ZLWCMAX = SNOW3LHOLD(PSNOWRHO(JJ,JST), PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) * &
              XRHOLW / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ))
        ELSE IF ( HSNOWHOLD == 'SPK' ) THEN
          ZLWCMAX = SNOWSPKHOLD(ZDRYDENSITY, PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) *  &
              XRHOLW / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ))
        ELSE IF ( HSNOWHOLD == 'O04' ) THEN
          ZLWCMAX = SNOWO04HOLD_0D(ZDRYDENSITY, PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ)) * &
              XRHOLW / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ))
        ENDIF
!
!       Condition below occurs for different ice fractions depending on SNOWHOLD option
        IF (ZLWCMAX < XUEPSI) THEN
            ZEC = 0.
        ELSE
            ZEC = ZSCW / ZLWCMAX
        ENDIF
!
!       ZSHE_SPH sphericity
!       Qualitatively ZSHE_SPH increases almost linearly from (spher,ZSHE_SPH)=(0,0.45),
!       i.e. faceted snow, to (spher,ZSHE_SPH)=(1.0,1.15), i.e. roundish snow.
!       The special case for histo = 3 or 5 and spher = 0.8 is to constrained the
!       sphericity below 0.8 for snow that has been angular and humid (relatively
!       resistant, even angular).
!
        IF ((PSNOWSPHER(JJ,JST) > 0.8).AND.&
           ((ISNOWHIST(JJ,JST)==3).OR.(ISNOWHIST(JJ,JST)==5))) THEN
          ZSHE_SPH = 1.05
        ELSEIF (PSNOWSPHER(JJ,JST) <= 0.25) THEN
          ZSHE_SPH = 0.45  + 0.7  * PSNOWSPHER(JJ,JST)
        ELSEIF (PSNOWSPHER(JJ,JST) <= 0.50) THEN
          ZSHE_SPH = 0.625 + 1.0 * (PSNOWSPHER(JJ,JST) - 0.25)
        ELSEIF (PSNOWSPHER(JJ,JST) <= 0.75) THEN
          ZSHE_SPH = 0.875 + 0.6 * (PSNOWSPHER(JJ,JST) - 0.50)
        ELSE
          ZSHE_SPH = 1.025 + 0.5 * (PSNOWSPHER(JJ,JST) - 0.75)
        ENDIF
!       ZSHE_DEN dendricity
!       Qualitatively, ZSHE_DEN decreases almost linearly from (dendr,ZSHE_DEN) = (0,1),
!       i.e. old snow, to (dendr,ZSHE_DEN) = (1.0,0.45), i.e. recent snow.
!
        IF (PSNOWDEND(JJ,JST) <= 0.25) THEN
          ZSHE_DEN = 1.0 - 0.4 *  PSNOWDEND(JJ,JST)
        ELSEIF (PSNOWDEND(JJ,JST) <= 0.50) THEN
          ZSHE_DEN = 0.9 - 0.4 * (PSNOWDEND(JJ,JST) - 0.25)
        ELSEIF (PSNOWDEND(JJ,JST) <= 0.75) THEN
          ZSHE_DEN = 0.8 - 0.8 * (PSNOWDEND(JJ,JST) - 0.50)
        ELSE
          ZSHE_DEN = 0.6 - 0.6 * (PSNOWDEND(JJ,JST) - 0.75)
        ENDIF
!
!       ZSHE_MTS
!       As expected, the grain size has no effect on dendritic snow.
!       Qualitatively, on non dendritic snow with a grain size larger than a kind of simplified optical radius (Z_TGM),
!       ZSHE_MTS accounts for the fact that the larger the grains, the smaller the strength, especially in case
!       of faceted snow types.
!
        Z_TGM = 0.0004 - 0.0001 * PSNOWSPHER(JJ,JST)
!
        IF ((PSNOWDEND(JJ,JST) > 0).OR.(PSNOWSIZE(JJ,JST) <= Z_TGM)) THEN
          ZSHE_MTS = 1
        ELSE
          ZSHE_MTS = 1 - (0.8 - 0.2 * PSNOWSPHER(JJ,JST)) * 530 * (PSNOWSIZE(JJ,JST) - Z_TGM)
        ENDIF
!
!       ZSHE_FE Wetting of snow
!       For low liquid water content (lower than thresold depending on density),
!       C_fe increases linearly with the liquid water content (max = 1.1). For larger values,
!       C_fe decreases linearly with the liquid water content until C_fe reaches about 0.62,
!       then it continues to decrease but more slowly.
!
        IF (ZEC <= 0.1) THEN
          ZSHE_FE = 1 + ZEC
        ELSEIF (ZEC <= 0.3) THEN
          ZSHE_FE = 1.335 - 2.35 * ZEC
        ELSEIF (ZEC <  0.9) THEN
          ZSHE_FE = 0.750 - 0.4  * ZEC
        ELSE
          ZSHE_FE = MAX(0.15,MIN(0.35,(PSNOWRHO(JJ,JST)-ZSCW)/1000))
        ENDIF
!
!       ZSHE_FRE takes into account the strengthening effect of melt/freeze cycles.
!       In case of never humid snow (histo = 0 or 1), C_fre = 1.0
!       In case of "has been humid snow" (histo = 2 or 3) (but no cycle),
!       C_fre = 1.0 if the snow is still humid, and C_fre about 2.0 if completely dry (ec=0).
!       In case of melt/freeze cycle snow, C_fre decreases from about 1.95 to 1 with
!       liquid water content, with 1/scw trend.
!
        IF ((ISNOWHIST(JJ,JST) <= 1).OR.(ZEC > 0.5)) THEN
          ZSHE_FRE = 1.
        ELSEIF (ZEC<XSURF_EPSILON) THEN !Changed from==0 to <1e-16
!         C_fre_sec = 1.5 * ((1.15 + 0.2 * (1.-spher)) / 1.15) * (1 + 0.2/C_mts)
          ZSHE_FRE = (1.7608695652 - 0.2608695652 * PSNOWSPHER(JJ,JST)) * (1. + 0.2/ZSHE_MTS)
        ELSEIF (ISNOWHIST(JJ,JST) <= 3) THEN
          ZSHE_FRE = 1
        ELSEIF (ZEC <= 0.1) THEN
          ZSHE_FRE = -2. * ZEC + 1.5
        ELSE
          ZSHE_FRE = -0.75 * ZEC + 1.375
        ENDIF
!
        PSNOWSHEAR(JJ,JST) = MAX(0.05, &
        ZSHE_SPH * ZSHE_DEN * ZSHE_MTS * ZSHE_FE * ZSHE_FRE *&
        (PSNOWRHO(JJ,JST)*PSNOWRHO(JJ,JST) / 10000. -0.6) + 0.12)
!
      ENDIF
!
!
!     ########################Strength stress ratio######################################################!
!     Computes the ratio PACC_RAT between the shear strength (PSNOWSHEAR kgf/dm2 = 0.981 kPa) of the
!     current layer and the shear stress due to a skier and to the overlying snow weight at the TOP of the layer
!     (= bottom of just above layer).
!     and the ratio PNAT_RAT between the shear strength (PSNOWSHEAR kgf/dm2 = 0.981 kPa) and the shear
!     stress due only to the weight of the layers at the BOTTOM of the current layer.
!
!     Stress skier takes into account that the additional stress induced by the skier
!     is more distributed (so lower) in the snowpack far below the skier than just below the skier.
!     The calculation is based on the simplification of the semi-infinite uniform elastic layer
!     theory developed by Boussinesq. Expert rules are used to account for
!     bridging effects of different snow types (calculation of beta).
!
!     Low value of beta indicates high bridging effect (i.e. decrease of max shear stress)
!     slab_thr = 1.5          #Threshold on shear strength (kgf/dm2)
!     beta_refrozen = 0.5     #(MF/RG, MF, MF/DH or MF/FC) and thermal_state <= 2
!     beta_humid = 1.1        #(MF/RG, MF, MF/DH or MF/FC) and thermal_state >  2
!     beta_evolved = 1.0      #Other snow types and shear_strength >  slab_thr
!     beta_recent_dry = 1.2   #Other snow types and shear_strength <= slab_thr
!
!     WARNING: two slightly different versions in "python snowtools" or "MEPRA fortran" of PACC_RAT.
!     The version in MEPRA fortran is used to calculate the accidental risk and the version
!     in snowtools just provides the value of 'MEPRA_ACCIDENTAL_RATIO' but is not used for
!     other calculations. Differences can be observed.
!
!     WARNING: not defined for slope angle = 0 and first top layer
!
!
!     Update of beta (bridging factor)
      IF (JST==1) THEN
        ZBETA(JJ) = 0
      ELSE
        IF ((PSNOWTYPEMEPRA(JJ,JST-1)==JP_MF_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST-1)==JP_RG_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST-1)==JP_MF_FC).OR.&
            (PSNOWTYPEMEPRA(JJ,JST-1)==JP_MF_DH)) THEN
          IF (PSNOWTEMP(JJ,JST-1) < 272.96) THEN
            ZBETA(JJ) = ZBETA(JJ) + 0.5 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ELSE
            ZBETA(JJ) = ZBETA(JJ) + 1.1 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ENDIF
        ELSE
          IF (PSNOWSHEAR(JJ,JST-1) > 1.5) THEN
            ZBETA(JJ) = ZBETA(JJ) + 1.0 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ELSE
            ZBETA(JJ) = ZBETA(JJ) + 1.2 * PSNOWDZ(JJ,JST-1) / PDIRCOSZW(JJ)
          ENDIF
        ENDIF
      ENDIF
!
      IF((PDIRCOSZW(JJ)<1).AND.(JST>1)) THEN
!      ACC_RAT is defined only for slope > 0 and not for the top layer
!      total_stress(jst) = weight_stress(jst-1) + beta(jst) * skier_stress(jst-1)
        PACC_RAT(JJ,JST) = ZWEIGH_STRESS(JJ) + ZBETA(JJ) / ZSNOW_DEPTH(JJ) * ZSKIER_STRESS(JJ)
        PACC_RAT(JJ,JST) = PSNOWSHEAR(JJ,JST) / PACC_RAT(JJ,JST)
      ELSE
        PACC_RAT(JJ,JST) = XUNDEF !-1 in previous versions
      ENDIF
!
!     Update of snow depth (the update is here and not before!)
      ZSNOW_DEPTH(JJ) = ZSNOW_DEPTH(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
!     Update of skier_stress
      IF (    ZSNOW_DEPTH(JJ) < 0.10) THEN
        ZSKIER_STRESS(JJ) = -15.0 * ZSNOW_DEPTH(JJ) + 4.00
      ELSEIF (ZSNOW_DEPTH(JJ) < 0.15) THEN
        ZSKIER_STRESS(JJ) = -10.0 * ZSNOW_DEPTH(JJ) + 3.50
      ELSEIF (ZSNOW_DEPTH(JJ) < 0.20) THEN
        ZSKIER_STRESS(JJ) = - 8.0 * ZSNOW_DEPTH(JJ) + 3.20
      ELSEIF (ZSNOW_DEPTH(JJ) < 0.35) THEN
        ZSKIER_STRESS(JJ) = - 4.0 * ZSNOW_DEPTH(JJ) + 2.40
      ELSEIF (ZSNOW_DEPTH(JJ) < 0.50) THEN
        ZSKIER_STRESS(JJ) = - 2.0 * ZSNOW_DEPTH(JJ) + 1.70
      ELSEIF (ZSNOW_DEPTH(JJ) < 0.80) THEN
        ZSKIER_STRESS(JJ) = - 1.5 * ZSNOW_DEPTH(JJ) + 1.45
      ELSE
        ZSKIER_STRESS(JJ) = 0
      ENDIF
!
      ZSKIER_STRESS(JJ) = 1.4 * ZSKIER_STRESS(JJ)
!
!     Update of shear stress due to overlying layers
      ZWEIGH_STRESS(JJ) = ZWEIGH_STRESS(JJ) + &
      PSNOWRHO(JJ,JST) * PSNOWDZ(JJ,JST) * SQRT(1-PDIRCOSZW(JJ)*PDIRCOSZW(JJ)) / 100.
!
      IF((PDIRCOSZW(JJ)<1).AND.(JST>1)) THEN
!     ACC_RAT is not defined for slope = 0 degrees and for the top layer
        PNAT_RAT(JJ,JST) = PSNOWSHEAR(JJ,JST) / ZWEIGH_STRESS(JJ)
      ELSE
        PNAT_RAT(JJ,JST) = XUNDEF!in previous -1
      ENDIF
!
!     ##########################################Accidental risk###################################
!     The accidental risk is later (see below) combined with the natural (spontaneous) risk.
!
!     Depth of top of current layer
      ZSNOW_DEPTH_TOP = ZSNOW_DEPTH(JJ) - PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
!     Update of total thickness of crusts above current layer (included)
!     In practice, it is not necessary to calculate it for the whole snowpack
!     when a high instability is already found
      IF((PSNOWTEMP(JJ,JST) < 272.96).AND.&
         ((PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
          (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF).OR.&
          (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
          (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC))) THEN
          ZCRUST_THICKNESS(JJ) = ZCRUST_THICKNESS(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
      ENDIF
!
!     If no slab structure was found yet (GBELOW_SLAB = to be under a slab)
      IF(.NOT.GBELOW_SLAB(JJ)) THEN
!
!       Update of GPREV_ISNOT_MF(JJ) = True if layer just below is MF-like snow
        IF((PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC)) THEN
!
          GPREV_ISNOT_MF(JJ) = .FALSE.
!
!       Update of GBELOW_SLAB(JJ) = True if current layer is below a slab
        ELSEIF((PSNOWSHEAR(JJ,JST).GT.XACC_SLA_STR).AND.&
               (GPREV_ISNOT_MF(JJ))                .AND.&
               ((PSNOWTYPEMEPRA(JJ,JST)==JP_DF_RG).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_RG).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_FC).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_DF))) THEN
!
          GBELOW_SLAB(JJ) = .TRUE.
!
        ELSE
!       Update of GPREV_ISNOT_MF(JJ) = True if layer just above is MF-like snow
          GPREV_ISNOT_MF(JJ) = .TRUE.
        ENDIF
!
!     If a slab structure was found and but no high instability -> searching for a weak layer
      ELSEIF ((ZSNOW_DEPTH_TOP>= 0.01).AND.&
              (ZSNOW_DEPTH_TOP<  1.00).AND.&
              (ZACC_HIG_DEP(JJ) == XUNDEF)) THEN
!
!       Permanent weak layer composed of facets or depth hoar
        IF((PSNOWTYPEMEPRA(JJ,JST)==JP_FC_FC).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_FC_DH).OR.&
           (PSNOWTYPEMEPRA(JJ,JST)==JP_DH_DH)) THEN
!
!         If skier_ratio < threshold1
          IF((PACC_RAT(JJ,JST).LT.XACC_RAT_HIG)) THEN
            ZACC_HIG_DEP(JJ) = ZSNOW_DEPTH(JJ) - 0.5 * PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
!           Condition on the crust thickness above weak layer
            IF(ZCRUST_THICKNESS(JJ) <= 0.01) THEN
              PACC_LEV(JJ) = JPACC_HIG
            ELSE
              PACC_LEV(JJ) = JPACC_MOD
            ENDIF
!
!         Else if no moderate was found yet and skier_ratio < XACC_RAT_MOD
          ELSEIF((ZACC_MOD_DEP(JJ) == XUNDEF).AND.(PACC_RAT(JJ,JST).LT.XACC_RAT_MOD)) THEN
            ZACC_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ) - 0.5 * PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
!           Condition on the crust thickness above weak layer
            IF(ZCRUST_THICKNESS(JJ) <= 0.01) THEN
              PACC_LEV(JJ) = JPACC_MOD
            ELSE
              PACC_LEV(JJ) = JPACC_LOW
            ENDIF
!
          ENDIF
!
!       Else if not already found looking for a
!       temporary weak layer composed of precipitation particles or decomposed snow
!       Note that there are no condition on ratio_acc
        ELSEIF((.NOT.GACC_MOD_TEMP(JJ)).AND.&
               ((PSNOWTYPEMEPRA(JJ,JST)==JP_PP_PP).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_PP_DF).OR.&
                (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_DF))) THEN
!
          ZACC_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ) - 0.5 * PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
          GACC_MOD_TEMP(JJ) = .TRUE.
!
!         Condition on the crust thickness above weak layer
          IF(ZCRUST_THICKNESS(JJ) <= 0.01) THEN
            PACC_LEV(JJ) = JPACC_MOD
          ELSE
            PACC_LEV(JJ) = JPACC_LOW
          ENDIF
!
        ENDIF
      ENDIF
!
!
!     #######################Superior profile detection and classification#########################!
!
!     Note that there are more superior profile classes in original Mepra but they are never used
!     in practice. In consequence, they are not considered here. Three types of sup. profile are
!     considered: NEW (new snow), FRO (top frozen), WET (top wet).
!
!     Profile type NEW
!     The superior profile, in case of the presence recent snow, is composed of
!     the recent snow layers (even not at surface) and all above layers and possibly directly
!     below humid layers (thermal_state>2), where thin ice crusts (connected thickness < 1 cm and
!     composed of dry MF-like snow layers) are "ignored".!
!
!     Profile types WET or FRO
!     Roughly, these profiles are composed of connected layers of MF-like snow if the top of
!     this block of MF layers is close to the snow surface. There are additional constraints to
!     define these profiles but they do not make clear sense to me.
!
!
!     MF-like snow
      GMF = (PSNOWTYPEMEPRA(JJ,JST)==JP_RG_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_MF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_DH).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_MF_FC)
!
!     Recent snow
      GPP = (PSNOWTYPEMEPRA(JJ,JST)==JP_PP_PP).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_PP_DF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_DF).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_RG).OR.&
            (PSNOWTYPEMEPRA(JJ,JST)==JP_DF_FC)
!
      IF(GPP) THEN
!       When recent snow is found, then
!       the superior profile type is, in all cases, NEW,
!       the superior profile reaches at least this layer,
!       what is direclty below potentially belongs to the superior profile and
!       the thickness of the connex crusts is re-initiliazed to zero.
        PPRO_SUP_TYP (JJ) = JPPRO_SUP_NEW
        IPRO_SUP_LIM (JJ) = JST
        PDEP_SUP     (JJ) = ZSNOW_DEPTH(JJ)
        GCOULD_BE_NEW(JJ) = .TRUE.
        ZPRO_CRUST   (JJ) = 0
!
      ELSEIF(GCOULD_BE_NEW(JJ)) THEN
!     if the current layer potentially belongs to the sup. profile
!
        IF((.NOT.GTHERMSTATE).AND.(ZPRO_CRUST(JJ) < XPRO_SUP_CRU)) THEN
!         if the current layer is humid and the above crust thickness < 1 cm
!         then the current layer and potentially what is below, belongs to the sup. profile and
!         the crust thickness is re-initialized to zero.
          IPRO_SUP_LIM(JJ) = JST
          PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)
          ZPRO_CRUST  (JJ) = 0
!
        ELSEIF(GMF.AND.GTHERMSTATE) THEN
!       elseif the current layer is a crust
!       then we increase the crust thickness accordingly
          ZPRO_CRUST(JJ) = ZPRO_CRUST(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
        ELSE
!       in all other cases, the layer and what is below do not belong to the sup. profile.
          GCOULD_BE_NEW(JJ) = .FALSE.

          IF((ZPRO_CRUST(JJ) < XPRO_SUP_CRU).AND.(JST.GT.1)) THEN
            IPRO_SUP_LIM(JJ) = JST-1
            PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ) - PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
          ENDIF
!
        ENDIF
!
      ELSEIF(PPRO_SUP_TYP(JJ).NE.JPPRO_SUP_NEW) THEN
!     if a sup. profile of type NEW was not found yet. COULD BE IMPROVED FOR CLARITY ...
!
        IF(.NOT.GFOUND(JJ)) THEN
!       if not found yet, we are looking for the first layer of type MF or being at a
!       depth > 3 cm. This layer thermal state determines whether the profile is FRO or WET
!       (not consistent...). In case of this layer being not MF, we are not sure yet that the
!       sup. profile will be FRO/WET or NAN.
!
          IF(GMF.OR.(ZSNOW_DEPTH(JJ).GT.XPRO_SUP_DEP)) THEN
            GFOUND      (JJ) = .TRUE.
            !IPRO_SUP_LIM(JJ) = JST
            !PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)
!
            IF(GTHERMSTATE) THEN
              PPRO_SUP_TYP(JJ) = JPPRO_SUP_FRO
            ELSE
              PPRO_SUP_TYP(JJ) = JPPRO_SUP_WET
            ENDIF
!
            GIKNOTMF(JJ) = .NOT.GMF
!
            IF(GIKNOTMF(JJ)) THEN
              PPRO_SUP_TYP(JJ) = JPPRO_SUP_NAN
            ELSE
              IPRO_SUP_LIM(JJ) = JST
              PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)
            ENDIF
!
          ENDIF
!
        ELSEIF(.NOT.GIKNOTMF(JJ)) THEN
          IF(GMF) THEN
!         if the current layer is of type MF
!         then we increase the sup. profile and
!         we do not care anymore on the type of the "found" layer
            IPRO_SUP_LIM(JJ) = JST
            PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)
          ELSE
            GIKNOTMF(JJ)=.TRUE.
          ENDIF
!
        ENDIF
      ENDIF
!
!
!     #######################Other cumulative quantities##########################################!
!     Compute depth and SWE of snow with age < X days (X in {1,3,5,7})
!     WARNING: there is some rounding on age that may difer between fortran and python snowtools
!     WARNNG: PSNOWSWE and PSNOWDEPTH are slope perpendicular. The projection is done later in surfex
      IF(PSNOWAGE(JJ,JST) <= 7) THEN
        PSNOWDEPTH_7DAYS(JJ) = PSNOWDEPTH_7DAYS(JJ) + PSNOWDZ (JJ,JST)
        PSNOWSWE_7DAYS  (JJ) = PSNOWSWE_7DAYS  (JJ) + PSNOWSWE(JJ,JST)
        IF(PSNOWAGE(JJ,JST) <= 5) THEN
          PSNOWDEPTH_5DAYS(JJ) = PSNOWDEPTH_5DAYS(JJ) + PSNOWDZ (JJ,JST)
          PSNOWSWE_5DAYS  (JJ) = PSNOWSWE_5DAYS  (JJ) + PSNOWSWE(JJ,JST)
!
          IF(PSNOWAGE(JJ,JST) <= 3) THEN
            PSNOWDEPTH_3DAYS(JJ) = PSNOWDEPTH_3DAYS(JJ) + PSNOWDZ (JJ,JST)
            PSNOWSWE_3DAYS  (JJ) = PSNOWSWE_3DAYS  (JJ) + PSNOWSWE(JJ,JST)
!
            IF(PSNOWAGE(JJ,JST) <= 1) THEN
              PSNOWDEPTH_1DAYS(JJ) = PSNOWDEPTH_1DAYS(JJ) + PSNOWDZ (JJ,JST)
              PSNOWSWE_1DAYS  (JJ) = PSNOWSWE_1DAYS  (JJ) + PSNOWSWE(JJ,JST)

              IF(PSNOWAGE(JJ,JST) <= 0.5) THEN
                PSNOWDEPTH_12H(JJ) = PSNOWDEPTH_12H(JJ) + PSNOWDZ (JJ,JST)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
!
!     Ramsonde top penetration
      IF ((GRAM(JJ)).AND.(PSNOWRAM(JJ,JST)<=2.)) THEN
        PSNOWRAM_SONDE(JJ)   = PSNOWRAM_SONDE(  JJ) + PSNOWDZ (JJ,JST)
      ELSE
        GRAM(JJ)=.FALSE.
      ENDIF
!
!     Depth of top wet snow
      IF ((GWET(JJ)).AND.(PSNOWLIQ(JJ,JST).GT.0)) THEN
        PSNOW_WETTHICKNESS(JJ) = PSNOW_WETTHICKNESS(JJ) + PSNOWDZ(JJ,JST)
      ELSE
        GWET(JJ)=.FALSE.
      ENDIF
!     Depth of top refrozen snow
      IF (GREFROZEN(JJ).AND.(ISNOWHIST(JJ,JST)>=2).AND.(PSNOWTEMP(JJ,JST)<273.15)) THEN
        PSNOW_REFTHICKNESS(JJ) = PSNOW_REFTHICKNESS(JJ) + PSNOWDZ(JJ,JST)
      ELSE
        GREFROZEN(JJ)=.FALSE.
      ENDIF
!
!     Specific surface area
      PSNOWSSA(JJ,JST) = 6. / (XRHOLI*PSNOWDIAMOPT(JJ,JST))
      !
      !SnowImpurity Concentration
      
      DO JIMP=1,NIMPUR !Modif the cond
        PSNOWIMP_CONC(JJ,JST,JIMP)=PSNOWIMP(JJ,JST,JIMP)/(1000*PSNOWSWE(JJ,JST)) 
      END DO 
!
    ENDIF
  END DO
END DO

!     #######################Re-initialization####################################################!
DO JJ=1,SIZE(PSNOWSWE,1)
! Some weird condition to add the bottom crust of the sup profile new if it is thin enough
  IF((ZPRO_CRUST   (JJ).LT.XPRO_SUP_CRU ).AND.&
     (PPRO_SUP_TYP (JJ).EQ.JPPRO_SUP_NEW).AND.&
     (GCOULD_BE_NEW(JJ))) THEN
    IPRO_SUP_LIM(JJ) = JST
    PDEP_SUP    (JJ) = ZSNOW_DEPTH(JJ)
  ENDIF
  PDEP_TOT   (JJ) = ZSNOW_DEPTH(JJ)
  ZSNOW_DEPTH(JJ) = 0
  ZDEP_INF   (JJ) = PDEP_TOT(JJ) - PDEP_SUP(JJ)
  GWET       (JJ) = .TRUE.
END DO
!     #######################Second layer-point loop##############################################!
DO JST=1,IMAX_USE
  DO JJ=1,SIZE(PSNOWSWE,1)
!
    IF (JST <= INLVLS_USE(JJ)) THEN
!   Do something only in case of non-empty layer
!
!     Update of snow depth (depth of bottom of current layer).
      ZSNOW_HEIGHT_TOP = PDEP_TOT(JJ) - ZSNOW_DEPTH(JJ)
      ZSNOW_DEPTH(JJ)  = ZSNOW_DEPTH(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
      IF(JST<=IPRO_SUP_LIM(JJ)) THEN
!     Inside superior profile
!
        ZSCW = XRHOLW * PSNOWLIQ(JJ,JST) / (PSNOWDZ(JJ,JST) * PDIRCOSZW(JJ))    !lwc (kg/m3)
        GTHERMSTATE = (PSNOWTEMP(JJ,JST) < 272.96).OR.(ZSCW < 5)                !thermstate <= 2
!       Storing the current GTHERMSTATE to get the one of the bottom layer of sup profile
        GTHERMSTATE_BOT(JJ) = GTHERMSTATE
!!!!!!!!!       Storing the current layer thickness to get the one of the bottom layer
!!!!!!!!!!        ZTHICK_BOT(JJ) = PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
!
!       Calculations of the height of the uppest continuous block of humid snow in the sup.
!       profile (slightly different from SNOWWETTHICKNESS)
        IF(GWET(JJ).AND.(.NOT.GTHERMSTATE)) THEN
          PDEP_HUM(JJ) = PDEP_HUM(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
        ELSEIF(PDEP_HUM(JJ).GT.0) THEN
          GWET(JJ) = .FALSE.
        ENDIF
!
!       Determination of natural risk based on the stress/strength ratio
        IF(ZSNOW_HEIGHT_TOP.GE.XNAT_HEI_MIN) THEN
!       Only searching above a certain height
!
          IF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW) THEN
!         Case of NEW sup. profile
!
            IF(PNAT_RAT(JJ,JST).LE.XNAT_RAT_HIG) THEN
              ZNAT_HIG_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ELSEIF(PNAT_RAT(JJ,JST).LE.(XNAT_RAT_MOD + 0.05)) THEN !verrue d'origine inconnue
              ZNAT_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ENDIF
!
          ELSE
!         Case of WET or FRO profiles (code not accesible for nan profile)
!
            IF(PNAT_RAT(JJ,JST).LE.XNAT_RAT_HIG) THEN
              ZNAT_HIG_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ELSEIF(PNAT_RAT(JJ,JST).LE.XNAT_RAT_MOD) THEN
              ZNAT_MOD_DEP(JJ) = ZSNOW_DEPTH(JJ)
            ENDIF
!
          ENDIF
!
!         Calculations used for avalanche type determination
!         GHUM and GDRY used only for profile sup. NEW
          IF(GTHERMSTATE) THEN
            GHUM(JJ) = .FALSE.
          ELSE
            GDRY(JJ) = .FALSE.
          ENDIF
!         GMEL_GRO used only for profiles sup. WET and NAN
!         GMEL_GRO = There is no dry layer at a depth > 10 cm and height > 10 cm in the sup profile
          IF((ZSNOW_DEPTH(JJ).GT.XNAT_HEI_MIN).AND.GTHERMSTATE) THEN
            GMEL_GRO(JJ) = .FALSE.
          ENDIF
!
!         ZHUMTHICK used for profiles sup. FRO and NAN. Non-sense ...
          ZHUMTHICK(JJ) = ZHUMTHICK(JJ) + PSNOWDZ(JJ,JST) / PDIRCOSZW(JJ)
        ENDIF
!
      ELSE
!     inside inferior profile to determine its type
        IF(PPRO_INF_TYP(JJ).NE.JPPRO_INF_HAR) THEN
          IF(ZSNOW_HEIGHT_TOP.GT.(XPRO_INF_COE * ZDEP_INF(JJ))) THEN !NOOOOOOOOOOOOOOT SURE of top
            IF(PSNOWRAM(JJ,JST).LT.XPRO_INF_RAM) THEN
              PPRO_INF_TYP(JJ) = JPPRO_INF_SOF
            ELSE
              PPRO_INF_TYP(JJ) = JPPRO_INF_HAR
            ENDIF
          ENDIF
        ENDIF
!
      ENDIF
!
!
    ENDIF
!
  ENDDO
ENDDO

!First loop only on points
DO JJ=1,SIZE(PSNOWSWE,1)
!
! #######################Natural risk level determination##########################################!
!
  SELECT CASE(INT(PPRO_SUP_TYP(JJ)))
    CASE(JPPRO_SUP_NEW)
      IPRO_CLASS = 0
    CASE(JPPRO_SUP_WET)
      IPRO_CLASS = 15
    CASE(JPPRO_SUP_FRO)
      IPRO_CLASS = 30
    CASE DEFAULT
      IPRO_CLASS = 0
  END SELECT
!
  IF    (ZNAT_HIG_DEP(JJ) > XNAT_HEI_HIG) THEN
    INAT_LEV_HIG = JPNAT_TAB(15 + IPRO_CLASS)
  ELSEIF(ZNAT_HIG_DEP(JJ) > XNAT_HEI_MOD) THEN
    INAT_LEV_HIG = JPNAT_TAB(12 + IPRO_CLASS)
  ELSEIF(ZNAT_HIG_DEP(JJ) > XNAT_HEI_LOW) THEN
    INAT_LEV_HIG = JPNAT_TAB(9  + IPRO_CLASS)
  ELSEIF(ZNAT_HIG_DEP(JJ) > 0           ) THEN
    INAT_LEV_HIG = JPNAT_TAB(6  + IPRO_CLASS)
  ELSE
    INAT_LEV_HIG = JPNAT_TAB(3  + IPRO_CLASS)
  ENDIF
!
  IF    (ZNAT_MOD_DEP(JJ) > XNAT_HEI_HIG) THEN
    INAT_LEV_MOD = JPNAT_TAB(14 + IPRO_CLASS)
  ELSEIF(ZNAT_MOD_DEP(JJ) > XNAT_HEI_MOD) THEN
    INAT_LEV_MOD = JPNAT_TAB(11 + IPRO_CLASS)
  ELSEIF(ZNAT_MOD_DEP(JJ) > XNAT_HEI_LOW) THEN
    INAT_LEV_MOD = JPNAT_TAB(8  + IPRO_CLASS)
  ELSEIF(ZNAT_MOD_DEP(JJ) > 0)            THEN
    INAT_LEV_MOD = JPNAT_TAB(5  + IPRO_CLASS)
  ELSE
    INAT_LEV_MOD = JPNAT_TAB(2  + IPRO_CLASS)
  ENDIF
!
  IF    (PDEP_SUP(JJ) > XNAT_HEI_HIG) THEN
    INAT_LEV_LOW = JPNAT_TAB(13 + IPRO_CLASS)
  ELSEIF(PDEP_SUP(JJ) > XNAT_HEI_MOD) THEN
    INAT_LEV_LOW = JPNAT_TAB(10 + IPRO_CLASS)
  ELSEIF(PDEP_SUP(JJ) > XNAT_HEI_LOW) THEN
    INAT_LEV_LOW = JPNAT_TAB(7  + IPRO_CLASS)
  ELSEIF(PDEP_SUP(JJ) > 0)            THEN
    INAT_LEV_LOW = JPNAT_TAB(4  + IPRO_CLASS)
  ELSE
    INAT_LEV_LOW = JPNAT_TAB(1  + IPRO_CLASS)
  ENDIF
!
  IF(INAT_LEV_HIG.NE.JPNAT_MOA) THEN
    PNAT_LEV(JJ) = MAX(INT(INAT_LEV_HIG),MAX(INT(INAT_LEV_MOD),INT(INAT_LEV_LOW)))
  ELSE
    PNAT_LEV(JJ) = JPNAT_MOA
  ENDIF
!
! #######################Avalanche type determination##############################################!
!
! Does not make sense, whatever...
!!!!!!!!!!!!!!!  ZHUMTHICK(JJ) = ZHUMTHICK(JJ) + ZTHICK_BOT(JJ)
!
  IF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW) THEN
! Profile sup. NEW
!
    IF(GHUM(JJ)) THEN
      PAVA_TYP(JJ) = JPAVA_NEW_WET
    ELSEIF(GDRY(JJ)) THEN
      PAVA_TYP(JJ) = JPAVA_NEW_DRY
    ELSE
      PAVA_TYP(JJ) = JPAVA_NEW_MIX
    ENDIF
!
  ELSEIF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NAN) THEN
! Profile sup. NAN
    PAVA_TYP(JJ) = JPAVA_NAN
!
  ELSEIF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_WET) THEN
! Profile sup. WET
!
    IF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_NAN) THEN
!   Profile inf NAN
      IF(GMEL_GRO(JJ)) THEN
        PAVA_TYP(JJ) = JPAVA_MEL_GRO
      ELSE
        PAVA_TYP(JJ) = JPAVA_MEL_SUR
      ENDIF
!
    ELSEIF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_HAR) THEN
!   Profile inf HAR
      PAVA_TYP(JJ) = JPAVA_MEL_SUR
    ELSE
!   Profile inf SOF
      PAVA_TYP(JJ) = JPAVA_MEL_GRO
    ENDIF
!
  ELSEIF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_FRO) THEN
! Profile sup. WET
    IF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_NAN) THEN
!   Profile inf NAN
      IF((.NOT.GTHERMSTATE_BOT(JJ)).AND.(ZHUMTHICK(JJ).GT.(PDEP_TOT(JJ)/3.0))) THEN
        PAVA_TYP(JJ) = JPAVA_MEL_GRO
      ELSE
        PAVA_TYP(JJ) = JPAVA_MEL_SUR
      ENDIF

    ELSEIF(PPRO_INF_TYP(JJ).EQ.JPPRO_INF_HAR) THEN
!   Profile inf HAR
      PAVA_TYP(JJ) = JPAVA_MEL_SUR
    ELSE
!   Profile inf SOF
      PAVA_TYP(JJ) = JPAVA_MEL_GRO
    ENDIF
!
    IF(PNAT_LEV(JJ).EQ.JPNAT_VLO) PAVA_TYP(JJ) = JPAVA_NAN
!
  ENDIF
!

! #######################Risk natural actualization################################################!
!
  ZNAT_LEV_TMP = PNAT_LEV(JJ)
! only in case, there is not a NAN risk in the previous step
  IF(ZNAT_LEV_PRE(JJ).NE.JPNAT_NAN) THEN
!
!   For type sup. NEW
    IF((  PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW   )   .AND.&
       (  PAVA_TYP    (JJ).EQ.ZAVA_TYP_PRE(JJ))   .AND.&
       (  PDEP_TOT    (JJ).LT.ZDEP_TOT_PRE(JJ))   .AND.&
       (  PDEP_SUP    (JJ).LE.ZDEP_SUP_PRE(JJ))   .AND.&
       (((PAVA_TYP(JJ).EQ.JPAVA_NEW_MIX).AND.(PDEP_HUM(JJ).LE.ZDEP_HUM_PRE(JJ))).OR.&
        (PAVA_TYP(JJ).NE.JPAVA_NEW_MIX))) THEN
!
      ZNAT_LEV_TMP = JPNAT_ACT(INT(PNAT_LEV(JJ) + ZNAT_LEV_PRE(JJ)*7 + 1))
!
!     Weird addtional instruction
      IF(( PAVA_TYP(JJ).NE.JPAVA_NEW_MIX).AND.&
         ((PNAT_LEV(JJ).EQ.JPNAT_HIG).OR.(PNAT_LEV(JJ).EQ.JPNAT_VHI)).AND.&
         ((ZNAT_LEV_PRE(JJ).EQ.JPNAT_MOD).OR.(ZNAT_LEV_PRE(JJ).EQ.JPNAT_HIG).OR.&
          (ZNAT_LEV_PRE(JJ).EQ.JPNAT_VHI))) THEN
         ZNAT_LEV_TMP = JPNAT_MOD
      ENDIF
!
    ENDIF
!
    IF((PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NEW   )   .AND.&
       (PAVA_TYP    (JJ).EQ.ZAVA_TYP_PRE(JJ))   .AND.&
       (PDEP_TOT    (JJ).LT.ZDEP_TOT_PRE(JJ))   .AND.&
       (PNAT_LEV    (JJ).EQ.JPNAT_MOA)) THEN
!
        IF(ZNAT_LEV_PRE(JJ).EQ.JPNAT_HIG) ZNAT_LEV_TMP = JPNAT_MOD
        IF(ZNAT_LEV_PRE(JJ).EQ.JPNAT_VHI) ZNAT_LEV_TMP = JPNAT_HIG
!
    ENDIF
!
    PNAT_LEV(JJ) = ZNAT_LEV_TMP
!
!   For type sup. WET or FRO
    IF(((PPRO_SUP_TYP    (JJ).EQ.JPPRO_SUP_WET).OR.(PPRO_SUP_TYP    (JJ).EQ.JPPRO_SUP_FRO)).AND.&
       ((ZPRO_SUP_TYP_PRE(JJ).EQ.JPPRO_SUP_WET).OR.(ZPRO_SUP_TYP_PRE(JJ).EQ.JPPRO_SUP_FRO)).AND.&
       (ZDEP_INF(JJ).GE.(ZDEP_TOT_PRE(JJ)-ZDEP_SUP_PRE(JJ)-0.05))) THEN
!
      IF(PNAT_LEV(JJ).EQ.JPNAT_MOD) THEN
        PNAT_LEV(JJ) = JPNAT_LOW
      ENDIF
!
      IF((PPRO_SUP_TYP    (JJ).EQ.JPPRO_SUP_WET).AND.&
         ((PNAT_LEV(JJ).EQ.JPNAT_HIG).OR.(PNAT_LEV(JJ).EQ.JPNAT_VHI))) THEN
!
        IF((ZNAT_LEV_PRE(JJ).EQ.JPNAT_MOD).OR.&
           (ZNAT_LEV_PRE(JJ).EQ.JPNAT_HIG).OR.&
           (ZNAT_LEV_PRE(JJ).EQ.JPNAT_VHI)) THEN
          PNAT_LEV(JJ) = JPNAT_MOD
        ELSEIF(ZNAT_LEV_PRE(JJ).EQ.JPNAT_LOW) THEN
          PNAT_LEV(JJ) = JPNAT_LOW
        ENDIF
      ENDIF
    ENDIF
  ENDIF
!
! For no profil sup. (JPPRO_SUP_NAN)
  IF(PPRO_SUP_TYP(JJ).EQ.JPPRO_SUP_NAN) THEN
    PNAT_LEV(JJ) = JPNAT_VLO
  ENDIF
!
!
! Case with no snow (most of the time we do not even enter snowcrodiag)
  IF(PDEP_TOT(JJ).LE.0.0000001) THEN
    PNAT_LEV(JJ) = JPNAT_NAN
  ENDIF
!
! Combination of natural and accidental risk levels
!
  IF    ((PNAT_LEV(JJ).EQ.JPNAT_VHI).OR.(PNAT_LEV(JJ).EQ.JPNAT_HIG)) THEN
    IACC_FROMNAT = JPACC_HIG
  ELSEIF((PNAT_LEV(JJ).EQ.JPNAT_MOD).OR.(PNAT_LEV(JJ).EQ.JPNAT_MOA)) THEN
    IACC_FROMNAT = JPACC_MOD
  ELSE
    IACC_FROMNAT = JPACC_LOW
  ENDIF
!
  IF(PDEP_TOT(JJ).LE.XNAT_HEI_LOW) THEN
    PACC_LEV(JJ) = JPACC_LOW
  ENDIF
!
!
! Accidental risk
  IF(ZACC_HIG_DEP(JJ) == XUNDEF) THEN
    IF(ZACC_MOD_DEP(JJ) == XUNDEF) THEN
      ZACC_HIG_DEP(JJ) = -1
      ZACC_MOD_DEP(JJ) = -1
    ELSE
      ZACC_HIG_DEP(JJ) = -1
      ZACC_MOD_DEP(JJ) = ZACC_MOD_DEP(JJ)
    ENDIF
  ELSE
    ZACC_HIG_DEP(JJ) = ZACC_HIG_DEP(JJ)
    ZACC_MOD_DEP(JJ) = -1
  ENDIF
!
  IF((IACC_FROMNAT.GE.PACC_LEV(JJ)).OR.(PDIRCOSZW(JJ).GT.0.8)) THEN
    PACC_LEV(JJ) = IACC_FROMNAT
    PDEP_HIG(JJ) = ZNAT_HIG_DEP(JJ)
    PDEP_MOD(JJ) = ZNAT_MOD_DEP(JJ)
  ELSE
    PDEP_HIG(JJ) = ZACC_HIG_DEP(JJ)
    PDEP_MOD(JJ) = ZACC_MOD_DEP(JJ)
  ENDIF
!
  IF((PDIRCOSZW(JJ).GT.0.8).AND.(PACC_LEV(JJ).EQ.JPACC_LOW)) THEN
    PACC_LEV(JJ) = JPACC_NUL
  ENDIF
!
! Verrue additionelle pour le cas sans pente
! Does not make sense, since these variables could be also defined for flat terrain
  IF(PDIRCOSZW(JJ).GT.0.99) THEN
    PDEP_SUP(JJ) = 0
    PDEP_HUM(JJ) = 0
    PDEP_HIG(JJ) = XUNDEF
    PDEP_MOD(JJ) = XUNDEF
    PAVA_TYP(JJ) = 0
    PNAT_LEV(JJ) = JPNAT_VLO
    PACC_LEV(JJ) = JPACC_NUL
  ENDIF
!
  IF(PDEP_HIG(JJ).EQ.0) THEN
    PDEP_HIG(JJ) = XUNDEF
  ENDIF
!
  IF(PDEP_MOD(JJ).EQ.0) THEN
    PDEP_MOD(JJ) = XUNDEF
  ENDIF
ENDDO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO_DIAG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCRO_DIAG
