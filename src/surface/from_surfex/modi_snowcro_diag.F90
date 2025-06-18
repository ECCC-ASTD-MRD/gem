!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!auto_modi:spll_snowcro_diag.D
MODULE MODI_SNOWCRO_DIAG
INTERFACE
SUBROUTINE SNOWCRO_DIAG(HSNOWHOLD, HSNOWMETAMO, PSNOWDZ, PSNOWSWE, PSNOWRHO, PSNOWDIAMOPT, PSNOWSPHERI, PSNOWAGE, &
                        PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PDIRCOSZW, PSNOWIMP, PSNOWDEND, PSNOWSPHER, &
                        PSNOWSIZE, PSNOWSSA, PSNOWTYPEMEPRA, PSNOWRAM, PSNOWSHEAR, &
                        PACC_RAT, PNAT_RAT, &
                        PSNOWDEPTH_12H, PSNOWDEPTH_1DAYS, PSNOWDEPTH_3DAYS, PSNOWDEPTH_5DAYS, PSNOWDEPTH_7DAYS,&
                        PSNOWSWE_1DAYS, PSNOWSWE_3DAYS, PSNOWSWE_5DAYS,PSNOWSWE_7DAYS,&
                        PSNOWRAM_SONDE, PSNOW_WETTHICKNESS, PSNOW_REFTHICKNESS,PSNOWIMP_CONC,&
                        PDEP_HIG, PDEP_MOD, PDEP_SUP, PDEP_TOT, PDEP_HUM,&
                        PACC_LEV, PNAT_LEV, PPRO_SUP_TYP, PPRO_INF_TYP, PAVA_TYP)
IMPLICIT NONE
CHARACTER(3)        , INTENT(IN)    :: HSNOWHOLD           ! liquid water retention option
CHARACTER(3)        , INTENT(IN)    :: HSNOWMETAMO         ! metamorphism option
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ             ! slope-perpendicular thickness (m) (projection in reproj_diag_isban)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSWE            ! mass (kg/m2)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO            ! density (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDIAMOPT        ! grain morphology variable 1 (optical diameter)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWSPHERI         ! grain morphology variable 2 (sphericity)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE            ! age since snowfall (day)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWHIST           ! historical parameter (-) in {0-5}
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP           ! temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWLIQ            ! vertical liquid water content (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PDIRCOSZW           ! cosine of slope angle (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDEND           ! dendricity (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSPHER          ! sphericity (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSIZE           ! grain size (m)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSSA            ! specific surface area (m2/kg)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWTYPEMEPRA      ! snow type (-)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWRAM            ! ram penetration strength (kgf = 9.81 N)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWSHEAR          ! shear strength (kgf/dm2 = 0.981 kPa)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PACC_RAT            ! accidental ratio shear strength/stress
REAL, DIMENSION(:,:), INTENT(OUT)   :: PNAT_RAT            ! natural ratio shear strength/stress
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
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PSNOWIMP_CONC
REAL, DIMENSION(:,:,:), INTENT(IN) :: PSNOWIMP
END SUBROUTINE SNOWCRO_DIAG
END INTERFACE
END MODULE MODI_SNOWCRO_DIAG
