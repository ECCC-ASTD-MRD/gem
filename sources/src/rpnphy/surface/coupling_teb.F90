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
!     ###############################################################################
SUBROUTINE COUPLING_TEB2(PTSTEP, KYEAR, KMONTH, KDAY, PTIME, PTSUN, PZENITH, PAZIM,      &
             PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PRAIN, PSNOW, PLW, PDIR_SW,     &
             PSCA_SW, PSW_BANDS, PPS, PPA, PSFTQ, PSFTH, PSFU, PSFV,                     &
             PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PLAT                                      )
!     ###############################################################################
!
use sfclayer_mod, only: sl_prelim

!****  *COUPLING_TEB$n * - Driver for TEB 
!
!    PURPOSE
!    -------
!
!**  METHOD
!    ------
!
!    REFERENCE
!    ---------
!      
!
!    AUTHOR
!    ------
!     V. Masson 
!
!    MODIFICATIONS
!    -------------
!      Original    01/2004
!---------------------------------------------------------------
!
!
USE MODD_CSTS,       ONLY : XRD, XCPD, XP00
!
USE MODD_TEB,        ONLY :                                                     &
!                      TTIME,                                                   &
                       XT_CANYON, XQ_CANYON, XTI_BLD,                           &
                       XT_ROOF, XT_ROAD, XT_WALL, XWS_ROOF, XWS_ROAD,           &
                       XWSNOW_ROOF,  XTSNOW_ROOF, XRSNOW_ROOF,                  &
                       XASNOW_ROOF, XESNOW_ROOF, XTSSNOW_ROOF,                  &
                       XWSNOW_ROAD,  XTSNOW_ROAD, XRSNOW_ROAD,                  &
                       XASNOW_ROAD, XESNOW_ROAD, XTSSNOW_ROAD,                  &
                       XH_TRAFFIC, XLE_TRAFFIC, XH_INDUSTRY, XLE_INDUSTRY,      &
                       XZ0_TOWN, XZ0_ROOF, XZ0_ROAD,                            &
                       XBLD, XBLD_HEIGHT, XWALL_O_HOR, XCAN_HW_RATIO,           &
                       XALB_ROOF, XEMIS_ROOF, XHC_ROOF,XTC_ROOF, XD_ROOF,       &
                       XALB_ROAD, XEMIS_ROAD, XHC_ROAD,XTC_ROAD, XD_ROAD,       &
                       XALB_WALL, XEMIS_WALL, XHC_WALL,XTC_WALL, XD_WALL,       &
                       XSVF_ROAD, XSVF_WALL
!
USE MODD_TOWN,       ONLY :                                                     &
                       NNI, XTOWN,                                              &
                       XU_CANYON,                                               &
                       XRN_ROOF,XH_ROOF,XLE_ROOF,XLES_ROOF,                     &
                       XGFLUX_ROOF,XRUNOFF_ROOF,                                &
                       XRN_ROAD,XH_ROAD,XLE_ROAD,XLES_ROAD,                     &
                       XGFLUX_ROAD,XRUNOFF_ROAD,                                &
                       XRN_WALL,XH_WALL,XLE_WALL,XGFLUX_WALL,                   &
                       XRNSNOW_ROOF,XHSNOW_ROOF,XLESNOW_ROOF,                   &
                       XGSNOW_ROOF,XMELT_ROOF,                                  &
                       XRNSNOW_ROAD,XHSNOW_ROAD,XLESNOW_ROAD,                   &
                       XGSNOW_ROAD,XMELT_ROAD,                                  &
                       XRN,XH,XLE,XGFLUX,XEVAP,XRUNOFF,                         &
                       XCH,XRI,XUSTAR,                                          &
                      XTRAD_IN,XTRAD_SUN,XTRAD_SHADE,XTRAD_RFSUN,XTRAD_RFSHADE, &
                      XTGLOBE_SUN,XTGLOBE_SHADE,XTGLOBE_RFSUN,XTGLOBE_RFSHADE,  &
                      XTWETB,XTWETB_ROOF,                                       &
                XUTCI_IN,XUTCI_OUTSUN,XUTCI_OUTSHADE,XUTCI_RFSUN,XUTCI_RFSHADE, &
                    XWBGT_OUTSUN, XWBGT_OUTSHADE, XWBGT_RFSUN, XWBGT_RFSHADE,   &
            XUTCIC_IN,XUTCIC_OUTSUN,XUTCIC_OUTSHADE,XUTCIC_RFSUN,XUTCIC_RFSHADE &
           ,XTRFZT,XTRDZT,XURDZU  				                                   &
           ,XQ1,XQ2,XQ3,XQ4,XQ5,XQ6,XQ7,XQ8,XQ9,XQ10,XQ11,XQ12,XQ13
!
USE MODI_TEB2
! 
implicit none
#include <arch_specific.hf>
!
!
REAL    :: XUNDEF    ! undefined value
PARAMETER( XUNDEF = 999. )  
!

!*      0.1    declarations of arguments
INTEGER,            INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,            INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,            INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,               INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,               INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(:), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(:), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(:), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(:), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(:), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(:,:),INTENT(IN) :: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN) :: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH   ! zenithal angle       (radian from the vertical)
REAL, DIMENSION(:), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(:), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                            !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(:), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PLAT      ! latitude                              (deg)
!
REAL, DIMENSION(:), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(:), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
!
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(:,:),INTENT(OUT):: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(OUT):: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
!*      0.2    declarations of local variables
!
INTEGER                     :: STAT   ! status
INTEGER                     :: JSWB   ! loop counter on shortwave spectral bands
!         
REAL, DIMENSION(SIZE(PTA))  :: ZQA    ! specific humidity                 (kg/kg)
REAL, DIMENSION(SIZE(PTA))  :: ZEXNA  ! Exner function at forcing level
REAL, DIMENSION(SIZE(PTA))  :: ZEXNS  ! Exner function at surface level
REAL, DIMENSION(SIZE(PTA))  :: ZWIND  ! wind speed
REAL, DIMENSION(SIZE(PTA))  :: ZDIR   ! wind direction
!
REAL, DIMENSION(SIZE(PTA))  :: ZDIR_ALB      ! direct albedo of town
REAL, DIMENSION(SIZE(PTA))  :: ZSCA_ALB      ! diffuse albedo of town
!
REAL, DIMENSION(SIZE(PTA))  :: ZH_TRAFFIC    ! anthropogenic sensible
!                                            ! heat fluxes due to traffic
REAL, DIMENSION(SIZE(PTA))  :: ZLE_TRAFFIC   ! anthropogenic latent
!                                            ! heat fluxes due to traffic
REAL, DIMENSION(SIZE(PTA))  :: ZRESA_TOWN    ! aerodynamical resistance
!
REAL                        :: ZBEGIN_TRAFFIC_TIME ! start traffic time (solar time, s)
REAL                        :: ZEND_TRAFFIC_TIME   ! end traffic time   (solar time, s)
REAL, DIMENSION(SIZE(PTA))  :: ZDIR_SW       ! total direct SW
REAL, DIMENSION(SIZE(PTA))  :: ZSCA_SW       ! total diffuse SW
!-------------------------------------------------------------------------------------
! Preliminaries:
!-------------------------------------------------------------------------------------
!
! specific humidity (conversion from kg/m3 to kg/kg)
!
ZQA(:) = PQA(:) 
!
! Exner functions
!
ZEXNS(:)     = (PPS(:)/XP00)**(XRD/XCPD)
ZEXNA(:)     = (PPA(:)/XP00)**(XRD/XCPD)
!
!! scalar fluxes
!!
!PSFTS(:,:) = 0.
!
! broadband radiative fluxes
!
ZDIR_SW(:) = PDIR_SW(:,1)
ZSCA_SW(:) = PSCA_SW(:,1)
!
! wind
!
stat = sl_prelim(PTA,PQA,PU,PV,PPS,PUREF,spd_air=ZWIND,dir_air=ZDIR)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Anthropogenic fluxes (except building heating)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZBEGIN_TRAFFIC_TIME = 21600.
ZEND_TRAFFIC_TIME   = 64800.
!
WHERE(       PTSUN>ZBEGIN_TRAFFIC_TIME   &
      .AND.  PTSUN<ZEND_TRAFFIC_TIME     )
  ZH_TRAFFIC  (:) = XH_TRAFFIC   (:)
  ZLE_TRAFFIC (:) = XLE_TRAFFIC  (:)
ELSEWHERE
  ZH_TRAFFIC  (:) = 0.
  ZLE_TRAFFIC (:) = 0.   
END WHERE
!
!--------------------------------------------------------------------------------------
! Over Urban surfaces/towns:
!--------------------------------------------------------------------------------------
!

CALL TEB2  (XT_CANYON, XQ_CANYON, XU_CANYON,                          &
     XTI_BLD,                                                         &
     XT_ROOF, XT_ROAD, XT_WALL, XWS_ROOF,XWS_ROAD,                    &
     XWSNOW_ROOF, XTSNOW_ROOF,                                        &
     XRSNOW_ROOF, XASNOW_ROOF,                                        &
     XTSSNOW_ROOF, XESNOW_ROOF,                                       &
     XWSNOW_ROAD, XTSNOW_ROAD,                                        &
     XRSNOW_ROAD, XASNOW_ROAD,                                        &
     XTSSNOW_ROAD, XESNOW_ROAD,                                       &
     PPS, PPA, ZEXNS, ZEXNA, PTA, ZQA, PRHOA,                         &
     PLW, ZDIR_SW, ZSCA_SW, PZENITH,                                  &
     PRAIN, PSNOW, PZREF, PUREF, ZWIND, ZDIR,                         &
     ZH_TRAFFIC, ZLE_TRAFFIC, XH_INDUSTRY, XLE_INDUSTRY,              &
     PTSTEP,                                                          &
     XZ0_TOWN, XZ0_ROOF, XZ0_ROAD,                                    &
     XBLD, XBLD_HEIGHT, XWALL_O_HOR, XCAN_HW_RATIO,                   &
     XALB_ROOF, XEMIS_ROOF,                                           &
     XHC_ROOF,XTC_ROOF,XD_ROOF,                                       &
     XALB_ROAD, XEMIS_ROAD, XSVF_ROAD,                                &
     XHC_ROAD,XTC_ROAD,XD_ROAD,                                       &
     XALB_WALL, XEMIS_WALL, XSVF_WALL, PLAT,                          &
     XHC_WALL,XTC_WALL,XD_WALL,                                       &
     XRN_ROOF, XH_ROOF, XLE_ROOF, XLES_ROOF, XGFLUX_ROOF,             &
     XRUNOFF_ROOF,                                                    &
     XRN_ROAD, XH_ROAD, XLE_ROAD, XLES_ROAD, XGFLUX_ROAD,             &
     XRUNOFF_ROAD,                                                    &
     XRN_WALL, XH_WALL, XLE_WALL, XGFLUX_WALL,                        &
     XRNSNOW_ROOF, XHSNOW_ROOF, XLESNOW_ROOF, XGSNOW_ROOF,            &
     XMELT_ROOF,                                                      &
     XRNSNOW_ROAD, XHSNOW_ROAD, XLESNOW_ROAD, XGSNOW_ROAD,            &
     XMELT_ROAD,                                                      &
     XRN, XH, XLE, XGFLUX, XEVAP, XRUNOFF,                            &
     XUSTAR, XCH, XRI,                                                &
     PTRAD, PEMIS, ZDIR_ALB, ZSCA_ALB, ZRESA_TOWN ,                   &
     XTRAD_IN, XTRAD_SUN, XTRAD_SHADE, XTRAD_RFSUN, XTRAD_RFSHADE,    &
     XTGLOBE_SUN, XTGLOBE_SHADE, XTGLOBE_RFSUN, XTGLOBE_RFSHADE,      &
     XTWETB,XTWETB_ROOF,                                              &
     XUTCI_IN, XUTCI_OUTSUN, XUTCI_OUTSHADE, XUTCI_RFSUN, XUTCI_RFSHADE, &
      XWBGT_OUTSUN, XWBGT_OUTSHADE, XWBGT_RFSUN, XWBGT_RFSHADE,          &
     XUTCIC_IN,XUTCIC_OUTSUN,XUTCIC_OUTSHADE,XUTCIC_RFSUN,XUTCIC_RFSHADE &
      ,XTRFZT,XTRDZT,XURDZU &
     ,XQ1,XQ2,XQ3,XQ4,XQ5,XQ6,XQ7,XQ8,XQ9,XQ10,XQ11,XQ12,XQ13   )
!
!
!-------------------------------------------------------------------------------------
! Outputs:
!-------------------------------------------------------------------------------------
!
! Momentum fluxes
!
PSFU = 0.
PSFV = 0.
WHERE (ZWIND(:)>0.)
  PSFU(:) = - PRHOA(:) * XUSTAR(:)**2 * PU(:) / ZWIND(:)
  PSFV(:) = - PRHOA(:) * XUSTAR(:)**2 * PV(:) / ZWIND(:)
END WHERE
!
! Heat and CO2 fluxes
!
PSFTH(:)        = XH(:)
PSFTQ(:)        = XEVAP(:)
!
DO JSWB=1,SIZE(PSW_BANDS)
  PDIR_ALB(:,JSWB) = ZDIR_ALB(:)
  PSCA_ALB(:,JSWB) = ZSCA_ALB(:)
END DO
!
END SUBROUTINE COUPLING_TEB2
