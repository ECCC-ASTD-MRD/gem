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
!   ##########################################################################
    SUBROUTINE URBAN_SNOW_EVOL(                                               &
                     PT_CANYON, PQ_CANYON, PU_CANYON,                         &
                     PTS_ROOF,PTS_ROAD,PTS_WALL,                              &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPS, PTA, PQA, PRHOA,                                    &
                     PLW_RAD,                                                 &
                     PSR, PZREF, PUREF, PVMOD,                                &
                     PTSTEP,                                                  &
                     PBLD_HEIGHT,                                             &
                     PEMIS_ROAD, PSVF_ROAD,                                   &
                     PEMIS_WALL, PSVF_WALL,                                   &
                     PDN_ROOF, PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,          &
                     PDF_ROAD, PDN_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,&
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PLW_S_TO_N                                               )
!   ##########################################################################
!
!!****  *URBAN_SNOW_EVOL*  
!!
!!    PURPOSE
!!    -------
!
!     
!!**  METHOD
!     ------
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SNOW_PAR, ONLY : XANSMIN, XANSMAX, XANS_TODRY, XRHOSMIN, XRHOSMAX, &
                          XZ0SN, XZ0HSN, XWCRN, SWE_CRIT
USE MODD_CSTS,     ONLY : XSTEFAN
!
USE MODE_SURF_SNOW_FRAC
!
USE MODI_SNOW_COVER_1LAYER
USE MODI_URBAN_LW_COEF
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT_CANYON  ! canyon air temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PQ_CANYON  ! canyon air specific humidity
REAL, DIMENSION(:),   INTENT(IN)    :: PU_CANYON  ! canyon hor. wind
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_ROOF   ! roof surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_ROAD   ! road surface temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PTS_WALL   ! wall surface temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PTSNOW_ROOF ! snow layers temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PRSNOW_ROOF ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROOF ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROOF ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROOF! snow surface temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PTSNOW_ROAD ! snow layers temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PRSNOW_ROAD ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROAD ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROAD ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROAD! snow surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PPS      ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PTA      ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA      ! specific humidity
                                              ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD    ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA    ! air density at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD  ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PSR      ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PZREF    ! reference height of the first
                                              ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF    ! reference height of the first
                                              ! atmospheric level (wind)
                                              ! at first atmospheric level
REAL,               INTENT(IN)    :: PTSTEP   ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT   ! buildings h
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD ! road emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD  ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL  ! wall sky view factor
!
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROOF ! SW absorbed by roof snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROOF ! LW absorbed by roof snow
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD          ! snow-free    road frac.
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD          ! snow-covered road frac.
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROAD ! SW absorbed by road snow
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROAD ! LW absorbed by road snow
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROOF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROAD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROAD   ! snow melt
!
REAL, DIMENSION(:), INTENT(OUT)   :: PLW_S_TO_N   ! coefficient for sky LW
!                                                 ! contribution to snow on road
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW1_ROAD   ! independant from
REAL, DIMENSION(SIZE(PTA)) :: ZLW1_ROOF   ! surface temperature
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW2_ROAD   ! to be multiplied by
REAL, DIMENSION(SIZE(PTA)) :: ZLW2_ROOF   ! 4th power of
!                                         ! surface temperature
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_W
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_W
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_W
REAL, DIMENSION(SIZE(PTA)) :: ZLW_N_TO_W
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_N
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_N
REAL, DIMENSION(SIZE(PTA)) :: ZLW_N_TO_N

REAL, DIMENSION(SIZE(PTA)) :: ZSR_ROOF    ! snow fall on roof snow (kg/s/m2 of snow)
REAL, DIMENSION(SIZE(PTA)) :: ZSR_ROAD    ! snow fall on road snow (kg/s/m2 of snow)
REAL, DIMENSION(SIZE(PTA)) :: UCAN_MOD    ! module of canyon wind

REAL  :: WCRN_SMALLER ,  ANS_BIGGER  
!
! flags to call to snow routines
!
LOGICAL :: GSNOW_ROOF, GSNOW_ROAD
!
!
!-------------------------------------------------------------------------------
!
PRNSNOW_ROOF(:)=0.
PHSNOW_ROOF (:)=0.
PLESNOW_ROOF(:)=0.
PGSNOW_ROOF (:)=0.
PMELT_ROOF  (:)=0.
PRNSNOW_ROAD(:)=0.
PHSNOW_ROAD (:)=0.
PLESNOW_ROAD(:)=0.
PGSNOW_ROAD (:)=0.
PMELT_ROAD  (:)=0.
!
PABS_LW_SNOW_ROOF(:) = 0.
PABS_LW_SNOW_ROAD(:) = 0.
!
PLW_S_TO_N(:) = 0.
!
!-------------------------------------------------------------------------------
!
!!$GSNOW_ROOF = ANY( PSR(:)>0. .OR. PWSNOW_ROOF(:)>SWE_CRIT )
!!$GSNOW_ROAD = ANY( PSR(:)>0. .OR. PWSNOW_ROAD(:)>SWE_CRIT )
GSNOW_ROOF = .true.
GSNOW_ROAD = .true.

!
!-------------------------------------------------------------------------------
!
!*      5.     Snow mantel model
!              -----------------
!
!*      5.1    roofs
!              -----
!
IF ( GSNOW_ROOF ) THEN
!
!* initializes LW radiative coefficients
!
  ZLW1_ROOF(:) =   PESNOW_ROOF(:) * PLW_RAD(:)
  ZLW2_ROOF(:) = - PESNOW_ROOF(:) * XSTEFAN
!
!* LW absorbed by snow on roof
!
  PABS_LW_SNOW_ROOF(:) = ZLW1_ROOF(:)+ZLW2_ROOF(:)*PTSSNOW_ROOF(:)**4 
!
!* The global amount of snow on roofs is supposed located on a
!  fraction of the roof surface. All computations are then
!  done only for each m2 of snow, and not for each m2 of roof.
!
 
  WHERE (PDN_ROOF(:)>0.) PWSNOW_ROOF(:) = PWSNOW_ROOF(:) / PDN_ROOF(:)
 
  ZSR_ROOF=0.
  WHERE (PDN_ROOF(:)>0.) ZSR_ROOF   (:) = PSR   (:) / PDN_ROOF(:)
!
!* call to snow mantel scheme
!

  WCRN_SMALLER= XWCRN * 0.1
 
  ANS_BIGGER = XANS_TODRY * 5.
  CALL SNOW_COVER_1LAYER(PTSTEP, XANSMIN, XANSMAX,ANS_BIGGER,             &
                         XRHOSMIN, XRHOSMAX, 0.01, .TRUE.,                    &
                         0.,WCRN_SMALLER,                                     &
                         XZ0SN,XZ0HSN,                                        &
                         PTSNOW_ROOF, PASNOW_ROOF,                            &
                         PRSNOW_ROOF, PWSNOW_ROOF, PTSSNOW_ROOF,              &
                         PESNOW_ROOF,                                         &
                         PTS_ROOF, PABS_SW_SNOW_ROOF,                         &
                         ZLW1_ROOF, ZLW2_ROOF,                                &
                         PTA, PQA, PVMOD, PPS, PRHOA, ZSR_ROOF, PZREF, PUREF, &
                         PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,&
                         PMELT_ROOF                                           )
!
!* The global amount of snow on roofs is reported to total roof surface.
!
  
  PWSNOW_ROOF(:) = PWSNOW_ROOF(:) * PDN_ROOF(:)
  
!           
END IF
!
!*      5.2    roads
!              -----
!
IF ( GSNOW_ROAD ) THEN
!
!* call to urban_lw_coef by inverting the snow-free and snow-covered
!  characteristics (emissivities, fractions, outputs)
!
  CALL URBAN_LW_COEF(PESNOW_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,       &
                     PDF_ROAD, PDN_ROAD, PEMIS_ROAD,                      &
                     ZLW_W_TO_W, ZLW_N_TO_W, ZLW_W_TO_N, ZLW_N_TO_N,      &
                     ZLW_S_TO_W, PLW_S_TO_N, ZLW_R_TO_W, ZLW_R_TO_N       )
!
  ZLW1_ROAD(:) =  ZLW_W_TO_N(:) * PTS_WALL(:)**4 &
                + ZLW_R_TO_N(:) * PTS_ROAD(:)**4 &
                + PLW_S_TO_N(:) * PLW_RAD (:)
  ZLW2_ROAD(:) =  ZLW_N_TO_N(:)
!
!* LW absorbed by snow on road
!
  PABS_LW_SNOW_ROAD(:) =  PLW_S_TO_N(:)*PLW_RAD     (:)             &
                        + ZLW_N_TO_N(:)*PTSSNOW_ROAD(:)**4          &
                        + ZLW_W_TO_N(:)*PTS_WALL    (:)**4          &
                        + ZLW_R_TO_N(:)*PTS_ROAD    (:)**4
!
!* The global amount of snow on roads is supposed located on a
!  fraction of the road surface. All computations are then
!  done only for each m2 of snow, and not for each m2 of road.
!
  
  WHERE (PDN_ROAD(:)>0.) PWSNOW_ROAD(:) = PWSNOW_ROAD(:) / PDN_ROAD(:)
  
  ZSR_ROAD=0.
  WHERE (PDN_ROAD(:)>0.) ZSR_ROAD   (:) = PSR   (:) / PDN_ROAD(:)
!
!
!* call to snow mantel scheme
!
  
  UCAN_MOD(:)=MAX(ABS(PU_CANYON(:)),1.0)

  WCRN_SMALLER = XWCRN * 0.1 
  ANS_BIGGER = XANS_TODRY  * 10.
  CALL SNOW_COVER_1LAYER(PTSTEP, 0.2, XANSMAX, ANS_BIGGER,                &
                         XRHOSMIN, XRHOSMAX, 0.1, .FALSE.,                    &
                         1., WCRN_SMALLER,                                    &
                         XZ0SN,XZ0HSN,                                        &
                         PTSNOW_ROAD, PASNOW_ROAD,                            &
                         PRSNOW_ROAD, PWSNOW_ROAD, PTSSNOW_ROAD,              &
                         PESNOW_ROAD,                                         &
                         PTS_ROAD, PABS_SW_SNOW_ROAD, ZLW1_ROAD, ZLW2_ROAD,   &
                         PT_CANYON, PQ_CANYON,UCAN_MOD, PPS, PRHOA,         &
                         ZSR_ROAD, PBLD_HEIGHT/2., PBLD_HEIGHT/2.,            &
                         PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,&
                         PMELT_ROAD                                           )
!
!* The global amount of snow on roads is reported to total road surface.
!
  PWSNOW_ROAD(:) = PWSNOW_ROAD(:) * PDN_ROAD(:)
!
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_SNOW_EVOL
