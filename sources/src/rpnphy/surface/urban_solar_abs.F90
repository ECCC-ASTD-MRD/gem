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
    SUBROUTINE URBAN_SOLAR_ABS(PDIR_SW, PSCA_SW, PZENITH,                    &
                               PBLD, PWALL_O_HOR, PCAN_HW_RATIO,             &
                               PALB_ROOF,                                    &
                               PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
                               PASNOW_ROOF, PASNOW_ROAD,                     &
                               PDN_ROOF, PDF_ROOF, PDN_ROAD, PDF_ROAD,       &
                               PABS_SW_ROOF, PABS_SW_ROAD, PABS_SW_WALL,     &
                               PABS_SW_SNOW_ROOF, PABS_SW_SNOW_ROAD,         &
                               PDIR_ALB_TOWN, PSCA_ALB_TOWN ,                &
                               PREF_SW_ROOF, PREF_SW_ROAD, PREF_SW_WALL,     &
                               PREF_SW_SNOW_ROOF, PREF_SW_SNOW_ROAD         )
!   ##########################################################################
!
!!****  *URBAN_SOLAR_ABS*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the solar radiation flux absorbed by roofs, roads and walls.
!     The absorption by roofs is trivial.
!         
!     
!!**  METHOD
!     ------
!
!
!        computation of input solar radiation on each surface
!        ****************************************************
!
!    direct fluxes:
!    -------------
!
!    dir_Rg_road (Wm-2) =   S * 2*theta0/pi
!                         - S *2/tan(zen) * h/W /pi * (1-cos(theta0))
!
!    dir_Rg_wall (Wm-2) =   S / tan(zen) /pi * (1-cos(theta0))
!                         + S * W/h * (1/2 -theta0/pi)
!
!   where zen      is the zenithal angle, from horizon
!         h/W      is the aspect ratio of the canyon
!         S        is the direct solar radiation flux on a horizontal surface
!
!         theta0 = arcsin(min(W/h * tan(zen),1))
!
!   The surfaces will keep (1-a) times these fluxes, and reflect the
!   remaining
!
!    scattered fluxes:
!    ----------------
!
!   sca_Rg_road = sca_Rg * SVF_road
!
!   sca_Rg_wall = sca_Rg * SVF_wall
!
!
!    solar flux and isotropic reflections :
!    ------------------------------------
!
!  after 0 reflection, the absorbed part of the flux is:
!
!      ARg_r(0) = (1-a_r) (sca_Rg_road + dir_Rg_road)
!
!      ARg_w(0) = (1-a_w) (sca_Rg_wall + dir_Rg_wall)
!  
!    and the reflected parts are
!
!      RRg_r(0) = a_r (sca_Rg_road + dir_Rg_road)
!
!      RRg_w(0) = a_w (sca_Rg_wall + dir_Rg_wall)
!
!  after n reflection:
!
!      ARg_r(n) = ARg_r(n-1) + RRg_w(n-1) * (1-  SVF_r)(1-a_r)
!
!      ARg_w(n) = ARg_w(n-1) + RRg_r(n-1) *      SVF_w (1-a_w)
!                            + RRg_w(n-1) * (1-2*SVF_w)(1-a_w)
!
!      RRg_r(n) = (1- SVF_r) a_r RRg_w(n-1)
!
!      RRg_w(n) =     SVF_w  a_w RRg_r(n-1)
!                +(1-2SVF_w) a_w RRg_w(n-1)
!
!
!   i.e.
!                                               n-1
!      ARg_r(n) = ARg_r(0) + (1-  SVF_r)(1-a_r) SUM RRg_w(k)
!                                               k=0
!
!                                               n-1
!      ARg_w(n) = ARg_w(0) +      SVF_w (1-a_w) SUM RRg_r(k)
!                                               k=0
!                                               n-1
!                          + (1-2*SVF_w)(1-a_w) SUM RRg_w(k)
!                                               k=0
!
! with
!
!     n                             n-1
!    SUM RRg_r(k) = (1-  SVF_r) a_r SUM RRg_w(k)      +  RRg_r(0)
!    k=0                            k=0
!
!     n                             n-1
!    SUM RRg_w(k) =      SVF_w  a_w SUM RRg_r(k) 
!    k=0                            k=0
!                                   n-1
!                  +(1-2*SVF_w) a_w SUM RRg_w(k)      +  RRg_w(0)
!                                   k=0
!
!
!   Then
!
!     n                                        n-1
!    SUM RRg_w(k) =  (1-2*SVF_w)       a_w     SUM RRg_w(k)
!    k=0                                       k=0
!                                              n-2
!                  + (1-  SVF_r) SVF_w a_w a_r SUM RRg_w(k) 
!                                              k=0
!
!                  + RRg_w(0) + SVF_w a_w RRg_r(0)
!
!
!
!
!  solving this system, lead after an infinity of reflections/absorptions:
!
!    inf                      RRg_w(0) + SVF_w a_w RRg_r(0)
!    SUM RRg_w(k) = ----------------------------------------------------
!    k=0             1 - (1-2*SVF_w) a_w - (1-  SVF_r) SVF_w a_w a_r
!
!
!    inf            (1-  SVF_r) a_r ( a_w SVF_w RRg_r(0) + RRg_w(0) )
!    SUM RRg_r(k) = ------------------------------------------------------------ + RRg_r(0)
!    k=0             1 - (1-2*SVF_w) a_w - (1-  SVF_r) SVF_w a_w a_r
!
!
! ARg_r(n) and ARg_w(n) follow
!
!
! If snow is present, the albedos in all these formulae (and only these,
! not the final net radiation budget) are modified by the albedo of the
! snow-covered surface.
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
!!                  21/11/00 (V. Masson)  bug in reflections for roads
!!                     12/02 (A. Lemonsu) bug in diagnostic of albedo
!!                                        bug in comp. of ztanzen
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS, ONLY : XPI
!
implicit none
#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW           ! incoming direct solar radiation
                                                       ! on an horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH           ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW           ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PBLD              ! buildings fraction
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO     ! canyon    h/W
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR       ! wall surf. / hor. surf
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROOF         ! roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROAD         ! road albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD         ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_WALL         ! wall albedo
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL         ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PASNOW_ROOF       ! roof snow albedo
REAL, DIMENSION(:), INTENT(IN)    :: PASNOW_ROAD       ! road snow albedo
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF          ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD          ! snow-covered road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD          ! snow-free    road fraction
!
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_SW_ROOF      ! solar radiation absorbed
!                                                      ! by roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_SW_ROAD      ! solar radiation absorbed
!                                                      ! by roads
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_SW_WALL      ! solar radiation absorbed
!                                                      ! by walls
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_SW_SNOW_ROOF ! solar radiation absorbed
!                                                      ! by snow-covered roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_SW_SNOW_ROAD ! solar radiation absorbed
!                                                      ! by snow-covered roads
REAL, DIMENSION(:), INTENT(OUT)   :: PDIR_ALB_TOWN     ! town direct albedo
REAL, DIMENSION(:), INTENT(OUT)   :: PSCA_ALB_TOWN     ! town diffuse albedo
!
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_ROAD      ! reflected solar radiations
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_WALL      ! by roads and walls and roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_ROOF      ! 
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_SNOW_ROOF      ! reflected solar radiations
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_SNOW_ROAD      ! by snow over road and roof
!
!*      0.2    declarations of local variables
!
!                                                           
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW       ! direct and diffuse incoming radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW       ! with a minimum to compute albedo
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZTANZEN       ! tangente of solar zenithal angle
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZTHETA0       ! canyon angle for
!                                               ! which solar
!                                               ! radiation
!                                               ! reaches the road
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZAALB_ROAD    ! averaged albedo
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_ROAD  ! direct radiation reaching
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_WALL  ! road and wall
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_ROAD  ! diffuse radiation reaching
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_WALL  ! road and wall
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZROOF_SW   ! roof, wall and
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZWALL_SW   ! road fractions of SW
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZROAD_SW   ! interacting surf.
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZDIR_SW_UP ! upward solar rad.
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZSCA_SW_UP ! from direct or diffuse rad.
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_ROOF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_ROAD      ! absorbed by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_WALL      ! road, wall, and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_SNOW_ROOF ! over roof and wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_DIR_SW_SNOW_ROAD ! coming from direct rad.
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_ROOF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_ROAD      ! absorbed by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_WALL      ! road, wall, and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_SNOW_ROOF ! over roof and wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZABS_SCA_SW_SNOW_ROAD ! coming from diffuse rad.
!
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_DIR_SW_ROOF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_DIR_SW_ROAD      ! reflected by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_DIR_SW_WALL      ! road, wall, and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_DIR_SW_SNOW_ROOF ! over roof and wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_DIR_SW_SNOW_ROAD ! coming from direct rad.
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_SCA_SW_ROOF      ! solar radiation
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_SCA_SW_ROAD      ! reflected by roofs,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_SCA_SW_WALL      ! road, wall, and snow
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_SCA_SW_SNOW_ROOF ! over roof and wall,
REAL, DIMENSION(SIZE(PDIR_SW)) :: ZREF_SCA_SW_SNOW_ROAD ! coming from diffuse rad.
!
INTEGER                        :: JI                    ! loop index
!-------------------------------------------------------------------------------
!
PABS_SW_ROOF(:) = 0.
PABS_SW_ROAD(:) = 0.
PABS_SW_WALL(:) = 0.
PABS_SW_SNOW_ROOF(:) = 0.
PABS_SW_SNOW_ROAD(:) = 0.
!
PREF_SW_ROOF(:) = 0.
PREF_SW_ROAD(:) = 0.
PREF_SW_WALL(:) = 0.
PREF_SW_SNOW_ROOF(:) = 0.
PREF_SW_SNOW_ROAD(:) = 0.
!
ZDIR_SW = MAX(PDIR_SW,1.E-3)
ZSCA_SW = MAX(PSCA_SW,1.E-3)
!
!-------------------------------------------------------------------------------
!
!*      1.     SOLAR RADIATIONS FOR ROOFS
!              --------------------------
!
ZABS_DIR_SW_ROOF     (:) = ZDIR_SW(:) * (1. - PALB_ROOF  (:))
ZABS_DIR_SW_SNOW_ROOF(:) = ZDIR_SW(:) * (1. - PASNOW_ROOF(:))
ZABS_SCA_SW_ROOF     (:) = ZSCA_SW(:) * (1. - PALB_ROOF  (:))
ZABS_SCA_SW_SNOW_ROOF(:) = ZSCA_SW(:) * (1. - PASNOW_ROOF(:))
!
ZREF_DIR_SW_ROOF     (:) = ZDIR_SW(:) * PALB_ROOF  (:)
ZREF_DIR_SW_SNOW_ROOF(:) = ZDIR_SW(:) * PASNOW_ROOF(:)
ZREF_SCA_SW_ROOF     (:) = ZSCA_SW(:) * PALB_ROOF  (:)
ZREF_SCA_SW_SNOW_ROOF(:) = ZSCA_SW(:) * PASNOW_ROOF(:)
!
!-------------------------------------------------------------------------------
!
!*      2.     SOLAR RADIATIONS FOR ROADS AND WALLS
!              ------------------------------------
!
DO JI=1,SIZE(PZENITH(:))
!
  IF (ABS(0.5*XPI-PZENITH(JI)) <  1.E-6) THEN
    IF(0.5*XPI-PZENITH(JI) > 0.)  ZTANZEN(JI)=TAN(0.5*XPI-1.E-6)
    IF(0.5*XPI-PZENITH(JI) <= 0.) ZTANZEN(JI)=TAN(0.5*XPI+1.E-6)
  ELSEIF (ABS(PZENITH(JI)) <  1.E-6) THEN
    ZTANZEN(JI)=SIGN(1.,PZENITH(JI))*TAN(1.E-6)
  ELSE
    ZTANZEN(JI) = TAN(PZENITH(JI))
  ENDIF

ENDDO

!*      2.1    radiation coefficients
!              ----------------------
!
ZTHETA0(:) = ASIN( MIN(ABS( 1./ZTANZEN(:))/PCAN_HW_RATIO(:), 1. ) )
!
!*      2.2    direct solar radiation received by roads
!               ---------------------------------------
!
ZDIR_SW_ROAD(:) =  ZDIR_SW(:) * 2. * ZTHETA0(:) / XPI                  &
                 - ZDIR_SW(:) * 2. * ZTANZEN(:) / XPI                  &
                              * PCAN_HW_RATIO(:) * (1.-COS(ZTHETA0(:)))
!
!*      2.3    direct solar radiation received by walls
!              ----------------------------------------
!
ZDIR_SW_WALL(:) = ( ZDIR_SW(:) - ZDIR_SW_ROAD(:) )  &
                   * 0.5 / PCAN_HW_RATIO(:)
!
!*      2.4    diffuse solar radiation received by roads
!              -----------------------------------------
!
ZSCA_SW_ROAD(:) = ZSCA_SW(:) * PSVF_ROAD(:)
!
!*      2.5    diffuse solar radiation received by walls
!              -----------------------------------------
!
ZSCA_SW_WALL(:) = ZSCA_SW(:) * PSVF_WALL(:)
!
!*      2.6    averaged albedos when snow is present
!              -------------------------------------
!
ZAALB_ROAD(:) =  PDF_ROAD(:) * PALB_ROAD  (:) &
               + PDN_ROAD(:) * PASNOW_ROAD(:)
!
!*      2.8    absorption of direct incoming solar radiation
!              ---------------------------------------------
!
CALL SOLAR_REFLECTIONS(ZDIR_SW_ROAD,ZDIR_SW_WALL,                               &
                       ZABS_DIR_SW_ROAD, ZABS_DIR_SW_SNOW_ROAD,ZABS_DIR_SW_WALL,&
                       ZREF_DIR_SW_ROAD,ZREF_DIR_SW_SNOW_ROAD,ZREF_DIR_SW_WALL  )
!
!
!*      2.9    absorption of diffuse incoming solar radiation
!              ----------------------------------------------
!
CALL SOLAR_REFLECTIONS(ZSCA_SW_ROAD,ZSCA_SW_WALL,                               &
                       ZABS_SCA_SW_ROAD, ZABS_SCA_SW_SNOW_ROAD,ZABS_SCA_SW_WALL,&
                       ZREF_SCA_SW_ROAD,ZREF_SCA_SW_SNOW_ROAD,ZREF_SCA_SW_WALL )
!
!-------------------------------------------------------------------------------
!
!*      3.     Town albedo
!              -----------
!
!*      3.1    direct albedo
!              -------------
!
CALL TOWN_ALBEDO(ZDIR_SW,ZABS_DIR_SW_ROOF,ZABS_DIR_SW_SNOW_ROOF,            &
                 ZABS_DIR_SW_ROAD, ZABS_DIR_SW_SNOW_ROAD,ZABS_DIR_SW_WALL,  &
                 PDIR_ALB_TOWN                                              )
!
!*      3.2    direct albedo
!              -------------
!
CALL TOWN_ALBEDO(ZSCA_SW,ZABS_SCA_SW_ROOF,ZABS_SCA_SW_SNOW_ROOF,            &
                 ZABS_SCA_SW_ROAD, ZABS_SCA_SW_SNOW_ROAD,ZABS_SCA_SW_WALL,  &
                 PSCA_ALB_TOWN                                              )
!
!-------------------------------------------------------------------------------
!
!*      4.     Trivial cases
!              -------------
!
WHERE(PDIR_SW(:)==0.)
  ZABS_DIR_SW_ROOF     (:) = 0.
  ZABS_DIR_SW_ROAD     (:) = 0.
  ZABS_DIR_SW_WALL     (:) = 0.
  ZABS_DIR_SW_SNOW_ROOF(:) = 0.
  ZABS_DIR_SW_SNOW_ROAD(:) = 0.
  ZREF_DIR_SW_ROOF     (:) = 0.
  ZREF_DIR_SW_ROAD     (:) = 0.
  ZREF_DIR_SW_WALL     (:) = 0.
  ZREF_DIR_SW_SNOW_ROOF(:) = 0.
  ZREF_DIR_SW_SNOW_ROAD(:) = 0.
END WHERE
!
WHERE(PSCA_SW(:)==0.)
  ZABS_SCA_SW_ROOF     (:) = 0.
  ZABS_SCA_SW_ROAD     (:) = 0.
  ZABS_SCA_SW_WALL     (:) = 0.
  ZABS_SCA_SW_SNOW_ROOF(:) = 0.
  ZABS_SCA_SW_SNOW_ROAD(:) = 0.
  ZREF_SCA_SW_ROOF     (:) = 0.
  ZREF_SCA_SW_ROAD     (:) = 0.
  ZREF_SCA_SW_WALL     (:) = 0.
  ZREF_SCA_SW_SNOW_ROOF(:) = 0.
  ZREF_SCA_SW_SNOW_ROAD(:) = 0.
END WHERE
!
!-------------------------------------------------------------------------------
!
!*      5.     Total solar radiation absorbed by each surface
!              ----------------------------------------------
!
! solar radiation absorbed by roofs
!
PABS_SW_ROOF     (:) = ZABS_DIR_SW_ROOF     (:) + ZABS_SCA_SW_ROOF     (:)
!
! solar radiation absorbed by roads
!
PABS_SW_ROAD     (:) = ZABS_DIR_SW_ROAD     (:) + ZABS_SCA_SW_ROAD     (:)
!
! solar radiation absorbed by walls
!
PABS_SW_WALL     (:) = ZABS_DIR_SW_WALL     (:) + ZABS_SCA_SW_WALL     (:)
!
! solar radiation absorbed by snow on roofs
!
PABS_SW_SNOW_ROOF(:) = ZABS_DIR_SW_SNOW_ROOF(:) + ZABS_SCA_SW_SNOW_ROOF(:)
!
! solar radiation absorbed by snow on roads
!
PABS_SW_SNOW_ROAD(:) = ZABS_DIR_SW_SNOW_ROAD(:) + ZABS_SCA_SW_SNOW_ROAD(:)
!
!-------------------------------------------------------------------------------
!
!*      6.     Total solar radiation reflected by each surface
!              ----------------------------------------------
!
! solar radiation reflected by roofs
!
PREF_SW_ROOF     (:) = ZREF_DIR_SW_ROOF     (:) + ZREF_SCA_SW_ROOF     (:)
!
! solar radiation reflected by roads
!
PREF_SW_ROAD     (:) = ZREF_DIR_SW_ROAD     (:) + ZREF_SCA_SW_ROAD     (:)
!
! solar radiation reflected by walls
!
PREF_SW_WALL     (:) = ZREF_DIR_SW_WALL     (:) + ZREF_SCA_SW_WALL     (:)
!
! solar radiation reflected by snow on roofs
!
PREF_SW_SNOW_ROOF(:) = ZREF_DIR_SW_SNOW_ROOF(:) + ZREF_SCA_SW_SNOW_ROOF(:)
!
! solar radiation reflected by snow on roads
!
PREF_SW_SNOW_ROAD(:) = ZREF_DIR_SW_SNOW_ROAD(:) + ZREF_SCA_SW_SNOW_ROAD(:)
!
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
SUBROUTINE SOLAR_REFLECTIONS(ZSW_ROAD,ZSW_WALL,ZABS_SW_ROAD,ZABS_SW_SNOW,ZABS_SW_WALL, &
			     ZSREF_SW_ROAD,ZSREF_SW_SNOW,ZSREF_SW_WALL)
!
REAL, DIMENSION(:), INTENT(IN) :: ZSW_ROAD     ! solar radiation received by road
REAL, DIMENSION(:), INTENT(IN) :: ZSW_WALL     ! and wall before reflection
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_ROAD ! solar radiation absorbed by
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_SNOW ! solar radiation absorbed by
REAL, DIMENSION(:), INTENT(OUT):: ZABS_SW_WALL ! road, snow over road, and wall 
!
REAL, DIMENSION(:), INTENT(OUT):: ZSREF_SW_ROAD ! sum of all reflections
REAL, DIMENSION(:), INTENT(OUT):: ZSREF_SW_SNOW ! sum of all reflections
REAL, DIMENSION(:), INTENT(OUT):: ZSREF_SW_WALL ! against road and wall
!
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZREF0_SW_ROAD ! first solar reflection
REAL, DIMENSION(SIZE(ZSW_ROAD)) :: ZREF0_SW_WALL ! against road and wall
!
!
!*      A.     first solar radiation reflection
!              --------------------------------
!
!
  ZREF0_SW_ROAD(:) = ZAALB_ROAD(:) * ZSW_ROAD(:)
!
  ZREF0_SW_WALL(:) = PALB_WALL (:) * ZSW_WALL(:)
!
!*      B.     sum of solar radiation reflected
!              --------------------------------
!

  ZSREF_SW_WALL(:) = (                                  ZREF0_SW_WALL(:) &
                        + PSVF_WALL(:) * PALB_WALL(:) * ZREF0_SW_ROAD(:))&
                     / (1. - (1.-2.*PSVF_WALL(:)) * PALB_WALL(:)         &
                           - (1.-   PSVF_ROAD(:)) * PSVF_WALL(:)         &
                                                  * PALB_WALL(:)         &
                                                  * ZAALB_ROAD(:)        )
!
  ZSREF_SW_ROAD(:) = ( (1.-   PSVF_ROAD(:)) * ZAALB_ROAD(:)              &
                                            * ZREF0_SW_WALL(:)           &
                      +(1.-   PSVF_ROAD(:)) * ZAALB_ROAD(:)              &
                                            * PSVF_WALL(:)               &
                                            * PALB_WALL(:)               &
                                            * ZREF0_SW_ROAD(:)          )&
                     / (1. - (1.-2.*PSVF_WALL(:)) * PALB_WALL(:)         &
                           - (1.-   PSVF_ROAD(:)) * PSVF_WALL(:)         &
                                                  * PALB_WALL(:)         &
                                                  * ZAALB_ROAD(:)       )&
                    + ZREF0_SW_ROAD(:)
!
!
!*      C.     total solar radiation received by roads
!              ---------------------------------------
!
  ZABS_SW_ROAD(:) = (1.-PALB_ROAD(:))                              &
                    * (   ZSW_ROAD(:)                              &
                        + ZSREF_SW_WALL(:) * (1.-   PSVF_ROAD(:)) )
!
  ZABS_SW_SNOW(:) = (1.-PASNOW_ROAD  (:))                          &
                    * (   ZSW_ROAD(:)                              &
                        + ZSREF_SW_WALL(:) * (1.-   PSVF_ROAD(:)) )
!
!
!*      D.     total solar radiation received by walls
!              ---------------------------------------
!
  ZABS_SW_WALL(:) = (1.-PALB_WALL(:))                              &
                    * (   ZSW_WALL(:)                              &
                        + ZSREF_SW_ROAD(:) *        PSVF_WALL(:)   &
                        + ZSREF_SW_WALL(:) * (1.-2.*PSVF_WALL(:)) )
!
!*      E.     total solar radiation reflected by roads
!              ---------------------------------------
!
  ZSREF_SW_ROAD(:) = PALB_ROAD(:)                                   &
                    * (   ZSW_ROAD(:)                              &
                        + ZSREF_SW_WALL(:) * (1.-   PSVF_ROAD(:)) )
!
  ZSREF_SW_SNOW(:) = PASNOW_ROAD  (:)                               &
                    * (   ZSW_ROAD(:)                              &
                        + ZSREF_SW_WALL(:) * (1.-   PSVF_ROAD(:)) )
!
!
!*      F.     total solar radiation reflected by walls
!              ---------------------------------------
!
  ZSREF_SW_WALL(:) = PALB_WALL(:)                                   &
                    * (   ZSW_WALL(:)                              &
                        + ZSREF_SW_ROAD(:) *        PSVF_WALL(:)   &
                        + ZSREF_SW_WALL(:) * (1.-2.*PSVF_WALL(:)) )
!
END SUBROUTINE SOLAR_REFLECTIONS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TOWN_ALBEDO(ZSW,ZABS_SW_ROOF,ZABS_SW_SNOW_ROOF,            &
                       ZABS_SW_ROAD, ZABS_SW_SNOW_ROAD,ZABS_SW_WALL,  &
                       ZALBEDO                                        )
!
REAL, DIMENSION(:), INTENT(IN) :: ZSW               ! incoming solar radiation
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_ROOF      ! solar radiation absorbed by roofs
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_ROAD      ! solar radiation absorbed by roads
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_WALL      ! solar radiation absorbed by walls
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_SNOW_ROOF ! solar radiation absorbed by roof snow
REAL, DIMENSION(:), INTENT(IN) :: ZABS_SW_SNOW_ROAD ! solar radiation absorbed by road snow
REAL, DIMENSION(:), INTENT(OUT):: ZALBEDO           ! town averaged albedo

REAL, DIMENSION(SIZE(ZSW))     :: ZSW_UP            ! outgoing solar radiation

ZSW_UP(:) = ZSW(:)                                          &
          - ( PBLD(:)*PDF_ROOF(:)*ZABS_SW_ROOF     (:)      &
             +PBLD(:)*PDN_ROOF(:)*ZABS_SW_SNOW_ROOF(:)      &
             +(1.-PBLD(:))*PDF_ROAD(:)*ZABS_SW_ROAD     (:) &
             +(1.-PBLD(:))*PDN_ROAD(:)*ZABS_SW_SNOW_ROAD(:) &
             +PWALL_O_HOR(:)*ZABS_SW_WALL(:))
!
ZALBEDO(:)  = ZSW_UP(:) / ZSW(:)
!
END SUBROUTINE TOWN_ALBEDO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_SOLAR_ABS
