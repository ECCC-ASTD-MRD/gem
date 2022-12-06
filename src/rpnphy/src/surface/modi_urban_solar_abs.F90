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
    MODULE MODI_URBAN_SOLAR_ABS
!
!
!
INTERFACE
!
!
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
                               PREF_SW_SNOW_ROOF, PREF_SW_SNOW_ROAD          )
!
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
!                                                      ! by snow-free roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_SW_ROAD      ! solar radiation absorbed
!                                                      ! by snow-free roads
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
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_SNOW_ROOF ! reflected solar radiations
REAL, DIMENSION(:), INTENT(OUT)   :: PREF_SW_SNOW_ROAD ! by snow over road and roof
!
END SUBROUTINE URBAN_SOLAR_ABS
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_URBAN_SOLAR_ABS
