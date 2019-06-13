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
    MODULE MODI_TEB2
!
!
!
INTERFACE
!
!
    SUBROUTINE TEB2 (PT_CANYON, PQ_CANYON, PU_CANYON,                         &
                     PTI_BLD,                                                 &
                     PT_ROOF, PT_ROAD, PT_WALL, PWS_ROOF,PWS_ROAD,            &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA,                 &
                     PLW_RAD, PDIR_SW_RAD, PSCA_SW_RAD, PZENITH,              &
                     PRR, PSR,                                                &
                     PZREF, PUREF, PVMOD, PVDIR,                              &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PTSTEP,                                                  &
                     PZ0_TOWN, PZ0_ROOF, PZ0_ROAD,                            &
                     PBLD,PBLD_HEIGHT,PWALL_O_HOR,PCAN_HW_RATIO,              &
                     PALB_ROOF, PEMIS_ROOF,                                   &
                     PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
                     PALB_ROAD, PEMIS_ROAD, PSVF_ROAD,                        &
                     PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
                     PALB_WALL, PEMIS_WALL, PSVF_WALL, PLAT,                  &
                     PHC_WALL,PTC_WALL,PD_WALL,                               &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLES_ROOF, PGFLUX_ROOF,     &
                     PRUNOFF_ROOF,                                            &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLES_ROAD, PGFLUX_ROAD,     &
                     PRUNOFF_ROAD,                                            &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,    &
                     PRUNOFF_TOWN, PUSTAR_TOWN, PCH_TOWN, PRI_TOWN,           &
                     PTS_TOWN, PEMIS_TOWN, PDIR_ALB_TOWN, PSCA_ALB_TOWN, PRESA_TOWN , &
                     PTRAD_IN, PTRAD_SUN, PTRAD_SHADE, PTRAD_RFSUN, PTRAD_RFSHADE,    &
            PTGLOBE_SUN, PTGLOBE_SHADE, PTGLOBE_RFSUN, PTGLOBE_RFSHADE,       &
            PTWETB,PTWETB_ROOF,                                               &
            PUTCI_IN, PUTCI_OUTSUN, PUTCI_OUTSHADE, PUTCI_RFSUN, PUTCI_RFSHADE, &
            WBGT_SUN, WBGT_SHADE, WBGT_RFSUN, WBGT_RFSHADE,                     &
            PUTCIC_IN,PUTCIC_OUTSUN,PUTCIC_OUTSHADE,PUTCIC_RFSUN,PUTCIC_RFSHADE &
            ,PTRFZT,PTRDZT,PURDZU                                &
           ,PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,PQ8,PQ9,PQ10,PQ11,PQ12,PQ13             )
!
!
!*      0.1    declarations of arguments 
!
!
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON     ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON     ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(OUT)   :: PU_CANYON     ! canyon hor. wind
REAL, DIMENSION(:), INTENT(INOUT) :: PTI_BLD       ! inside building temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF     ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD     ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL     ! wall layers temperatures
REAL, DIMENSION(:),   INTENT(INOUT) :: PWS_ROOF    ! roof water reservoir
REAL, DIMENSION(:),   INTENT(INOUT) :: PWS_ROAD    ! road water reservoir
REAL, DIMENSION(:),   INTENT(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSNOW_ROOF ! snow layers temperature
REAL, DIMENSION(:),   INTENT(INOUT) :: PRSNOW_ROOF ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROOF ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROOF ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROOF! snow surface temperature
REAL, DIMENSION(:),   INTENT(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSNOW_ROAD ! snow layers temperature
REAL, DIMENSION(:),   INTENT(INOUT) :: PRSNOW_ROAD ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROAD ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROAD ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROAD! snow surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PPS           ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PPA           ! pressure at the first atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS         ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA           ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA           ! specific humidity
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD         ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PVDIR         ! direction of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA         ! exner function
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA         ! air density
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD       ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW_RAD   ! incoming direct solar radiation
                                                   ! on an horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH       ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW_RAD   ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PRR           ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PSR           ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY   ! anthropogenic sensible
!                                                  ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PZREF         ! reference height of the first
                                                   ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF         ! reference height of the first
                                                   ! atmospheric level (wind)
REAL,               INTENT(IN)    :: PTSTEP        ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN      ! town roughness length
                                                   ! for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ROOF      ! roof roughness length
                                                   ! for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ROAD      ! road roughness length
                                                   ! for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PBLD          ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT   ! buildings h
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO ! canyon    h/W
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROOF     ! roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF    ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF      ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF      ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF       ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROAD     ! road albedo
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD    ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROAD      ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROAD      ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROAD       ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD     ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_WALL     ! wall albedo
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL    ! wall emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL      ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL      ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL       ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL     ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PLAT          ! latitude

REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF      ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROOF       ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF      ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLES_ROOF     ! latent heat flux over roof (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF   ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROOF  ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD      ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROAD       ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD      ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLES_ROAD     ! latent heat flux over road (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD   ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD  ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL      ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL       ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL      ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL   ! flux through the wall
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROOF  ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROOF   ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROOF  ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROOF   ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROOF    ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROAD  ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROAD   ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROAD  ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROAD   ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROAD    ! snow melt
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_TOWN      ! net radiation over town
REAL, DIMENSION(:), INTENT(OUT)   :: PH_TOWN       ! sensible heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_TOWN      ! latent heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_TOWN   ! flux through the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PEVAP_TOWN    ! evaporation flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_TOWN  ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN   ! friction velocity over town
REAL, DIMENSION(:), INTENT(OUT)   :: PCH_TOWN      ! town averaged heat transfer
!                                                  ! coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PRI_TOWN      ! town averaged Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PTS_TOWN      ! town surface temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIS_TOWN    ! town equivalent emissivity
REAL, DIMENSION(:), INTENT(OUT)   :: PDIR_ALB_TOWN ! town equivalent direct albedo
REAL, DIMENSION(:), INTENT(OUT)   :: PSCA_ALB_TOWN ! town equivalent diffuse albedo
REAL, DIMENSION(:), INTENT(OUT)   :: PRESA_TOWN    ! town aerodynamical resistance
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTRAD_IN       ! body MRT inside building (K)
REAL, DIMENSION(:), INTENT(OUT)   :: PTRAD_SUN      ! body MRT in the exposed street (K)
REAL, DIMENSION(:), INTENT(OUT)   :: PTRAD_SHADE    ! body MRT in the shaded street (K)
REAL, DIMENSION(:), INTENT(OUT)   :: PTRAD_RFSUN      ! body MRT on the exposed roof (K)
REAL, DIMENSION(:), INTENT(OUT)   :: PTRAD_RFSHADE    ! body MRT on the shaded roof (K)
REAL, DIMENSION(:), INTENT(OUT)   :: PTWETB          ! wet-bulb temperature over the ground (K)
REAL, DIMENSION(:), INTENT(OUT)   :: PTWETB_ROOF     ! wet-bulb temperature over the roof (K)

REAL, DIMENSION(:), INTENT(OUT)   :: PUTCI_IN       ! UTCI inside building
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCI_OUTSUN   ! UTCI in the exposed street
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCI_OUTSHADE ! UTCI in the shaded street
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCI_RFSUN    ! UTCI on the exposed roof
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCI_RFSHADE  ! UTCI on the shaded roof

REAL, DIMENSION(:), INTENT(OUT)   :: WBGT_SUN         ! WBGT  wet bulb globe temperature in the street (C)
REAL, DIMENSION(:), INTENT(OUT)   :: WBGT_SHADE       ! WBGT  wet bulb globe temperature in the street (C)
REAL, DIMENSION(:), INTENT(OUT)   :: WBGT_RFSUN       ! WBGT  wet bulb globe temperature on the roof (C)
REAL, DIMENSION(:), INTENT(OUT)   :: WBGT_RFSHADE     ! WBGT  wet bulb globe temperature on the roof (C)

REAL, DIMENSION(:), INTENT(OUT)   :: PUTCIC_IN       ! UTCI inside building
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCIC_OUTSUN   ! UTCI in the exposed street
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCIC_OUTSHADE ! UTCI in the shaded street
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCIC_RFSUN   ! UTCI in the exposed street
REAL, DIMENSION(:), INTENT(OUT)   :: PUTCIC_RFSHADE ! UTCI in the shaded street
!
 REAL, DIMENSION(:), INTENT(OUT)   :: PTGLOBE_SUN      ! Globe Temperature in the exposed street (K)
 REAL, DIMENSION(:), INTENT(OUT)   :: PTGLOBE_SHADE    ! Globe Temperature in the shaded street (K)
 REAL, DIMENSION(:), INTENT(OUT)   :: PTGLOBE_RFSUN    ! Globe Temperature  on the exposed roof (K)
 REAL, DIMENSION(:), INTENT(OUT)   :: PTGLOBE_RFSHADE  ! Globe Temperature  on the shaded roof  (K)
!
 REAL, DIMENSION(:), INTENT(OUT)   :: PTRFZT
 REAL, DIMENSION(:), INTENT(OUT)   :: PTRDZT
 REAL, DIMENSION(:), INTENT(OUT)   :: PURDZU
!
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ1
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ2
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ3
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ4
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ5
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ6
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ7
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ8
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ9
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ10
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ11
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ12
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ13
!
!
END SUBROUTINE TEB2
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_TEB2
