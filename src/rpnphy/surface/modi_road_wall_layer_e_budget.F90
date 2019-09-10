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
    MODULE MODI_ROAD_WALL_LAYER_E_BUDGET
!
!
!
INTERFACE
!
!
    SUBROUTINE ROAD_WALL_LAYER_E_BUDGET(PT_ROAD, PT_WALL, PQSAT_ROAD,          &
                                      PT_CANYON, PQ_CANYON,                    &
                                      PTA, PQA, PPS,                           &
                                      PLW_RAD,  PTSTEP,                        &
                                      PH_TRAFFIC, PLE_TRAFFIC,                 &
                                      PBLD, PWALL_O_ROAD,                      &
                                      PEMIS_ROAD, PSVF_ROAD,                   &
                                      PHC_ROAD,PTC_ROAD,PD_ROAD,               &
                                      PEMIS_WALL, PSVF_WALL,                   &
                                      PHC_WALL,PTC_WALL,PD_WALL,               &
                                      PTI_BLD, PAC_BLD,                        &
                                      PDELT_ROAD, PDN_ROAD,                    &
                                      PTSNOW_ROAD, PESNOW_ROAD,                &
                                      PGSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD,  &
                                      PRHOA, PAC_WALL,                         &
                                      PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,         &
                                      PABS_SW_ROAD, PABS_SW_WALL,              &
                                      PABS_LW_ROAD, PABS_LW_WALL,              &
                                      PLW_S_TO_R, PLW_S_TO_W                   )
!
!
!*      0.1    declarations of arguments 
!
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD    ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL    ! wall layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_ROAD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PQ_CANYON    ! and specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTA          ! atmospheric air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQA          ! and specific humidity at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL,                 INTENT(IN)    :: PTSTEP     ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC   ! anthropogenic sensible
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC  ! anthropogenic latent
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_ROAD ! wall Surf. / road Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD   ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROAD     ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROAD     ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROAD      ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD    ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL   ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL     ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL     ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL      ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL    ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD      ! inside building temperature
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
                                                  ! inside the building itself
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD   ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! road snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PESNOW_ROAD  ! road snow emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! road snow conduction
!                                                 ! heat fluxes at mantel
!                                                 ! base
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD  ! snow sensible heat
!                                                 ! fluxes at mantel top
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD ! snow latent heat
!                                                 ! fluxes at mantel top
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD     ! aerodynamical conductance
!                                                 ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT ! aerodynamical conductance
!                                                 ! between road and canyon
!                                                 ! (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_TOP      ! aerodynamical conductance
!                                                 ! between atmosphere and
!                                                 ! canyon top
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD ! absorbed infrared rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL ! absorbed infrared rad.
!
REAL, DIMENSION(:), INTENT(OUT)   :: PLW_S_TO_R    ! 
REAL, DIMENSION(:), INTENT(OUT)   :: PLW_S_TO_W    ! 
!
END SUBROUTINE ROAD_WALL_LAYER_E_BUDGET
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_ROAD_WALL_LAYER_E_BUDGET
