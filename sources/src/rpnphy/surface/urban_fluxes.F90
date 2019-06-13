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
    SUBROUTINE URBAN_FLUXES(                                                  &
                     PTS_TOWN, PEMIS_TOWN,                                    &
                     PT_CANYON, PQ_CANYON,                                    &
                     PTS_ROOF,PTS_ROAD,PTS_WALL,                              &
                     PTA, PQA, PRHOA,                                         &
                     PVMOD,                                                   &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PBLD,                                                    &
                     PEMIS_ROOF,                                              &
                     PESNOW_ROOF,                                             &
                     PLW_RAD, PLW_S_TO_R, PLW_S_TO_W, PLW_S_TO_N,             & 
                     PABS_SW_ROOF, PABS_LW_ROOF, PABS_SW_WALL, PABS_LW_WALL,  &
                     PABS_SW_ROAD, PABS_LW_ROAD,                              &
                     PABS_LW_SNOW_ROOF,                                       &
                     PABS_LW_SNOW_ROAD,                                       &
                     PAC_ROOF, PAC_ROOF_WAT,                                  &
                     PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PCD,                   &
                     PQSAT_ROOF, PQSAT_ROAD,                                  &
                     PDELT_ROOF, PDELT_ROAD,                                  &
                     PROOF, PWALL, PROAD, PTOTS_O_HORS,                       &
                     PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLES_ROOF, PGFLUX_ROOF,     &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLES_ROAD, PGFLUX_ROAD,     &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,    &
                     PUSTAR_TOWN                                              )
!   ##########################################################################
!
!!****  *URBAN_FLUXES* computes fluxes on urbanized surfaces  
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
!!                     12/02 (A. Lemonsu) modifications of emissivity and Tstown
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XLSTT, XSTEFAN
!
USE MODI_URBAN_LW_COEF
!
implicit none
!!!#include <arch_specific.hf>
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTS_TOWN     ! town surface temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIS_TOWN   ! town equivalent emissivity
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON    ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF     ! roof surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL     ! wall surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PTA          ! temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA          ! specific humidity
                                                  ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD        ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! air density
                                                  ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC   ! anthropogenic sensible
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC  ! anthropogenic latent
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY  ! anthropogenic sensible
!                                                 ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY ! anthropogenic latent
!                                                 ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF   ! roof emissivity
!
REAL, DIMENSION(:), INTENT(IN)    :: PESNOW_ROOF  !snow  roof emissivity
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! incoming longwave rad.
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R   ! contribution of LW
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W   ! from sky to Road,
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_N   ! walls, and road snow
!
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF ! absorbed SW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROOF ! absorbed LW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL ! absorbed SW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL ! absorbed LW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed SW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROAD ! absorbed LW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_SNOW_ROOF ! abs. LW rad. by snow
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_SNOW_ROAD ! abs. LW rad. by snow
!
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF     ! surface conductance
!                                                 ! for heat transfers
!                                                 ! above roofs
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT ! surface conductance
!                                                 ! for heat transfers
!                                                 ! above roofs (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! surface conductance
!                                                 ! for heat transfer
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD     ! surface conductance
!                                                 ! for heat transfers
!                                                 ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT ! surface conductance
!                                                 ! for heat transfers
!                                                 ! inside canyon (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PCD          ! drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROOF   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF   ! water fraction on snow-free
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD   ! roof and roads
REAL, DIMENSION(:), INTENT(IN)    :: PROOF        ! 
REAL, DIMENSION(:), INTENT(IN)    :: PWALL        ! roof, wall and road
REAL, DIMENSION(:), INTENT(IN)    :: PROAD        ! fractions of exchange surf.
REAL, DIMENSION(:), INTENT(IN)    :: PTOTS_O_HORS ! total canyon+roof surface
!                                                 ! over horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF     ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF     ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD     ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered road fraction
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLES_ROOF    ! latent heat flux of sublimation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROAD      ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLES_ROAD    ! latent heat flux of sublimation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL     ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL      ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL     ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL  ! flux through the wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! flux under the snow
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_TOWN     ! net radiation over town
REAL, DIMENSION(:), INTENT(OUT)   :: PH_TOWN      ! sensible heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_TOWN     ! latent heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_TOWN  ! flux through the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PEVAP_TOWN   ! evaporation (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friction velocity over town
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLW_RAD))    :: ZLW_UP ! upwards radiations

!-------------------------------------------------------------------------------
!
!*      8.     Fluxes at snow-free roofs
!              -------------------------
!
!                                            net radiation
!
PRN_ROOF(:) = PABS_SW_ROOF(:) + PABS_LW_ROOF(:)
!
!                                            sensible heat flux
!
PH_ROOF(:) =(PTS_ROOF(:) - PTA(:)) &
             * PAC_ROOF(:) * PRHOA(:) * XCPD
!
!                                            latent heat of evaporation from
!                                            the ground
!
PLE_ROOF(:) =(PQSAT_ROOF(:) - PQA(:))                     &
              * PAC_ROOF_WAT(:) * PDELT_ROOF(:) * XLVTT * PRHOA(:)
!
!-------------------------------------------------------------------------------
!
!*      9.     Fluxes at snow-free roads
!              -------------------------
!
!                                            net radiation
!
PRN_ROAD(:) = PABS_SW_ROAD(:) + PABS_LW_ROAD(:)
!
!                                            sensible heat flux
!
PH_ROAD(:) =  (PTS_ROAD(:) - PT_CANYON(:)) * PAC_ROAD(:)  * XCPD * PRHOA(:)
!
!                                            latent heat of evaporation from
!                                            the ground
!
PLE_ROAD(:) = ( PQSAT_ROAD(:) - PQ_CANYON(:) )                   &
              *   PAC_ROAD_WAT(:) * PDELT_ROAD(:) * XLVTT * PRHOA(:)
!
!
!-------------------------------------------------------------------------------
!
!*     10.     Fluxes at walls
!              ---------------
!
!                                            net radiation
!
PRN_WALL(:) = PABS_SW_WALL(:) + PABS_LW_WALL(:)
!
!                                            sensible heat flux
!
PH_WALL(:) =  (PTS_WALL(:) - PT_CANYON(:)) * PAC_WALL(:) * XCPD * PRHOA(:)
!
!                                            latent heat of evaporation from
!                                            the ground
!
PLE_WALL(:) = 0.
!
!                                            heat flux into the ground
!
PGFLUX_WALL(:) = PRN_WALL(:) - PH_WALL(:) - PLE_WALL(:)
!
!-------------------------------------------------------------------------------
!
!*     12.     Snow-free and snow-covered surfaces averaging
!              ---------------------------------------------
!
!*     12.1    roads
!              -----
!
!                                            heat flux into the ground
!
!
PGFLUX_ROAD (:) =  PDF_ROAD(:) * (PRN_ROAD(:) - PH_ROAD (:) - PLE_ROAD(:) ) &
                 + PDN_ROAD(:) *  PGSNOW_ROAD(:)
!
!
!                                            net radiation
!
PRN_ROAD(:) = PRN_ROAD(:) * PDF_ROAD(:) + PRNSNOW_ROAD(:) * PDN_ROAD(:)
!
!                                            sensible heat flux
!
PH_ROAD (:) = PH_ROAD (:) * PDF_ROAD(:) + PHSNOW_ROAD(:)  * PDN_ROAD(:)
!
!                                            latent heat of evaporation from
!                                            the ground
!
PLE_ROAD (:) = PLE_ROAD    (:) * PDF_ROAD(:)
PLES_ROAD(:) = PLESNOW_ROAD(:) * PDN_ROAD(:)
!
!*     12.2    roofs
!              -----
!
!                                            heat flux into the ground
!
!
PGFLUX_ROOF (:) =  PDF_ROOF(:) * (PRN_ROOF(:) - PH_ROOF (:) - PLE_ROOF(:) ) &
                 + PDN_ROOF(:) *  PGSNOW_ROOF(:)
!
!
!                                            net radiation
!
PRN_ROOF(:) = PRN_ROOF(:) * PDF_ROOF(:) + PRNSNOW_ROOF(:) * PDN_ROOF(:)
!
!                                            sensible heat flux
!
PH_ROOF (:) = PH_ROOF (:) * PDF_ROOF(:) + PHSNOW_ROOF(:)  * PDN_ROOF(:)
!
!                                            latent heat of evaporation from
!                                            the ground
!
PLE_ROOF (:) = PLE_ROOF    (:) * PDF_ROOF(:) 
PLES_ROOF(:) = PLESNOW_ROOF(:) * PDN_ROOF(:)
!
!-------------------------------------------------------------------------------
!
!*     13.     Momentum fluxes
!              ---------------
!
!
PUSTAR_TOWN(:) = SQRT( PCD(:) * PVMOD(:)**2 )
!
!-------------------------------------------------------------------------------
!
!*     14.     Averaged fluxes
!              ---------------
!
PRN_TOWN(:)= PROOF(:) * PRN_ROOF    (:) * PTOTS_O_HORS(:)  &
           + PROAD(:) * PRN_ROAD    (:) * PTOTS_O_HORS(:)  &
           + PWALL(:) * PRN_WALL    (:) * PTOTS_O_HORS(:)
!
PH_TOWN (:)= PROOF(:) * PH_ROOF    (:) * PTOTS_O_HORS(:)  &
           + PROAD(:) * PH_ROAD    (:) * PTOTS_O_HORS(:)  &
           + PWALL(:) * PH_WALL    (:) * PTOTS_O_HORS(:)  &
           + PH_TRAFFIC(:)                                &
           + PH_INDUSTRY(:)
!
PLE_TOWN(:)= PROOF(:) * (PLE_ROOF(:)+PLESNOW_ROOF(:)) * PTOTS_O_HORS(:) &
           + PROAD(:) * (PLE_ROAD(:)+PLESNOW_ROAD(:)) * PTOTS_O_HORS(:) &
           + PWALL(:) *  PLE_WALL(:)                  * PTOTS_O_HORS(:) &
           + PLE_TRAFFIC(:)                                             &
           + PLE_INDUSTRY(:)
!
PGFLUX_TOWN(:)= PROOF(:) * PGFLUX_ROOF(:) * PTOTS_O_HORS(:) &
              + PROAD(:) * PGFLUX_ROAD(:) * PTOTS_O_HORS(:) &
              + PWALL(:) * PGFLUX_WALL(:) * PTOTS_O_HORS(:)
!
!-------------------------------------------------------------------------------
!
!*     15.     Infra-red Radiative properties
!              ------------------------------
!
ZLW_UP(:) = PLW_RAD(:)                                  &
          - ( PROOF(:)*PDF_ROOF(:)*PABS_LW_ROOF     (:) &
             +PROOF(:)*PDN_ROOF(:)*PABS_LW_SNOW_ROOF(:) &
             +PROAD(:)*PDF_ROAD(:)*PABS_LW_ROAD     (:) &
             +PROAD(:)*PDN_ROAD(:)*PABS_LW_SNOW_ROAD(:) &
             +PWALL(:)            *PABS_LW_WALL     (:) &
            )*PTOTS_O_HORS(:)
!
!
!* town emissivity
!
PEMIS_TOWN(:) =   PBLD(:)     *PDF_ROOF(:)*PEMIS_ROOF(:)                            &
                + PBLD(:)     *PDN_ROOF(:)*PESNOW_ROOF(:)                           &
                + (1.-PBLD(:))*PDF_ROAD(:)                           *PLW_S_TO_R(:) &
                + (1.-PBLD(:))*PDN_ROAD(:)                           *PLW_S_TO_N(:) &
                +              PWALL(:)*PTOTS_O_HORS(:)              *PLW_S_TO_W(:)
!
!
!* town radiative surface temperature
!
PTS_TOWN(:)   = ((ZLW_UP(:) - PLW_RAD(:)*(1.-PEMIS_TOWN(:))) /PEMIS_TOWN(:)/XSTEFAN)**0.25
!
!-------------------------------------------------------------------------------
!
!*     16.     Averaged evaporative flux (kg/m2/s)
!              -----------------------------------
PEVAP_TOWN(:) = PROOF(:) * ((PLE_ROOF(:)/XLVTT)+(PLESNOW_ROOF(:)/XLSTT)) * PTOTS_O_HORS(:) &
              + PROAD(:) * ((PLE_ROAD(:)/XLVTT)+(PLESNOW_ROAD(:)/XLSTT)) * PTOTS_O_HORS(:) &
              + PWALL(:) *  (PLE_WALL(:)/XLVTT)                          * PTOTS_O_HORS(:) &
              + (PLE_TRAFFIC(:) /XLVTT)                                                    &
              + (PLE_INDUSTRY(:)/XLVTT)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_FLUXES
