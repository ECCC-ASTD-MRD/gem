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
    MODULE MODI_URBAN_DRAG2
!
!
!
INTERFACE
!
!
    SUBROUTINE URBAN_DRAG2(PTSTEP, PT_CANYON, PQ_CANYON,                    &
                          PTS_ROOF, PTS_ROAD, PTS_WALL, PDELT_SNOW_ROOF,    &
                          PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,               &
                          PZREF, PUREF, PVMOD, PVDIR, PLAT,                 &
                          PZ0_TOWN, PZ0_ROOF, PZ0_ROAD,                     &
                          PBLD, PBLD_HEIGHT, PCAN_HW_RATIO,                 &
                          PWALL_O_ROAD,                                     &
                          PWS_ROOF, PWS_ROAD,                               &
                          PWS_ROOF_MAX, PWS_ROAD_MAX,                       &
                          PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,   &
                          PCD, PAC_ROOF, PAC_ROOF_WAT,                      &
                          PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,        &
                          PU_CAN, PRI,                                      &
                          ZTRDZT,ZTRFZT,ZUZT,ZVZT                           )
!
!*      0.1    declarations of arguments 
!
!
REAL,               INTENT(IN)    :: PTSTEP         ! time-step
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON      ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON      ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL       ! surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_SNOW_ROOF! fraction of snow
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! specific humidity
                                                    ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD          ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PVDIR          ! direction of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PLAT           ! latitude
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
                                                    ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PZREF          ! reference height of the first
                                                    ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF          ! reference height of the first
                                                    ! atmospheric level (wind)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN       ! roughness length for momentum

REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ROOF       ! roughness length for momentum for roof
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ROAD       ! roughness length for momentum for road


REAL, DIMENSION(:), INTENT(IN)    :: PBLD           ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT    ! h
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO  ! h/W
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_ROAD   ! wall surf. / road surf.
!
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF       ! roof water content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD       ! road water content (kg/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROOF_MAX   ! maximum deepness of roof
REAL, DIMENSION(:), INTENT(IN)    :: PWS_ROAD_MAX   ! and water reservoirs (kg/m2)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD     ! qsat(Ts)
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF     ! water fraction on
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD     ! snow-free surfaces
REAL, DIMENSION(:), INTENT(OUT)   :: PCD            ! drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL       ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! walls 
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD       ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! roads
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT   ! aerodynamical conductance
!                                                   ! between canyon air and
!                                                   ! road (for water)
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP        ! aerodynamical conductance
!                                                   ! between canyon top and atm.
REAL, DIMENSION(:), INTENT(OUT)   :: PU_CAN         ! hor. wind in canyon
REAL, DIMENSION(:), INTENT(OUT)   :: PRI            ! Town Richardson number
!
!
REAL, DIMENSION(:), INTENT(OUT)   :: ztrfzt ! output for diasurf T,q at z=zt over roof
REAL, DIMENSION(:), INTENT(OUT)   :: ztrdzt ! output for diasurf T,q at z=zt over roof
REAL, DIMENSION(:), INTENT(OUT)   :: zuzt  ! output for diasurf u,v at z=zt over roof
REAL, DIMENSION(:), INTENT(OUT)   :: zvzt  ! output for diasurf u,v at z=zt over roof
!
END SUBROUTINE URBAN_DRAG2
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_URBAN_DRAG2
