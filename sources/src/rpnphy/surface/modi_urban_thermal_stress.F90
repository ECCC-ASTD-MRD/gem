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

!#TODO: transform the s/r into a module instead of making an interface module to avoid non matching interface

    MODULE MODI_URBAN_THERMAL_STRESS
!
!
!
INTERFACE
!
!
SUBROUTINE URBAN_THERMAL_STRESS(PT_CAN, PQ_CAN, PTI_BLD, PQI_BLD,        &
                    PU_CAN, PU10, PURF, PU10RF, PPS,PPA,                 &
                    PREF_SW_GRND, PREF_SW_FAC, PSCA_SW, PDIR_SW, PZENITH,&
                    PEMIT_LW_FAC, PEMIT_LW_GRND, PLW_RAD, PTRAD_IN,      &
	                 PTA, PQA, PREF_SW_ROOF, PEMIT_LW_ROOF,               &
                    PBLD, PBLD_HEIGHT, PWALL_O_HOR,                      &
                    PUTCI_IN, PUTCI_OUTSUN, PUTCI_OUTSHADE,              &
                    PUTCI_RFSUN, PUTCI_RFSHADE,                          &
                    WBGT_SUN,WBGT_SHADE,WBGT_RFSUN,WBGT_RFSHADE,         &
                    PTRAD_HSUN, PTRAD_HSHADE, PTRAD_HRFSUN,              &
                    PTRAD_HRFSHADE,                                      &
                    PTRAD_GSUN, PTRAD_GSHADE, PTRAD_GRFSUN,              &
                    PTRAD_GRFSHADE,                                      &
                    PTGLOBE_SUN, PTGLOBE_SHADE, PTGLOBE_RFSUN,           &
                    PTGLOBE_RFSHADE, PTWETB, PTWETB_ROOF,                &
                    PQ1_H,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H,PQ7_H,           &
                    PQ8_H,PQ9_H,PQ10_H,PQ11_H,PQ12_H,PQ13_H        )
!
!*      0.1    declarations of arguments
REAL, DIMENSION(:), INTENT(IN)  :: PBLD ! Building surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PBLD_HEIGHT ! Building surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PWALL_O_HOR ! Building surface fraction

REAL, DIMENSION(:), INTENT(IN)  :: PT_CAN  ! Air canyon temperature (K) 
REAL, DIMENSION(:), INTENT(IN)  :: PTA     ! Air temperature above the roof (K)  
REAL, DIMENSION(:), INTENT(IN)  :: PTI_BLD ! Indoor air temperature (K) 

REAL, DIMENSION(:), INTENT(IN)  :: PQ_CAN  ! Canyon specific humidity (kg/kg)
REAL, DIMENSION(:), INTENT(IN)  :: PQA     ! Air specific humidity over the roof (kg/kg)
REAL, DIMENSION(:), INTENT(IN)  :: PQI_BLD ! Indoor specific humidity (kg/kg) 

REAL, DIMENSION(:), INTENT(IN)  :: PPA     
REAL, DIMENSION(:), INTENT(IN)  :: PPS      

REAL, DIMENSION(:), INTENT(IN)  :: PU_CAN    !  Air canyon temperature (K) 
REAL, DIMENSION(:), INTENT(IN)  :: PU10      !  Canyon wind speed at 10m (m/s)
REAL, DIMENSION(:), INTENT(IN)  :: PURF      !  Air wind speed over the roof (m/s)
REAL, DIMENSION(:), INTENT(IN)  :: PU10RF    !  Air wind speed 10-m over the roof (m/s)

REAL, DIMENSION(:), INTENT(IN)  :: PZENITH   ! solar zenithal angle (rad from vert.)

REAL, DIMENSION(:), INTENT(IN)  :: PDIR_SW      !Direct solar radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PSCA_SW      !Diffuse solar radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_GRND !Solar radiation reflected by ground [road + garden] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_FAC  !Solar radiation reflected by facade [wall + glazing] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PREF_SW_ROOF !Solar radiation reflected by the roof (W/m²)

REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD       !Atmospheric longwave radiation (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_FAC  !Longwave radiation emitted by the facade [wall + glazing] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_GRND !Longwave radiation emitted by the ground [road + garden] (W/m²)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF !Longwave radiation emitted by the roof (W/m²)

REAL, DIMENSION(:), INTENT(IN)  :: PTRAD_IN      ! body MRT inside building (K)
REAL, DIMENSION(:), INTENT(OUT) :: PUTCI_IN      !UTCI for indoor person (°C)

REAL, DIMENSION(:), INTENT(OUT) :: PUTCI_OUTSUN   !UTCI for outdoor person at sun (°C)
REAL, DIMENSION(:), INTENT(OUT) :: PUTCI_OUTSHADE !UTCI for outdoor person in shade (°C)
REAL, DIMENSION(:), INTENT(OUT) :: PUTCI_RFSUN    !UTCI for outdoor person on the roof at sun (°C)
REAL, DIMENSION(:), INTENT(OUT) :: PUTCI_RFSHADE  !UTCI for outdoor person on the roof in shade (°C)
!
REAL, DIMENSION(:), INTENT(OUT) :: WBGT_SUN         ! WBGT  wet bulb globe temperature in the street (C)
REAL, DIMENSION(:), INTENT(OUT) :: WBGT_SHADE       ! WBGT  wet bulb globe temperature in the street (C)
REAL, DIMENSION(:), INTENT(OUT) :: WBGT_RFSUN       ! WBGT  wet bulb globe temperature on the roof (C)
REAL, DIMENSION(:), INTENT(OUT) :: WBGT_RFSHADE     ! WBGT  wet bulb globe temperature on the roof (C)

! optional to get out
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_HSUN 
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_HSHADE
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_HRFSUN
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_HRFSHADE
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_GSUN 
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_GSHADE
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_GRFSUN
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD_GRFSHADE

 REAL, DIMENSION(:), INTENT(OUT) :: PTGLOBE_SUN      ! Globe Temperature in the exposed street (K)
 REAL, DIMENSION(:), INTENT(OUT) :: PTGLOBE_SHADE    ! Globe Temperature in the shaded street (K)
 REAL, DIMENSION(:), INTENT(OUT) :: PTGLOBE_RFSUN    ! Globe Temperature  on the exposed roof (K)
 REAL, DIMENSION(:), INTENT(OUT) :: PTGLOBE_RFSHADE  ! Globe Temperature  on the shaded roof  (K)

 REAL, DIMENSION(:), INTENT(OUT) :: PTWETB          ! wet-bulb temperature in the street
 REAL, DIMENSION(:), INTENT(OUT) :: PTWETB_ROOF     ! wet-bulb temperature on the roof

 REAL, DIMENSION(:), INTENT(OUT)   :: PQ1_H  ! energy components for the standing standard clothed human
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ2_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ3_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ4_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ5_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ6_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ7_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ8_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ9_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ10_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ11_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ12_H
 REAL, DIMENSION(:), INTENT(OUT)   :: PQ13_H
!
!
!
END SUBROUTINE URBAN_THERMAL_STRESS
!
!
! 
END INTERFACE
!
!
!
END MODULE MODI_URBAN_THERMAL_STRESS
