!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END --------------------------

!#TODO: transform the s/r into a module instead of making an interface module to avoid non matching interface

module MODI_URBAN_THERMAL_STRESS

interface

subroutine URBAN_THERMAL_STRESS(PT_CAN, PQ_CAN, PTI_BLD, PQI_BLD,        &
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
                    PQ8_H,PQ9_H,PQ10_H,PQ11_H,PQ12_H,PQ13_H)
!
!*      0.1    declarations of arguments
real, dimension(:), intent(IN)  :: PBLD ! Building surface fraction
real, dimension(:), intent(IN)  :: PBLD_HEIGHT ! Building surface fraction
real, dimension(:), intent(IN)  :: PWALL_O_HOR ! Building surface fraction

real, dimension(:), intent(IN)  :: PT_CAN  ! Air canyon temperature (K)
real, dimension(:), intent(IN)  :: PTA     ! Air temperature above the roof (K)
real, dimension(:), intent(IN)  :: PTI_BLD ! Indoor air temperature (K)

real, dimension(:), intent(IN)  :: PQ_CAN  ! Canyon specific humidity (kg/kg)
real, dimension(:), intent(IN)  :: PQA     ! Air specific humidity over the roof (kg/kg)
real, dimension(:), intent(IN)  :: PQI_BLD ! Indoor specific humidity (kg/kg)

real, dimension(:), intent(IN)  :: PPA
real, dimension(:), intent(IN)  :: PPS

real, dimension(:), intent(IN)  :: PU_CAN    !  Air canyon temperature (K)
real, dimension(:), intent(IN)  :: PU10      !  Canyon wind speed at 10m (m/s)
real, dimension(:), intent(IN)  :: PURF      !  Air wind speed over the roof (m/s)
real, dimension(:), intent(IN)  :: PU10RF    !  Air wind speed 10-m over the roof (m/s)

real, dimension(:), intent(IN)  :: PZENITH   ! solar zenithal angle (rad from vert.)

real, dimension(:), intent(IN)  :: PDIR_SW      !Direct solar radiation (W/m²)
real, dimension(:), intent(IN)  :: PSCA_SW      !Diffuse solar radiation (W/m²)
real, dimension(:), intent(IN)  :: PREF_SW_GRND !Solar radiation reflected by ground [road + garden] (W/m²)
real, dimension(:), intent(IN)  :: PREF_SW_FAC  !Solar radiation reflected by facade [wall + glazing] (W/m²)
real, dimension(:), intent(IN)  :: PREF_SW_ROOF !Solar radiation reflected by the roof (W/m²)

real, dimension(:), intent(IN)  :: PLW_RAD       !Atmospheric longwave radiation (W/m²)
real, dimension(:), intent(IN)  :: PEMIT_LW_FAC  !Longwave radiation emitted by the facade [wall + glazing] (W/m²)
real, dimension(:), intent(IN)  :: PEMIT_LW_GRND !Longwave radiation emitted by the ground [road + garden] (W/m²)
real, dimension(:), intent(IN)  :: PEMIT_LW_ROOF !Longwave radiation emitted by the roof (W/m²)

real, dimension(:), intent(IN)  :: PTRAD_IN      ! body MRT inside building (K)
real, dimension(:), intent(OUT) :: PUTCI_IN      !UTCI for indoor person (°C)

real, dimension(:), intent(OUT) :: PUTCI_OUTSUN   !UTCI for outdoor person at sun (°C)
real, dimension(:), intent(OUT) :: PUTCI_OUTSHADE !UTCI for outdoor person in shade (°C)
real, dimension(:), intent(OUT) :: PUTCI_RFSUN    !UTCI for outdoor person on the roof at sun (°C)
real, dimension(:), intent(OUT) :: PUTCI_RFSHADE  !UTCI for outdoor person on the roof in shade (°C)

real, dimension(:), intent(OUT) :: WBGT_SUN         ! WBGT  wet bulb globe temperature in the street (C)
real, dimension(:), intent(OUT) :: WBGT_SHADE       ! WBGT  wet bulb globe temperature in the street (C)
real, dimension(:), intent(OUT) :: WBGT_RFSUN       ! WBGT  wet bulb globe temperature on the roof (C)
real, dimension(:), intent(OUT) :: WBGT_RFSHADE     ! WBGT  wet bulb globe temperature on the roof (C)

! optional to get out
real, dimension(:), intent(OUT) :: PTRAD_HSUN
real, dimension(:), intent(OUT) :: PTRAD_HSHADE
real, dimension(:), intent(OUT) :: PTRAD_HRFSUN
real, dimension(:), intent(OUT) :: PTRAD_HRFSHADE
real, dimension(:), intent(OUT) :: PTRAD_GSUN
real, dimension(:), intent(OUT) :: PTRAD_GSHADE
real, dimension(:), intent(OUT) :: PTRAD_GRFSUN
real, dimension(:), intent(OUT) :: PTRAD_GRFSHADE

 real, dimension(:), intent(OUT) :: PTGLOBE_SUN      ! Globe Temperature in the exposed street (K)
 real, dimension(:), intent(OUT) :: PTGLOBE_SHADE    ! Globe Temperature in the shaded street (K)
 real, dimension(:), intent(OUT) :: PTGLOBE_RFSUN    ! Globe Temperature  on the exposed roof (K)
 real, dimension(:), intent(OUT) :: PTGLOBE_RFSHADE  ! Globe Temperature  on the shaded roof  (K)

 real, dimension(:), intent(OUT) :: PTWETB          ! wet-bulb temperature in the street
 real, dimension(:), intent(OUT) :: PTWETB_ROOF     ! wet-bulb temperature on the roof

 real, dimension(:), intent(OUT)   :: PQ1_H  ! energy components for the standing standard clothed human
 real, dimension(:), intent(OUT)   :: PQ2_H
 real, dimension(:), intent(OUT)   :: PQ3_H
 real, dimension(:), intent(OUT)   :: PQ4_H
 real, dimension(:), intent(OUT)   :: PQ5_H
 real, dimension(:), intent(OUT)   :: PQ6_H
 real, dimension(:), intent(OUT)   :: PQ7_H
 real, dimension(:), intent(OUT)   :: PQ8_H
 real, dimension(:), intent(OUT)   :: PQ9_H
 real, dimension(:), intent(OUT)   :: PQ10_H
 real, dimension(:), intent(OUT)   :: PQ11_H
 real, dimension(:), intent(OUT)   :: PQ12_H
 real, dimension(:), intent(OUT)   :: PQ13_H

end subroutine URBAN_THERMAL_STRESS

end interface

end module MODI_URBAN_THERMAL_STRESS
