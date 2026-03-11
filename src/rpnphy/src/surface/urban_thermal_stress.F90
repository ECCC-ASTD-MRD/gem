
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
!-------------------------------------- LICENCE END ---------------------------

!/@*
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

!    PURPOSE       : Computes thermal stress indicators in the street and over rooftop
!    AUTHOR        :  S. Leroyer   (Original  10/2016),  based on CNRM/G. Pigeon  UTCI code
!    REFERENCE     :  Leroyer et al. (2018), urban climate
!    MODIFICATIONS :  S. Leroyer (2020): wetbulb in C instead of K
!    MODIFICATIONS :  S. Leroyer (2026): Modularity for computation
!    METHOD        :  For each case,
! (1) COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY
!           -- for UTCI  (standard clothed standing human) ZOPT_BODY=1
!           -- for WBGT  (globe sensor)  ZOPT_BODY=2
! (3) COMPUTE THE MEAN RADIANT TEMPERATURES
!           -- reference MRT and var. used for UTCI  comp. (standard clothed standing human)
!                * MRT that will be in the output as it  corresponds to standards
!           -- for WBGT  (globe sensor)
! (4) COMPUTE THE GLOBE TEMPERATURE - for a black globe sensor
! (5) compute the (psychometric) wet-bulb temperatures
! (6) compute the wet bulb globe temperature indices  (WBGT)
! (7) compute the Universal Thermal and Climate Index UTCI

!-------------------------------------------------------------------------------

!*       0.     DECLARATIONS
!               ------------

use sfc_options,   only : thermal_stress_roof, thermal_stress_utci, thermal_stress_shade
use MODD_CSTS, only : XTT
use MODI_UTCI_APPROX
use MODI_MRT_BODY
use MODI_TGLOBE_BODY
use MODI_URBAN_OUTQENV
use MODI_WETBULBT

! USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
! USE PARKIND1  ,ONLY : JPRB

implicit none
!!!#include <arch_specific.hf>

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
real, dimension(:), intent(OUT) :: PTRAD_HSUN        ! for a human body
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

 real, dimension(:), intent(OUT) :: PTWETB          ! wet-bulb temperature in the street (C)
 real, dimension(:), intent(OUT) :: PTWETB_ROOF     ! wet-bulb temperature on the roof   (C)

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

!  declarations of local variables
real, dimension(size(PTI_BLD)) :: ZEHPA !water vapour pressure (hPa)
real, dimension(size(PTI_BLD)) :: ZUIN !indoor air wind speed (m/s)
real, dimension(size(PTI_BLD)) :: ZUNDEF
! energy components for the globe
real, dimension(size(PTI_BLD)) :: PQ1_G
real, dimension(size(PTI_BLD)) :: PQ2_G
real, dimension(size(PTI_BLD)) :: PQ3_G
real, dimension(size(PTI_BLD)) :: PQ4_G
real, dimension(size(PTI_BLD)) :: PQ5_G
real, dimension(size(PTI_BLD)) :: PQ6_G
real, dimension(size(PTI_BLD)) :: PQ7_G
real, dimension(size(PTI_BLD)) :: PQ8_G
real, dimension(size(PTI_BLD)) :: PQ9_G
real, dimension(size(PTI_BLD)) :: PQ10_G
real, dimension(size(PTI_BLD)) :: PQ11_G
real, dimension(size(PTI_BLD)) :: PQ12_G
real, dimension(size(PTI_BLD)) :: PQ13_G

real :: ZEB_G = 0.957 !emissivity of a globe sensor          WBGT
real :: ZEB_H = 0.97  !emissivity of clothed human body      UTCI
real :: ZAB_H = 0.7   !absorption coef of solar radiation by human body        UTCI
real :: ZAB_G = 0.957 !absorption coef of solar radiation by globe sensor      WBGT
real :: ZHB_H = 1.7  !average height of human person (m)                       UTCI & WBGT
real :: ZHB_G=  1.7  !2.5 panam  ! average height of the globe sensor
real :: ZGD = 0.148  ! black globe sensor diameter in m (value given by Matt Wright for PanAm2015)
integer :: ZOPT
integer :: ZOPT_BODY

! REAL(KIND=JPRB) :: ZHOOK_HANDLE

! IF (LHOOK) CALL DR_HOOK('UTCI_TEB',0,ZHOOK_HANDLE)
ZUNDEF(:)=0.0

 ZOPT=2   ! street
!========================================================
! COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY (standard clothed standing human)
!========================================================
 ZOPT_BODY=1

     call URBAN_OUTQENV(PSCA_SW, PREF_SW_FAC, PREF_SW_GRND, ZUNDEF,  &
               PEMIT_LW_FAC, PEMIT_LW_GRND, ZUNDEF, PLW_RAD,         &
               PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,     &
               PQ1_H,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H,PQ7_H,            &
               PQ8_H,PQ9_H,PQ10_H,PQ11_H,PQ12_H,PQ13_H,              &
               ZOPT,ZOPT_BODY,ZEB_H,ZAB_H, ZHB_H  )

!========================================================
! COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY (globe sensor)
!========================================================
 ZOPT_BODY=2

     call URBAN_OUTQENV(PSCA_SW, PREF_SW_FAC, PREF_SW_GRND, ZUNDEF,  &
               PEMIT_LW_FAC, PEMIT_LW_GRND, ZUNDEF, PLW_RAD,         &
               PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,     &
               PQ1_G,PQ2_G,PQ3_G,PQ4_G,PQ5_G,PQ6_G,PQ7_G,            &
               PQ8_G,PQ9_G,PQ10_G,PQ11_G,PQ12_G,PQ13_G,              &
               ZOPT,ZOPT_BODY, ZEB_G, ZAB_G, ZHB_G )

!========================================================
! COMPUTE THE MEAN RADIANT TEMPERATURES (MRT)
!========================================================

! Standards and MRT output
PTRAD_HSUN     = MRT_BODY(ZEB_H,PQ1_H,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H,PQ7_H)

! MRT equivalent to the black globe sensor and used for WBGT comp.
PTRAD_GSUN    =MRT_BODY(ZEB_G,PQ1_G,PQ2_G,PQ3_G,PQ4_G,PQ5_G,PQ6_G,PQ7_G)

!========================================================
! COMPUTE THE GLOBE TEMPERATURE
!========================================================

PTGLOBE_SUN     = TGLOBE_BODY(PTRAD_GSUN, PT_CAN, PU_CAN, ZGD, ZEB_G)

!========================================================
! compute the (psychometric) wet-bulb temperatures
!========================================================

PTWETB       = WETBULBT(PPS, PT_CAN -XTT, PQ_CAN)

!========================================================
! compute the wet bulb globe temperature indices  (WBGT)
!========================================================

WBGT_SUN   = 0.2 * (PTGLOBE_SUN -XTT)     +   &
             0.7 *  PTWETB            +   &
             0.1 * (PT_CAN-XTT)

!========================================================
! COMPUTE VARS. IN THE SHADE
!========================================================
    IF_THERMAL_STRESS_SHADE: if ( thermal_stress_shade ) then

! MRT ref and used for UTCI
PTRAD_HSHADE   = MRT_BODY(ZEB_H,ZUNDEF,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H, &
                                                                    PQ7_H)
! MRT for WBGT
PTRAD_GSHADE  =MRT_BODY(ZEB_G,ZUNDEF,PQ2_G,PQ3_G,PQ4_G,PQ5_G,PQ6_G,  &
                          PQ7_G)

! Globe temperature
PTGLOBE_SHADE   = TGLOBE_BODY(PTRAD_GSHADE, PT_CAN, PU_CAN, ZGD, ZEB_G)

! WBGT - formula in the shade
WBGT_SHADE = 0.3 * (PTGLOBE_SHADE -XTT)   +   &
             0.7 *  PTWETB

    endif IF_THERMAL_STRESS_SHADE

!========================================================
! compute the Universal Thermal and Climate Index UTCI
!========================================================
    IF_THERMAL_STRESS_UTCI: if ( thermal_stress_utci ) then

! Inside buildings - UTCI_IN
ZEHPA = PQI_BLD * PPS /(0.622 + 0.378 * PQI_BLD) / 100.
ZUIN = 0.5
PUTCI_IN = UTCI_APPROX(PTI_BLD - XTT, ZEHPA, PTRAD_IN - XTT, ZUIN)

! 2-calculation of UTCI
 ZEHPA = PQ_CAN * PPS / (0.622 + 0.378 * PQ_CAN) /100.

PUTCI_OUTSUN = UTCI_APPROX(PT_CAN - XTT, ZEHPA, PTRAD_HSUN - XTT, PU10)

! SAME IN SHADE
    IF_THERMAL_STRESS_UTCI_SHADE: if ( thermal_stress_shade ) then

PUTCI_OUTSHADE = UTCI_APPROX(PT_CAN - XTT, ZEHPA, PTRAD_HSHADE - XTT, PU10)

    endif IF_THERMAL_STRESS_UTCI_SHADE

! INSIDE BUILDING - UTCI_IN
!    IF_THERMAL_STRESS_UTCI_IN: if ( thermal_stress_in ) then

ZEHPA = PQI_BLD * PPS /(0.622 + 0.378 * PQI_BLD) / 100.
ZUIN = 0.5
PUTCI_IN = UTCI_APPROX(PTI_BLD - XTT, ZEHPA, PTRAD_IN - XTT, ZUIN)

!    endif IF_THERMAL_STRESS_UTCI_IN

    endif IF_THERMAL_STRESS_UTCI


!========================================================
! COMPUTE VAR. OVER THE ROOF
!========================================================
    IF_THERMAL_STRESS_ROOF: if ( thermal_stress_roof ) then

 ZOPT=3   ! roof

! COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY (standard clothed standing human)
 ZOPT_BODY=1

     call URBAN_OUTQENV(PSCA_SW, PREF_SW_FAC, ZUNDEF, PREF_SW_ROOF,  &
               PEMIT_LW_FAC, ZUNDEF, PEMIT_LW_ROOF, PLW_RAD,         &
               PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,     &
               PQ1_H,PQ2_H,PQ3_H,PQ4_H,PQ5_H,PQ6_H,PQ7_H,            &
               PQ8_H,PQ9_H,PQ10_H,PQ11_H,PQ12_H,PQ13_H,              &
               ZOPT,ZOPT_BODY,ZEB_H,ZAB_H, ZHB_H  )

! COMPUTE THE ENERGY BUDGETS RECEIVED BY A BODY (globe sensor)
 ZOPT_BODY=2

     call URBAN_OUTQENV(PSCA_SW, PREF_SW_FAC, ZUNDEF, PREF_SW_ROOF,  &
               PEMIT_LW_FAC, ZUNDEF, PEMIT_LW_ROOF, PLW_RAD,         &
               PBLD, PBLD_HEIGHT, PWALL_O_HOR, PDIR_SW, PZENITH,     &
               PQ1_G,PQ2_G,PQ3_G,PQ4_G,PQ5_G,PQ6_G,PQ7_G,            &
               PQ8_G,PQ9_G,PQ10_G,PQ11_G,PQ12_G,PQ13_G,              &
               ZOPT,ZOPT_BODY, ZEB_G, ZAB_G, ZHB_G )

!  MRT - REF and for UTCI
PTRAD_HRFSUN   = MRT_BODY(ZEB_H,PQ1_H,PQ8_H,PQ9_H,PQ10_H,PQ11_H,      &
                          PQ12_H,PQ13_H)

!  MRT - for WBGT
PTRAD_GRFSUN  =MRT_BODY(ZEB_G,PQ1_G,PQ8_G,PQ9_G,PQ10_G,PQ11_G,      &
                          PQ12_G,PQ13_G)

! GLOBE T.
PTGLOBE_RFSUN   = TGLOBE_BODY(PTRAD_GRFSUN, PTA, PURF, ZGD, ZEB_G)

! T wet-bulb
PTWETB_ROOF  = WETBULBT(PPA, PTA -XTT, PQA)

! WBGT
WBGT_RFSUN = 0.2 * (PTGLOBE_RFSUN -XTT)   +   &
             0.7 *  PTWETB_ROOF       +   &
             0.1 * (PTA-XTT)

    IF_THERMAL_STRESS_ROOF_SHADE: if ( thermal_stress_shade ) then

! MRT - REF and for UTCI
PTRAD_HRFSHADE = MRT_BODY(ZEB_H,ZUNDEF,PQ8_H,PQ9_H,PQ10_H,PQ11_H, &
                          PQ12_H,PQ13_H)

! MRT - for a black globe sensor
PTRAD_GRFSHADE=MRT_BODY(ZEB_G,ZUNDEF,PQ8_G,PQ9_G,PQ10_G,PQ11_G,  &
                          PQ12_G,PQ13_G)

! COMPUTE THE GLOBE TEMPERATURE
PTGLOBE_RFSHADE = TGLOBE_BODY(PTRAD_GRFSHADE, PTA, PURF, ZGD, ZEB_G)

! (WBGT) formula in the shade
WBGT_RFSHADE = 0.3 * (PTGLOBE_RFSHADE -XTT) +   &
               0.7 *  PTWETB_ROOF

    endif IF_THERMAL_STRESS_ROOF_SHADE

!========================================================
! compute the Universal Thermal and Climate Index UTCI over ROOF
!========================================================
    IF_THERMAL_STRESS_UTCI_ROOF: if ( thermal_stress_utci ) then

! 4-calculation of UTCI_RFSUN
ZEHPA = PQA * PPA / (0.622 + 0.378 * PQA) /100.

PUTCI_RFSUN = UTCI_APPROX(PTA- XTT, ZEHPA, PTRAD_HRFSUN - XTT, PU10RF)

    IF_THERMAL_STRESS_ROOF_SHADE_UTCI: if ( thermal_stress_shade ) then
PUTCI_RFSHADE =UTCI_APPROX(PTA-XTT,ZEHPA,PTRAD_HRFSHADE- XTT, PU10RF)
    endif IF_THERMAL_STRESS_ROOF_SHADE_UTCI

    endif IF_THERMAL_STRESS_UTCI_ROOF

    endif IF_THERMAL_STRESS_ROOF

end subroutine URBAN_THERMAL_STRESS
