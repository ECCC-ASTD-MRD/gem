!-------------------------------------- LICENCE BEGIN ------------------------
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

subroutine TEB2(PT_CANYON, PQ_CANYON, PU_CANYON,                         &
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
                     PZ0_TOWN,PZ0_ROOF,PZ0_ROAD,                              &
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
                     PTS_TOWN, PEMIS_TOWN, PDIR_ALB_TOWN, PSCA_ALB_TOWN,      &
                     PRESA_TOWN,                                              &
            PTRAD_IN, PTRAD_SUN, PTRAD_SHADE, PTRAD_RFSUN, PTRAD_RFSHADE,     &
            PTGLOBE_SUN, PTGLOBE_SHADE, PTGLOBE_RFSUN, PTGLOBE_RFSHADE,       &
                    PTWETB,PTWETB_ROOF,                                       &
                    PUTCI_IN, PUTCI_OUTSUN, PUTCI_OUTSHADE,                   &
                    PUTCI_RFSUN, PUTCI_RFSHADE,                               &
                    WBGT_SUN, WBGT_SHADE, WBGT_RFSUN, WBGT_RFSHADE,           &
                    PUTCIC_IN,PUTCIC_OUTSUN,PUTCIC_OUTSHADE,                  &
                    PUTCIC_RFSUN,PUTCIC_RFSHADE,                              &
                    PTRFZT,PTRDZT,PURDZU,                                      &
                    PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,PQ8,PQ9,PQ10,PQ11,PQ12,PQ13)


!!****  *TEB*
!
!!    PURPOSE
!!    -------

!     Computes the evoultion of prognostic variables and the fluxes
!     over artificial surfaces as towns, taking into account the canyon like
!     geometry of urbanized areas.


!!**  METHOD
!     ------

!     The prognostic variables are:
!       - the surface temperature for roofs, roads, and walls
!       - the water reservoir, whose maximum value is 10mm


!    1 : Warning about snow
!        ******************

!     Except for snow mantel evolution, all other computation with snow
!   variables must be performed with these variables at previous time-step,
!   and NOT new time-step. This insure coherence between snow fractions
!   (computed at the begining) and other snow characteristics (albedo, Ts).


!    2 : computation of input solar radiation on each surface
!        ****************************************************

!    direct fluxes:
!    -------------

!    dir_Rg_road (Wm-2) =   S * 2*theta0/pi
!                         - S *2/tan(zen) * h/W /pi * (1-cos(theta0))

!    dir_Rg_wall (Wm-2) =   S / tan(zen) /pi * (1-cos(theta0))
!                         + S * W/h * (1/2 -theta0/pi)

!   where zen      is the zenithal angle, from horizon
!         h/W      is the aspect ratio of the canyon
!         S        is the direct solar radiation flux on a horizontal surface

!         theta0 = arcsin(min(W/h * tan(zen),1))

!   The surfaces will keep (1-a) times these fluxes, and reflect the
!   remaining

!    scattered fluxes:
!    ----------------

!   sca_Rg_road = sca_Rg * SVF_road

!   sca_Rg_wall = sca_Rg * SVF_wall


!    solar flux and isotropic reflections :
!    ------------------------------------

!  after 0 reflection, the absorbed part of the flux is:

!      ARg_r(0) = (1-a_r) (sca_Rg_road + dir_Rg_road)

!      ARg_w(0) = (1-a_w) (sca_Rg_wall + dir_Rg_wall)

!    and the reflected parts are

!      RRg_r(0) = a_r (sca_Rg_road + dir_Rg_road)

!      RRg_w(0) = a_w (sca_Rg_wall + dir_Rg_wall)

!  after n reflection:

!      ARg_r(n) = ARg_r(n-1) + RRg_w(n-1) * (1-  SVF_r)(1-a_r)

!      ARg_w(n) = ARg_w(n-1) + RRg_r(n-1) *      SVF_w (1-a_w)
!                            + RRg_w(n-1) * (1-2*SVF_w)(1-a_w)

!      RRg_r(n) = (1- SVF_r) a_r RRg_w(n-1)

!      RRg_w(n) =     SVF_w  a_w RRg_r(n-1)
!                +(1-2SVF_w) a_w RRg_w(n-1)


!   i.e.
!                                               n-1
!      ARg_r(n) = ARg_r(0) + (1-  SVF_r)(1-a_r) SUM RRg_w(k)
!                                               k=0

!                                               n-1
!      ARg_w(n) = ARg_w(0) +      SVF_w (1-a_w) SUM RRg_r(k)
!                                               k=0
!                                               n-1
!                          + (1-2*SVF_w)(1-a_w) SUM RRg_w(k)
!                                               k=0

! with

!     n                             n-1
!    SUM RRg_r(k) = (1-  SVF_r) a_r SUM RRg_w(k)      +  RRg_r(0)
!    k=0                            k=0

!     n                             n-1
!    SUM RRg_w(k) =      SVF_w  a_w SUM RRg_r(k)
!    k=0                            k=0
!                                   n-1
!                  +(1-2*SVF_w) a_w SUM RRg_w(k)      +  RRg_w(0)
!                                   k=0


!   Then

!     n                                        n-1
!    SUM RRg_w(k) =  (1-2*SVF_w)       a_w     SUM RRg_w(k)
!    k=0                                       k=0
!                                              n-2
!                  + (1-  SVF_r) SVF_w a_w a_r SUM RRg_w(k)
!                                              k=0

!                  + RRg_w(0) + SVF_w a_w RRg_r(0)




!  solving this system, lead after an infinity of reflections/absorptions:

!    inf                      RRg_w(0) + SVF_w a_w RRg_r(0)
!    SUM RRg_w(k) = ----------------------------------------------------
!    k=0             1 - (1-2*SVF_w) a_w + (1-  SVF_r) SVF_w a_w a_r


!    inf            (1-  SVF_r) a_r ( a_w SVF_w RRg_r(0) + RRg_w(0) ) + RRg_r(0)
!    SUM RRg_r(k) = ------------------------------------------------------------
!    k=0             1 - (1-2*SVF_w) a_w + (1-  SVF_r) SVF_w a_w a_r


! ARg_r(n) and ARg_w(n) follow


!    3 : drag coefficient for momentum
!        *****************************


!    4 : aerodynamical resistance for heat transfers
!        *******************************************


!    5 : equation for evolution of Ts_roof
!        *********************************


!       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )

!       H  = rho Cp CH V ( Ts (t+dt) - Tas )

!       LE = rho Lv CH V ( qs (t+dt) - qas )

!      where the as subscript denotes atmospheric values at ground level
!      (and not at first half level)


!    6 : equations for evolution of Ts_road and Ts_wall simultaneously
!        *************************************************************



!   Rn_w = abs_Rg_w
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * Rat
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)

!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * Rat
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)

!  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )

!  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )

!  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )

!  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )

! with again
!                AC_can * Swall/Sroad * Twall + AC_can * Troad + AC_top * Ta + H_traffic/Cp/rho/Sroad
!   Ta_canyon = -------------------------------------------------------------------------------------
!                AC_can * Swall/Sroad         + AC_can         + AC_top


!                 AC_can * delt_road * Hu_road * qsat(Troad) + AC_top * qa + LE_traffic/Lv/rho/Sroad
!   qa_canyon = ------------------------------------------------------------------------------------
!                 AC_can * delt_road                        + AC_top




!    7 : computation of fluxes for each surface type
!        *******************************************


!    8 : averaging of the fluxes
!        ***********************

!   This is done on the total exchange surface (roof + wall + road),
!  which is bigger than the horizontal surface (roof+road), leading
!  to bigger fluxes.

!   The fluxes due to industrial activity are directly added into the
!  atmosphere


!    9 : road reservoir evolution
!        ************************

!   The roof reservoir runoff goes directly into the road reservoir.

!   Runoff occurs for road reservoir (too much water), as well as drainage
!   (evacuation system, typical time scale: 1 day)




!!    EXTERNAL
!!    --------
!
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!
!!    MODD_CST
!
!
!!    REFERENCE
!!    ---------
!
!
!!    AUTHOR
!!    ------
!
!!  V. Masson           * Meteo-France *
!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98
!!     21 / 10 / 2003   P. Tulet    output aerodynamical resistance

!-------------------------------------------------------------------------------

!*       0.     DECLARATIONS
!               ------------

use MODD_CSTS,     only : XTT, XSTEFAN
use MODD_SNOW_PAR, only : XEMISSN, XANSMAX,SWE_CRIT
use SFC_OPTIONS,   only : THERMAL_STRESS

use MODE_THERMOS
use MODE_SURF_SNOW_FRAC

use MODI_URBAN_SOLAR_ABS
use MODI_URBAN_DRAG2
use MODI_URBAN_SNOW_EVOL
use MODI_ROOF_LAYER_E_BUDGET
use MODI_ROAD_WALL_LAYER_E_BUDGET
use MODI_URBAN_FLUXES
use MODI_URBAN_HYDRO
use MODI_BLD_E_BUDGET
use MODI_URBAN_THERMAL_STRESS

implicit none
!!!#include <arch_specific.hf>

!*      0.1    declarations of arguments


real, dimension(:), intent(INOUT) :: PT_CANYON     ! canyon air temperature
real, dimension(:), intent(INOUT) :: PQ_CANYON     ! canyon air specific humidity
real, dimension(:), intent(OUT)   :: PU_CANYON     ! canyon hor. wind
real, dimension(:), intent(INOUT) :: PTI_BLD       ! inside building temperature
real, dimension(:,:), intent(INOUT) :: PT_ROOF     ! roof layers temperatures
real, dimension(:,:), intent(INOUT) :: PT_ROAD     ! road layers temperatures
real, dimension(:,:), intent(INOUT) :: PT_WALL     ! wall layers temperatures
real, dimension(:), intent(INOUT) :: PWS_ROOF      ! roof water reservoir
real, dimension(:), intent(INOUT) :: PWS_ROAD      ! road water reservoir
real, dimension(:),   intent(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
real, dimension(:),   intent(INOUT) :: PTSNOW_ROOF ! snow layers temperature
real, dimension(:),   intent(INOUT) :: PRSNOW_ROOF ! snow layers density
real, dimension(:),   intent(INOUT) :: PASNOW_ROOF ! snow albedo
real, dimension(:),   intent(INOUT) :: PESNOW_ROOF ! snow emissivity
real, dimension(:),   intent(INOUT) :: PTSSNOW_ROOF! snow surface temperature
real, dimension(:),   intent(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
real, dimension(:),   intent(INOUT) :: PTSNOW_ROAD ! snow layers temperature
real, dimension(:),   intent(INOUT) :: PRSNOW_ROAD ! snow layers density
real, dimension(:),   intent(INOUT) :: PASNOW_ROAD ! snow albedo
real, dimension(:),   intent(INOUT) :: PESNOW_ROAD ! snow emissivity
real, dimension(:),   intent(INOUT) :: PTSSNOW_ROAD! snow surface temperature
real, dimension(:), intent(IN)    :: PPS           ! pressure at the surface
real, dimension(:), intent(IN)    :: PPA           ! pressure at the first atmospheric level
real, dimension(:), intent(IN)    :: PEXNS         ! surface exner function
real, dimension(:), intent(IN)    :: PTA           ! temperature at the lowest level
real, dimension(:), intent(IN)    :: PQA           ! specific humidity
                                                   ! at the lowest level
real, dimension(:), intent(IN)    :: PVMOD         ! module of the horizontal wind
real, dimension(:), intent(IN)    :: PVDIR         ! direction of the horizontal wind
real, dimension(:), intent(IN)    :: PEXNA         ! exner function
                                                   ! at the lowest level
real, dimension(:), intent(IN)    :: PRHOA         ! air density
                                                   ! at the lowest level
real, dimension(:), intent(IN)    :: PLW_RAD       ! atmospheric infrared radiation
real, dimension(:), intent(IN)    :: PDIR_SW_RAD   ! incoming direct solar radiation
                                                   ! on an horizontal surface
real, dimension(:), intent(IN)    :: PZENITH       ! solar zenithal angle
real, dimension(:), intent(IN)    :: PSCA_SW_RAD   ! scattered incoming solar rad.
real, dimension(:), intent(IN)    :: PRR           ! rain rate
real, dimension(:), intent(IN)    :: PSR           ! snow rate
real, dimension(:), intent(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic
real, dimension(:), intent(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic
real, dimension(:), intent(IN)    :: PH_INDUSTRY   ! anthropogenic sensible
!                                                  ! heat fluxes due to factories
real, dimension(:), intent(IN)    :: PLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories
real, dimension(:), intent(IN)    :: PZREF         ! reference height of the first
                                                   ! atmospheric level (temperature)
real, dimension(:), intent(IN)    :: PUREF         ! reference height of the first
                                                   ! atmospheric level (wind)
real,               intent(IN)    :: PTSTEP        ! time step
real, dimension(:), intent(IN)    :: PZ0_TOWN      ! town roughness length
                                                   ! for momentum
real, dimension(:), intent(IN)    :: PZ0_ROOF      ! roof roughness length
                                                   ! for momentum
real, dimension(:), intent(IN)    :: PZ0_ROAD      ! road roughness length
                                                   ! for momentum
real, dimension(:), intent(IN)    :: PBLD          ! fraction of buildings
real, dimension(:), intent(IN)    :: PBLD_HEIGHT   ! buildings h
real, dimension(:), intent(IN)    :: PWALL_O_HOR   ! wall surf. / hor. surf.
real, dimension(:), intent(IN)    :: PCAN_HW_RATIO ! canyon    h/W
real, dimension(:), intent(IN)    :: PALB_ROOF     ! roof albedo
real, dimension(:), intent(IN)    :: PEMIS_ROOF    ! roof emissivity
real, dimension(:,:), intent(IN)  :: PHC_ROOF      ! heat capacity for roof layers
real, dimension(:,:), intent(IN)  :: PTC_ROOF      ! thermal conductivity for roof layers
real, dimension(:,:), intent(IN)  :: PD_ROOF       ! depth of roof layers
real, dimension(:), intent(IN)    :: PALB_ROAD     ! road albedo
real, dimension(:), intent(IN)    :: PEMIS_ROAD    ! road emissivity
real, dimension(:,:), intent(IN)  :: PHC_ROAD      ! heat capacity for road layers
real, dimension(:,:), intent(IN)  :: PTC_ROAD      ! thermal conductivity for road layers
real, dimension(:,:), intent(IN)  :: PD_ROAD       ! depth of road layers
real, dimension(:), intent(IN)    :: PSVF_ROAD     ! road sky view factor
real, dimension(:), intent(IN)    :: PALB_WALL     ! wall albedo
real, dimension(:), intent(IN)    :: PEMIS_WALL    ! wall emissivity
real, dimension(:,:), intent(IN)  :: PHC_WALL      ! heat capacity for wall layers
real, dimension(:,:), intent(IN)  :: PTC_WALL      ! thermal conductivity for wall layers
real, dimension(:,:), intent(IN)  :: PD_WALL       ! depth of wall layers
real, dimension(:), intent(IN)    :: PSVF_WALL     ! wall sky view factor
real, dimension(:), intent(IN)    :: PLAT          ! latitude

real, dimension(:), intent(OUT)   :: PRN_ROOF     ! net radiation over roof
real, dimension(:), intent(OUT)   :: PH_ROOF      ! sensible heat flux over roof
real, dimension(:), intent(OUT)   :: PLE_ROOF     ! latent heat flux over roof
real, dimension(:), intent(OUT)   :: PLES_ROOF    ! latent heat flux over roof (snow)
real, dimension(:), intent(OUT)   :: PGFLUX_ROOF  ! flux through the roof
real, dimension(:), intent(OUT)   :: PRUNOFF_ROOF ! runoff over the ground
real, dimension(:), intent(OUT)   :: PRN_ROAD     ! net radiation over road
real, dimension(:), intent(OUT)   :: PH_ROAD      ! sensible heat flux over road
real, dimension(:), intent(OUT)   :: PLE_ROAD     ! latent heat flux over road
real, dimension(:), intent(OUT)   :: PLES_ROAD    ! latent heat flux over road (snow)
real, dimension(:), intent(OUT)   :: PGFLUX_ROAD  ! flux through the road
real, dimension(:), intent(OUT)   :: PRUNOFF_ROAD ! runoff over the ground
real, dimension(:), intent(OUT)   :: PRN_WALL     ! net radiation over wall
real, dimension(:), intent(OUT)   :: PH_WALL      ! sensible heat flux over wall
real, dimension(:), intent(OUT)   :: PLE_WALL     ! latent heat flux over wall
real, dimension(:), intent(OUT)   :: PGFLUX_WALL  ! flux through the wall

real, dimension(:), intent(OUT)   :: PRNSNOW_ROOF ! net radiation over snow
real, dimension(:), intent(OUT)   :: PHSNOW_ROOF  ! sensible heat flux over snow
real, dimension(:), intent(OUT)   :: PLESNOW_ROOF ! latent heat flux over snow
real, dimension(:), intent(OUT)   :: PGSNOW_ROOF  ! flux under the snow
real, dimension(:), intent(OUT)   :: PMELT_ROOF   ! snow melt
real, dimension(:), intent(OUT)   :: PRNSNOW_ROAD ! net radiation over snow
real, dimension(:), intent(OUT)   :: PHSNOW_ROAD  ! sensible heat flux over snow
real, dimension(:), intent(OUT)   :: PLESNOW_ROAD ! latent heat flux over snow
real, dimension(:), intent(OUT)   :: PGSNOW_ROAD  ! flux under the snow
real, dimension(:), intent(OUT)   :: PMELT_ROAD   ! snow melt

real, dimension(:), intent(OUT)   :: PRN_TOWN     ! net radiation over town
real, dimension(:), intent(OUT)   :: PH_TOWN      ! sensible heat flux over town
real, dimension(:), intent(OUT)   :: PLE_TOWN     ! latent heat flux over town
real, dimension(:), intent(OUT)   :: PGFLUX_TOWN  ! flux through the ground
real, dimension(:), intent(OUT)   :: PEVAP_TOWN   ! evaporation flux (kg/m2/s)
real, dimension(:), intent(OUT)   :: PRUNOFF_TOWN ! runoff over the ground
real, dimension(:), intent(OUT)   :: PUSTAR_TOWN  ! friciton velocity over town
real, dimension(:), intent(OUT)   :: PCH_TOWN     ! town averaged heat transfer
!                                                 ! coefficient
real, dimension(:), intent(OUT)   :: PRI_TOWN      ! town averaged Richardson number
real, dimension(:), intent(OUT)   :: PTS_TOWN      ! town surface temperature
real, dimension(:), intent(OUT)   :: PEMIS_TOWN    ! town equivalent emissivity
real, dimension(:), intent(OUT)   :: PDIR_ALB_TOWN ! town equivalent direct albedo
real, dimension(:), intent(OUT)   :: PSCA_ALB_TOWN ! town equivalent diffuse albedo
real, dimension(:), intent(OUT)   :: PRESA_TOWN    ! town aerodynamical resistance

real, dimension(:), intent(OUT)   :: PTRAD_IN       ! body MRT inside building (K)
real, dimension(:), intent(OUT)   :: PTRAD_SUN     ! body MRT in the exposed street (K)
real, dimension(:), intent(OUT)   :: PTRAD_SHADE   ! body MRT in the shaded street (K)
real, dimension(:), intent(OUT)   :: PTRAD_RFSUN   ! body MRT on the exposed roof (K)
real, dimension(:), intent(OUT)   :: PTRAD_RFSHADE ! body MRT on the shaded roof (K)
real, dimension(:), intent(OUT)   :: PUTCI_IN       ! UTCI inside building
real, dimension(:), intent(OUT)   :: PUTCI_OUTSUN   ! UTCI in the exposed street
real, dimension(:), intent(OUT)   :: PUTCI_OUTSHADE ! UTCI in the shaded street
real, dimension(:), intent(OUT)   :: PUTCI_RFSUN    ! UTCI on the exposed roof
real, dimension(:), intent(OUT)   :: PUTCI_RFSHADE  ! UTCI on the shaded roof

real, dimension(:), intent(OUT)   :: PUTCIC_IN       ! UTCI inside building
real, dimension(:), intent(OUT)   :: PUTCIC_OUTSUN   ! UTCI in the exposed street
real, dimension(:), intent(OUT)   :: PUTCIC_OUTSHADE ! UTCI in the shaded street
real, dimension(:), intent(OUT)   :: PUTCIC_RFSUN   ! UTCI in the exposed street
real, dimension(:), intent(OUT)   :: PUTCIC_RFSHADE ! UTCI in the shaded street

real, dimension(:), intent(OUT)   :: WBGT_SUN         ! WBGT  wet bulb globe temperature in the street (C)
real, dimension(:), intent(OUT)   :: WBGT_SHADE       ! WBGT  wet bulb globe temperature in the street (C)
real, dimension(:), intent(OUT)   :: WBGT_RFSUN       ! WBGT  wet bulb globe temperature on the roof (C)
real, dimension(:), intent(OUT)   :: WBGT_RFSHADE     ! WBGT  wet bulb globe temperature on the roof (C)

 real, dimension(:), intent(OUT)   :: PTGLOBE_SUN      ! Globe Temperature in the exposed street (K)
 real, dimension(:), intent(OUT)   :: PTGLOBE_SHADE    ! Globe Temperature in the shaded street (K)
 real, dimension(:), intent(OUT)   :: PTGLOBE_RFSUN    ! Globe Temperature  on the exposed roof (K)
 real, dimension(:), intent(OUT)   :: PTGLOBE_RFSHADE  ! Globe Temperature  on the shaded roof  (K)
 real, dimension(:), intent(OUT)   :: PTWETB           ! wet-bulb Temperature  on the ground (K)
 real, dimension(:), intent(OUT)   :: PTWETB_ROOF      ! wet-bulb temperature on the  roof  (K)

 real, dimension(:), intent(OUT)   :: PTRFZT           ! T at z=zt above the RooF
 real, dimension(:), intent(OUT)   :: PTRDZT           ! T at z=zt above the RoaD
 real, dimension(:), intent(OUT)   :: PURDZU           ! U at z=zu above the RoaD

 real, dimension(:), intent(OUT)   :: PQ1       ! energy fluxes received by a body (here for a human)
 real, dimension(:), intent(OUT)   :: PQ2
 real, dimension(:), intent(OUT)   :: PQ3
 real, dimension(:), intent(OUT)   :: PQ4
 real, dimension(:), intent(OUT)   :: PQ5
 real, dimension(:), intent(OUT)   :: PQ6
 real, dimension(:), intent(OUT)   :: PQ7
 real, dimension(:), intent(OUT)   :: PQ8
 real, dimension(:), intent(OUT)   :: PQ9
 real, dimension(:), intent(OUT)   :: PQ10
 real, dimension(:), intent(OUT)   :: PQ11
 real, dimension(:), intent(OUT)   :: PQ12
 real, dimension(:), intent(OUT)   :: PQ13


!*      0.2    declarations of local variables

logical                    :: GTI_EVOL       ! true --> internal temperature
!                                            !          of buildings evolves
!                                            ! false--> it is fixed

real, dimension(size(PTA)) :: PTRAD_GSUN     ! globe MRT in the exposed street (K)
real, dimension(size(PTA)) :: PTRAD_GSHADE   ! globe MRT in the shaded street (K)
real, dimension(size(PTA)) :: PTRAD_GRFSUN   ! globe MRT on the exposed roof (K)
real, dimension(size(PTA)) :: PTRAD_GRFSHADE ! globe MRT on the shaded roof (K)
real, dimension(size(PTA)) :: ZVMOD          ! wind
real, dimension(size(PTA)) :: ZWS_ROOF_MAX   ! maximum deepness of roof
real, dimension(size(PTA)) :: ZWS_ROAD_MAX   ! and road water reservoirs


real, dimension(size(PTA)) :: ZDELT_ROOF     ! water fraction
real, dimension(size(PTA)) :: ZDELT_ROAD     ! on the surface
real, dimension(size(PTA)) :: ZDN_ROOF       ! snow fraction
real, dimension(size(PTA)) :: ZDN_ROAD       ! on the surface
real, dimension(size(PTA)) :: ZDF_ROOF       ! free-snow fraction
real, dimension(size(PTA)) :: ZDF_ROAD       ! on the surface

real, dimension(size(PTA)) :: ZQSAT_ROAD     ! q_sat(Ts)
real, dimension(size(PTA)) :: ZQSAT_ROOF     ! q_sat(Ts)
real, dimension(size(PTA)) :: ZAC_ROOF       ! surface conductance
!                                            ! for heat transfers
!                                            ! above roofs
real, dimension(size(PTA)) :: ZAC_ROOF_WAT   ! surface conductance
!                                            ! for heat transfers
!                                            ! above roofs (for water)
real, dimension(size(PTA)) :: ZAC_WALL       ! surface conductance
!                                            ! for heat transfers
!                                            ! between wall and canyon air
real, dimension(size(PTA)) :: ZAC_ROAD       ! surface conductance
!                                            ! for heat transfers
!                                            ! between road and canyon air
real, dimension(size(PTA)) :: ZAC_ROAD_WAT   ! surface conductance
!                                            ! for heat transfers
!                                            ! inside canyon (for water)
real, dimension(size(PTA)) :: ZAC_TOP        ! surface conductance
!                                            ! for heat transfers
!                                            ! canyon top and atm.
real, dimension(size(PTA)) :: ZAC_BLD        ! surface conductance inside
                                             ! the building itself
real, dimension(size(PTA)) :: ZCD            ! drag coefficient
real, dimension(size(PTA)) :: ZTA            ! air temperature extrapolated at roof level
real, dimension(size(PTA)) :: ZQA            ! air humidity extrapolated at roof level

real, dimension(size(PTA)) :: ZROOF          ! roof, wall and
real, dimension(size(PTA)) :: ZWALL          ! road fractions
real, dimension(size(PTA)) :: ZROAD          ! of exchange surf.
real, dimension(size(PTA)) :: ZTOTS_O_HORS   ! total canyon+roof surface
!                                            ! over horizontal surface
real, dimension(size(PTA)) :: ZWALL_O_ROAD   ! wall surface over road surface


! absorbed solar and infra-red radiation by road, wall and roof

real, dimension(size(PTA)) :: ZABS_SW_ROAD
real, dimension(size(PTA)) :: ZABS_SW_WALL
real, dimension(size(PTA)) :: ZABS_SW_ROOF
real, dimension(size(PTA)) :: ZABS_SW_SNOW_ROAD
real, dimension(size(PTA)) :: ZABS_SW_SNOW_ROOF
real, dimension(size(PTA)) :: ZABS_LW_SNOW_ROAD
real, dimension(size(PTA)) :: ZABS_LW_SNOW_ROOF
real, dimension(size(PTA)) :: ZABS_LW_ROAD
real, dimension(size(PTA)) :: ZABS_LW_WALL
real, dimension(size(PTA)) :: ZABS_LW_ROOF

! reflected solar and infra-red radiation by road, wall and roof

real, dimension(size(PTA)) :: ZREF_SW_ROAD
real, dimension(size(PTA)) :: ZREF_SW_WALL
real, dimension(size(PTA)) :: ZREF_SW_ROOF
real, dimension(size(PTA)) :: ZREF_SW_SNOW_ROAD
real, dimension(size(PTA)) :: ZREF_SW_SNOW_ROOF
real, dimension(size(PTA)) :: ZREF_SW_TOT_ROAD
real, dimension(size(PTA)) :: ZREF_SW_TOT_ROOF
real, dimension(size(PTA)) :: ZEMI_LW_SNOW_ROAD
real, dimension(size(PTA)) :: ZEMI_LW_SNOW_ROOF
real, dimension(size(PTA)) :: ZEMI_LW_TOT_ROAD
real, dimension(size(PTA)) :: ZEMI_LW_TOT_ROOF
real, dimension(size(PTA)) :: ZEMI_LW_ROAD
real, dimension(size(PTA)) :: ZEMI_LW_WALL
real, dimension(size(PTA)) :: ZEMI_LW_ROOF

real, dimension(size(PTA)) :: ZU10           ! wind speed at 10 m AGL
real, dimension(size(PTA)) :: ZUSR           ! wind speed at sensor level at ~2.5 m ARL
real, dimension(size(PTA)) :: PURFZU         ! U at z=zu above the RooF (zu=zt in urban_drag)
real, dimension(size(PTA)) :: PVRFZU         ! V at z=zu above the RooF (zu=zt in urban_drag)

! coefficients for LW computations over snow (from previous time-step)

real, dimension(size(PTA)) :: ZTSSNOW_ROAD ! road snow temperature
!                                          ! at previous time-step
real, dimension(size(PTA)) :: ZESNOW_ROOF  ! snow emissivity
!                                          ! at previous time-step
real, dimension(size(PTA)) :: ZESNOW_ROAD  ! snow emissivity
!                                          ! at previous time-step
real, dimension(size(PTA)) :: ZTSSNOW_ROOF ! road snow temperature
!                                          ! at previous time-step

! coefficients for contributions of sky long wave to wall, road, snow on road

real, dimension(size(PTA)) :: ZLW_S_TO_W
real, dimension(size(PTA)) :: ZLW_S_TO_R
real, dimension(size(PTA)) :: ZLW_S_TO_N


integer :: IWALL           ! number of wall layers
integer :: IROOF           ! number of roof layers
!-------------------------------------------------------------------------------
integer :: i
      real UNDEF
      parameter (UNDEF = -999)
!-------------------------------------------------------------------------------

!*      1.     initializations
!              ---------------

!*      1.1    water reservoirs
!              ----------------

ZWS_ROOF_MAX =  1. ! (1mm) maximum deepness of roof water reservoir
ZWS_ROAD_MAX =  1. ! (1mm) maximum deepness of road water reservoir

!*      1.2    surfaces relative fractions
!              ---------------------------

ZTOTS_O_HORS= 1. + PWALL_O_HOR(:)

ZROOF=PBLD        /ZTOTS_O_HORS
ZWALL=PWALL_O_HOR /ZTOTS_O_HORS
ZROAD=(1.-PBLD)   /ZTOTS_O_HORS

ZWALL_O_ROAD = ZWALL / ZROAD

!*       1.3    option for internal temperature of buildings
!               --------------------------------------------

GTI_EVOL = .true.


!*       1.4    set wind to a minimum value
!               ---------------------------

ZVMOD(:) = max(1.,PVMOD(:))

!-------------------------------------------------------------------------------

!*      2.     snow-covered surfaces relative effects
!              --------------------------------------
!       Impose minimum snow value, so no calculation are done on tiny reservoir

where (PWSNOW_ROAD(:)<SWE_CRIT) PWSNOW_ROAD(:)=0.0
where (PWSNOW_ROOF(:)<SWE_CRIT) PWSNOW_ROOF(:)=0.0



!*      2.1    snow-covered surfaces relative fractions (at previous time-step)
!              ----------------------------------------


call SNOW_FRAC_ROAD(PWSNOW_ROAD(:),PSR(:)>0.,ZDN_ROAD,ZDF_ROAD)
call SNOW_FRAC_ROOF(PWSNOW_ROOF(:),PSR(:)>0.,ZDN_ROOF,ZDF_ROOF)

!* new snow albedo

where (PWSNOW_ROAD(:)==0. .and. PSR(:)>0.) PASNOW_ROAD(:) = XANSMAX
where (PWSNOW_ROOF(:)==0. .and. PSR(:)>0.) PASNOW_ROOF(:) = XANSMAX
        

!*      2.2    effects on water reservoirs
!              ---------------------------

ZWS_ROOF_MAX(:) = ZWS_ROOF_MAX(:) * ZDF_ROOF(:)
ZWS_ROAD_MAX(:) = ZWS_ROAD_MAX(:) * ZDF_ROAD(:)

!*      2.3    if snow was not present at previous time-step but is falling
!              ------------------------------------------------------------

where (PWSNOW_ROAD(:)==0. .and. PSR(:)>0.)
  PASNOW_ROAD(:) = XANSMAX
  PESNOW_ROAD(:) = XEMISSN
  PTSSNOW_ROAD(:)= min(PT_ROAD(:,1), XTT)
end where
where (PWSNOW_ROOF(:)==0. .and. PSR(:)>0.)
  PASNOW_ROOF(:) = XANSMAX
  PESNOW_ROOF(:) = XEMISSN
  PTSSNOW_ROOF(:)= min(PT_ROOF(:,1), XTT)
end where

!*      2.4    radiative snow variables at previous time-step
!              ----------------------------------------------

ZESNOW_ROAD (:)=PESNOW_ROAD (:)
ZTSSNOW_ROAD(:)=PTSSNOW_ROAD(:)

ZESNOW_ROOF (:)=PESNOW_ROOF (:)
ZTSSNOW_ROOF(:)=PTSSNOW_ROOF(:)

!-------------------------------------------------------------------------------

!*      3.     SOLAR RADIATIONS
!              ----------------

!* for sake of simplicity, the radiative properties of the vegetation (if any)
!  situated in the canyon (and then influencing h/w) is replaced by the road
!  properties. This influences only the radiative part due to reflections,
!  which are less important when h/w is low (i.e. large streets with garden),
!  than when h/w is large (from 0.5 and above approximately). This case
!  (h/w large) occurs in city downtowns, and then does not occur in presence
!  of significant vegetation area.

call URBAN_SOLAR_ABS(PDIR_SW_RAD, PSCA_SW_RAD, PZENITH,            &
                     PBLD, PWALL_O_HOR, PCAN_HW_RATIO,             &
                     PALB_ROOF,                                    &
                     PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
                     PASNOW_ROOF, PASNOW_ROAD,                     &
                     ZDN_ROOF, ZDF_ROOF, ZDN_ROAD, ZDF_ROAD,       &
                     ZABS_SW_ROOF, ZABS_SW_ROAD, ZABS_SW_WALL,     &
                     ZABS_SW_SNOW_ROOF, ZABS_SW_SNOW_ROAD,         &
                     PDIR_ALB_TOWN,PSCA_ALB_TOWN,                  &
                     ZREF_SW_ROOF, ZREF_SW_ROAD, ZREF_SW_WALL,     &
                     ZREF_SW_SNOW_ROOF, ZREF_SW_SNOW_ROAD          )

! added
where (PWSNOW_ROAD(:)<SWE_CRIT) ZREF_SW_SNOW_ROAD(:)=0.0
where (PWSNOW_ROOF(:)<SWE_CRIT) ZREF_SW_SNOW_ROOF(:)=0.0

ZREF_SW_TOT_ROAD(:) = ZDF_ROAD(:) * ZREF_SW_ROAD(:) + ZDN_ROAD(:) * ZREF_SW_SNOW_ROAD(:)
ZREF_SW_TOT_ROOF(:) = ZDF_ROOF(:) * ZREF_SW_ROOF(:) + ZDN_ROOF(:) * ZREF_SW_SNOW_ROOF(:)

!-------------------------------------------------------------------------------

!*      4.     drag accounting for the effets of town temperature and humidity
!              ---------------------------------------------------------------

call URBAN_DRAG2(PTSTEP, PT_CANYON, PQ_CANYON,                       &
                PT_ROOF(:,1), PT_ROAD(:,1), PT_WALL(:,1), ZDN_ROOF,  &
                PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,                  &
                PZREF, PUREF, ZVMOD, PVDIR, PLAT,                    &
                PZ0_TOWN,PZ0_ROOF,PZ0_ROAD,                          &
                PBLD, PBLD_HEIGHT, PCAN_HW_RATIO,                    &
                ZWALL_O_ROAD,                                        &
                PWS_ROOF, PWS_ROAD,                                  &
                ZWS_ROOF_MAX, ZWS_ROAD_MAX,                          &
                ZQSAT_ROOF, ZQSAT_ROAD, ZDELT_ROOF, ZDELT_ROAD,      &
                ZCD, ZAC_ROOF, ZAC_ROOF_WAT,                         &
                ZAC_WALL, ZAC_ROAD, ZAC_ROAD_WAT, ZAC_TOP,           &
                PU_CANYON, PRI_TOWN,                                 &
                PTRDZT,PTRFZT,PURFZU,PVRFZU                        )

!* area-averaged heat transfer coefficient

PCH_TOWN(:) =  (PBLD(:) * ZAC_ROOF(:) + (1.-PBLD(:)) * ZAC_TOP (:)) / ZVMOD(:)

!-------------------------------------------------------------------------------

!*      5.     Extrapolation of atmospheric T and q at roof level (for fluxes computation)
!              --------------------------------------------------

ZTA(:) = PTA(:) * PEXNS(:) / PEXNA(:)
ZQA(:) = PQA(:) * QSAT(PTA(:),PPS(:)) / QSAT(ZTA(:),PPA(:))

!-------------------------------------------------------------------------------

!*      6.     Evolution of interior building air temperature
!              ----------------------------------------------

call BLD_E_BUDGET(GTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,             &
                  PRHOA, PT_ROOF, PT_WALL, PTI_BLD, ZAC_BLD        )

!-------------------------------------------------------------------------------

!*      7.     Snow mantel model
!              -----------------



call URBAN_SNOW_EVOL(PT_CANYON, PQ_CANYON, PU_CANYON,                         &
                     PT_ROOF(:,1),PT_ROAD(:,1),PT_WALL(:,1),                  &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPS, ZTA, ZQA, PRHOA,                                    &
                     PLW_RAD,                                                 &
                     PSR, PZREF, PUREF, ZVMOD,                                &
                     PTSTEP,                                                  &
                     PBLD_HEIGHT,                                             &
                     PEMIS_ROAD, PSVF_ROAD,                                   &
                     PEMIS_WALL, PSVF_WALL,                                   &
                     ZDN_ROOF, ZABS_SW_SNOW_ROOF, ZABS_LW_SNOW_ROOF,          &
                     ZDF_ROAD, ZDN_ROAD, ZABS_SW_SNOW_ROAD, ZABS_LW_SNOW_ROAD,&
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     ZLW_S_TO_N                                               )


!-------------------------------------------------------------------------------

!*      8.     Roof Ts computation
!              -------------------

!* ts_roof and qsat_roof are updated

call ROOF_LAYER_E_BUDGET(PT_ROOF, ZQSAT_ROOF,                              &
                         ZTA, ZQA, PPS,                                    &
                         PLW_RAD, PTSTEP,                                  &
                         PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,          &
                         PTI_BLD, ZAC_BLD, ZDELT_ROOF,                     &
                         ZDN_ROOF, PGSNOW_ROOF,                            &
                         PRHOA, ZAC_ROOF, ZAC_ROOF_WAT,                    &
                         ZABS_SW_ROOF, ZABS_LW_ROOF                        )


!-------------------------------------------------------------------------------


!*      9.    Road and wall Ts computations
!             -----------------------------

!* ts_road, ts_wall, qsat_road, t_canyon and q_canyon are updated

call ROAD_WALL_LAYER_E_BUDGET(PT_ROAD,                                 &
                              PT_WALL, ZQSAT_ROAD,                     &
                              PT_CANYON, PQ_CANYON,                    &
                              ZTA, ZQA, PPS,                           &
                              PLW_RAD,  PTSTEP,                        &
                              PH_TRAFFIC, PLE_TRAFFIC,                 &
                              PBLD, ZWALL_O_ROAD,                      &
                              PEMIS_ROAD, PSVF_ROAD,                   &
                              PHC_ROAD, PTC_ROAD, PD_ROAD,             &
                              PEMIS_WALL, PSVF_WALL,                   &
                              PHC_WALL, PTC_WALL, PD_WALL,             &
                              PTI_BLD, ZAC_BLD,                        &
                              ZDELT_ROAD, ZDN_ROAD,                    &
                              ZTSSNOW_ROAD, ZESNOW_ROAD,               &
                              PGSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD,  &
                              PRHOA, ZAC_WALL,                         &
                              ZAC_ROAD, ZAC_ROAD_WAT, ZAC_TOP,         &
                              ZABS_SW_ROAD, ZABS_SW_WALL,              &
                              ZABS_LW_ROAD, ZABS_LW_WALL,              &
                              ZLW_S_TO_R, ZLW_S_TO_W                   )


!-------------------------------------------------------------------------------

!*     10.     Fluxes
!              ------

call URBAN_FLUXES   (PTS_TOWN, PEMIS_TOWN,                                    &
                     PT_CANYON, PQ_CANYON,                                    &
                     PT_ROOF(:,1),PT_ROAD(:,1),PT_WALL(:,1),                  &
                     ZTA, ZQA, PRHOA,                                         &
                     ZVMOD,                                                   &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PBLD,                                                    &
                     PEMIS_ROOF,                                              &
                     ZESNOW_ROOF,                                             &
                     PLW_RAD, ZLW_S_TO_R, ZLW_S_TO_W, ZLW_S_TO_N,             &
                     ZABS_SW_ROOF, ZABS_LW_ROOF, ZABS_SW_WALL, ZABS_LW_WALL,  &
                     ZABS_SW_ROAD, ZABS_LW_ROAD,                              &
                     ZABS_LW_SNOW_ROOF,                                       &
                     ZABS_LW_SNOW_ROAD,                                       &
                     ZAC_ROOF, ZAC_ROOF_WAT,                                  &
                     ZAC_WALL, ZAC_ROAD, ZAC_ROAD_WAT, ZCD,                   &
                     ZQSAT_ROOF, ZQSAT_ROAD,                                  &
                     ZDELT_ROOF, ZDELT_ROAD,                                  &
                     ZROOF, ZWALL, ZROAD, ZTOTS_O_HORS,                       &
                     ZDF_ROOF, ZDN_ROOF, ZDF_ROAD, ZDN_ROAD,                  &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLES_ROOF, PGFLUX_ROOF,     &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLES_ROAD, PGFLUX_ROAD,     &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,    &
                     PUSTAR_TOWN                                              )

!-------------------------------------------------------------------------------

!*     11.     Roof and road reservoirs evolution
!              ----------------------------------

call URBAN_HYDRO(ZWS_ROOF_MAX,ZWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,        &
                 PRR, PTSTEP, PBLD, PLE_ROOF, PLE_ROAD,                &
                 PRUNOFF_ROOF,                                         &
                 PRUNOFF_ROAD,                                         &
                 PRUNOFF_TOWN                                          )

!-------------------------------------------------------------------------------

!*     12.     Compute aerodynamical resistance
!              --------------------------------

PRESA_TOWN(:) = 1. / ( PBLD(:) * ZAC_ROOF(:)  + ( 1. - PBLD(:)) * ZAC_TOP (:))


!*     13.     Compute thermal indices
!              --------------------------------
   IF_THERMAL_STRESS: if ( thermal_stress ) then
! input: Wind module effectively at 10 m AGL 
! note ; we may take zu instead of 10 ? but indices built with 10-m wind speed...
! Compute wind at z=zt
! check wind at z=zt 
!  i = sl_sfclayer(th,hu,vmod,vdir,zzusl,zztsl,sst,qs,z0m,z0h,zdlat,zfcor, &
!        hghtm_diag=zt,hghtt_diag=zt,u_diag=zusurfzt,v_diag=zvsurfzt, &
!        tdiaglim=WATER_TDIAGLIM)

!   if (i /= SL_OK) then
!      print*, 'Aborting in wate teb() because of error returned by sl_sfclayer()'
!      stop
!   endif

! => wind, udiag either from sl_sfclayer (above canopy) or adapted for street morphology (urb_diagwind)
! => wind
! if bldh>15m => u_canyon 
 do i=1,size(pta)
 if (pbld_height(i) .ge. 15.) then 
        zu10(i) = pu_canyon(i) 
 endif
!      if(urb_diagwind) then
!  linear interpolation between two cases du/dz=Utop-Ucan / H-2H/3
!                           => U(z=10)= U10 = (30/H-2) Utop + (-30/H+3) Ucan
      if( pbld_height(i) .ge. 10. .and. pbld_height(i) .lt. 15. )  then 
        zu10(i) =  (30.0 /pbld_height(i) -2.) * pvmod(i)              &
           * log( (     pbld_height(i)/3.) / pz0_town (i))            &
           / log( (puref(i) + pbld_height(i)/3.) / pz0_town   (i))    &
           + (- 30. /pbld_height(i) +3.) * pu_canyon(i)
      elseif (pbld_height(i) .lt. 10.0 ) then
!                            => log law. above roof level 
  zu10(i) =  pvmod(i)                                                   &
           * log( ( 10. -2. * pbld_height(i)/3.) / pz0_town(i))   &
           / log( ( puref(i)  + pbld_height(i)/3.) / pz0_town(i))
       endif
 !endif
 enddo

!* IMPOSE MINIMUM OF 0.01 for PU10 
 PURDZU(:) = MAX(ZU10(:) ,0.01)

! wind at the sensor level above the roof level (2.5m)
zusr(:) = max(   pvmod(:)  &
            * log(      (2.5  + 1.* pbld_height(:)/3.) / pz0_town(:) )    &
            / log( (puref(:) + 1.* pbld_height(:)/3.) / pz0_town(:) ) ,0.01) 

! input: longwave infrared rad
 ZEMI_LW_WALL(:) = XSTEFAN * PEMIS_WALL(:) * PT_WALL(:,1)**4
 ZEMI_LW_ROAD(:) = XSTEFAN * PEMIS_ROAD(:) * PT_ROAD(:,1)**4
 ZEMI_LW_ROOF(:) = XSTEFAN * PEMIS_ROOF(:) * PT_ROOF(:,1)**4
! add snow contribution
 ZEMI_LW_SNOW_ROAD(:) = XSTEFAN * ZESNOW_ROAD(:) * ZTSSNOW_ROAD(:)**4
 ZEMI_LW_SNOW_ROOF(:) = XSTEFAN * ZESNOW_ROOF(:) * ZTSSNOW_ROOF(:)**4
! total contrib ?
  ZEMI_LW_TOT_ROAD(:) = ZDF_ROAD(:) * ZEMI_LW_ROAD(:) + ZDN_ROAD(:) * ZEMI_LW_SNOW_ROAD(:)
  ZEMI_LW_TOT_ROOF(:) = ZDF_ROOF(:) * ZEMI_LW_ROOF(:) + ZDN_ROOF(:) * ZEMI_LW_SNOW_ROOF(:)

! input:  radiative temperature inside buildings (no BEM version) check 19
!       ZTS_FLOOR(:) = 19. + XTT
 IWALL = size(PT_WALL,2)
 IROOF = size(PT_ROOF,2)

      PTRAD_IN(:) = (PWALL_O_HOR(:) / PBLD(:) * PT_WALL(:,IWALL)+ &
                     PT_ROOF(:,IROOF) + (23. + XTT)  )             &
                    / (PWALL_O_HOR(:) / PBLD(:) + 1. + 1.)

 call URBAN_THERMAL_STRESS(PT_CANYON, PQ_CANYON, PTI_BLD, PQ_CANYON,    &
             PU_CANYON,PURDZU, ZUSR, ZVMOD, PPS, PPA,                   &
             ZREF_SW_TOT_ROAD, ZREF_SW_WALL, PSCA_SW_RAD, PDIR_SW_RAD,  &
             PZENITH, ZEMI_LW_WALL, ZEMI_LW_TOT_ROAD,PLW_RAD,PTRAD_IN,  &
              ZTA, ZQA, ZREF_SW_TOT_ROOF, ZEMI_LW_TOT_ROOF,             &
             PBLD, PBLD_HEIGHT, PWALL_O_HOR,                            &
             PUTCI_IN, PUTCI_OUTSUN, PUTCI_OUTSHADE,                    &
             PUTCI_RFSUN, PUTCI_RFSHADE,                                &
             WBGT_SUN,WBGT_SHADE,WBGT_RFSUN,WBGT_RFSHADE,               &
             PTRAD_SUN, PTRAD_SHADE, PTRAD_RFSUN, PTRAD_RFSHADE,        &
             PTRAD_GSUN, PTRAD_GSHADE, PTRAD_GRFSUN, PTRAD_GRFSHADE,    &
             PTGLOBE_SUN, PTGLOBE_SHADE,                                &
             PTGLOBE_RFSUN,PTGLOBE_RFSHADE,PTWETB, PTWETB_ROOF,         &
             PQ1,PQ2,PQ3,PQ4,PQ5,PQ6,PQ7,PQ8,PQ9,PQ10,PQ11,PQ12,PQ13 )

!not used for now
!CALL UTCIC_STRESS(PTSTEP, PUTCI, PUTCIC )
!CALL UTCIC_STRESS(PTSTEP, PUTCI_IN, PUTCIC_IN )
!CALL UTCIC_STRESS(PTSTEP, PUTCI_OUTSUN, PUTCIC_OUTSUN )
!CALL UTCIC_STRESS(PTSTEP, PUTCI_OUTSHADE, PUTCIC_OUTSHADE )
!CALL UTCIC_STRESS(PTSTEP, PUTCI_RFSUN, PUTCIC_RFSUN )
!CALL UTCIC_STRESS(PTSTEP, PUTCI_RFSHADE, PUTCIC_RFSHADE)

PUTCIC_IN(:)       = PUTCI_IN(:)
PUTCIC_OUTSUN(:)   = PUTCI_OUTSUN(:)
PUTCIC_OUTSHADE(:) = PUTCI_OUTSHADE(:)
PUTCIC_RFSUN(:)    = PUTCI_RFSUN(:)
PUTCIC_RFSHADE(:)  = PUTCI_RFSHADE(:)

 endif IF_THERMAL_STRESS
!-------------------------------------------------------------------------------

end subroutine TEB2
