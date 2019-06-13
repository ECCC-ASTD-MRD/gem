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

subroutine SOIL_FLUXES(DT, &
         WSATC, KSATC, PSISAT, BCOEF, ETR_GRID, WD, &
         F, WDT, DWD, OVRSHT, KHC, PSI, N)

        use sfc_options
        use svs_configs
      implicit none
!!!#include <arch_specific.hf>

      integer N,K
      ! **NUMERICAL** SATURATION MINIMUM of 1% TO AVOID UNDERFLOW/OVERFLOW
      ! PROBLEMS in SOIL PARAMETER PARAMETRIZATIONS
      real, parameter :: NUM_MIN=0.01

      ! input
      real :: dt
      real, dimension(n,nl_svs) :: wsatc, ksatc, psisat, bcoef, etr_grid, wd
      ! input/output
      real, dimension(n,nl_svs+1) :: f
      ! output
      real, dimension(n,nl_svs) :: wdt, dwd, ovrsht
      real, dimension(n,nl_svs-1):: khc, psi


!Author
!           V. Fortin (October 2016)

!Object
!     Calculate the flux of water between soil layers
!     based on soil characteristics and liquid water content
!     using the Brooks and Corey model and updates the water content
!     this function uses a Euler forward 1st order method to solve
!     for water content at the end of the time step but can be embedded
!     in a Runge-Kutta algorithm to increase the order of the scheme

!Arguments

!          - INPUT -

! DT                model time step used in the computation of F [s]

!          --- Soil characteristics (near saturation) for each layer ---

! WSATC (NL_SVS)    volumetric water content at soil saturation per layer [m3/m3]
! KSATC (NL_SVS)    vertical hydraulic conductivity at saturation per layer [m/s]
! PSISAT(NL_SVS)    value of soil water suction at air-entry (near saturation) per layer [m]
! BCOEF (NL_SVS)    slope of the retention curve per layer

!          --- Source/sink terms ---

! ETR_GRID(NL_SVS)  evapotranspiration rate for each layer, grid box average [m/s]

!          --- Current water content for each layer ---

! WD (NL_SVS)       soil volumetric water content (per layer) [m3/m3]

!          -  INPUT/OUTPUT  -

! F (NL_SVS+1)      total flux of water between each layer and at top and bottom boundaries during time step [m]
!                   initially, contains three values:
!                   F(1):       infiltration coming in at the top of the soil column
!                   F(KDP+1):   baseflow along impervious layer (if KDP<NL_SVS)
!                   F(NL_SVS+1):baseflow at the bottom of the soil column
!                   these three values are unchanged by this subroutine, which fills in the other values
!                   * THESE THREE VALUES OF F NEED TO BE PROPERLY INITIALIZED BY THE CALLING SUBROUTINE *

!          -  OUTPUT  -

!          --- Water content and flux of water between layers (volumetric) ---

! WDT (NL_SVS)      water content in each layer at the end of the time step [m3/m3]
! DWD (NL_SVS)      change in storage in each layer based on flux F for a time step DT [m3/m3]
! OVRSHT (NL_SVS)   overshoot correction required to keep WDT within physical bounds [m3/m3]
!                   this amount of water needs to be reallocated to maintain the water balance

!          ---  Diagnostic Soil Parameters (used only for testing SVS) ---

! KHC (NL_SVS-1)    soil hydraulic conductivity at soil boundaries (of layer) [m/s]
! PSI (NL_SVS-1)    soil water potential at the soil boundaries (of layer) [m]

!          -  DIMENSIONS  -

! N               number of grid cells


      integer I

      ! local arrays

      real, dimension(n,nl_svs-1) :: wsatbnd, dwddz_wdbnd, wdbnd, bbnd, ksatbnd, psisatbnd, satbnd

!***********************************************************************

!Compute water fluxes QZ between soil layers, find KSAT, PSISAT, K AND PSI at the boundaries
! do it for all layers except KDP and NL_SVS (water flux is computed above from watdrain).

      do K=1,NL_SVS-1
         do I=1,N
           !wsat at soil boundaries
            WSATBND(I,K)=(WSATC(I,K)+WSATC(I,K+1))/2.
           !WD at soil boundaries
            WDBND(I,K)= max((WD(I,K)+WD(I,K+1))/2.0,CRITWATER)
           !gradient of soil water content divided by water content
            DWDDZ_WDBND(I,K)=(WD(I,K+1)/WDBND(I,K)-1.)/DELZ(K+1)+ &
                 (1.-WD(I,K)/WDBND(I,K))/DELZ(K)
           !b-coefficient at the boundaries
            BBND(I,K)=(BCOEF(I,K)+BCOEF(I,K+1))/2.
           !ksat at soil boundaries
            KSATBND(I,K)=KSATC(I,K)*KSATC(I,K+1)*(DELZ(K)+&
                 DELZ(K+1))/(KSATC(I,K)*DELZ(K+1)+KSATC(I,K+1)*DELZ(K))
           !psisat at soil boundaries
            PSISATBND(I,K)=PSISAT(I,K)**(DELZ(K)/(DELZ(K)+&
                 DELZ(K+1)))*PSISAT(I,K+1)**(DELZ(K+1)/(DELZ(K)+DELZ(K+1)))
           !saturation at boundaries, impose numerical minimum to avoid underflow
            SATBND(I,K) =  max( WDBND(I,K)/WSATBND(I,K) , NUM_MIN )
           !soil hydraulic conductivity at soil boundaries
            KHC(I,K)= min(KSATBND(I,K)*SATBND(I,K)**(2.*BBND(I,K)+3.0),KSATBND(I,K))
           !soil water potential at the soil boundries
            PSI(I,K)= max(PSISATBND(I,K)*SATBND(I,K)**(-BBND(I,K)),PSISATBND(I,K))
           !vertical flux of water between soil layers
            if(K.ne.KDP)then
               F(I,K+1)=DT*KHC(I,K)*(-BBND(I,K)*PSI(I,K)*DWDDZ_WDBND(I,K)+1.0)
            end if
         end do
      end do

! Compute change in storage for each layer layer for a time step of DT seconds
! and update the water content

      do K=1,NL_SVS
         do I=1,N
            DWD(I,K)=(F(I,K)-F(I,K+1))/DELZ(K)-DT*ETR_GRID(I,K)
            ! Make sure DWD doesn't lead to WDT exceeding [CRITWATER,WSATC] range
            ! but keep track of the overshoot in order to be able to close the water balance
            OVRSHT(I,K)= min(max(DWD(I,K),CRITWATER-WD(I,K)),WSATC(I,K)-WD(I,K)) - DWD(I,K)
            DWD(I,K) = OVRSHT(I,K) + DWD(I,K)
            WDT(I,K)=WD(I,K)+DWD(I,K)
         end do
      end do

return
end subroutine SOIL_FLUXES
