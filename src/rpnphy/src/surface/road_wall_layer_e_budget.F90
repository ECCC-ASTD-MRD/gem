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
!   ##########################################################################
!
!!****  *ROAD_WALL_LAYER_E_BUDGET*
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of roads and walls surface temperatures
!
!
!!**  METHOD
!     ------
!
!    6 : equations for evolution of Ts_road and Ts_wall simultaneously
!        *************************************************************
!
!     dTw_k(t) / dt = 1/(dw_k*Cw_k) * (- 2*Kw_k-1*(Tw_k-Tw_k-1)/(dw_k-1 +dw_k)
!                                      - 2*Kw_k  *(Tw_k-Tw_k+1)/(dw_k+1 +dw_k) )
!
!     dTw_1(t) / dt = 1/(dw_1*Cw_1) * (  Rn_w - H_w - LE_w
!                                      - 2*Kw_1*(Tw_1-Tw_2)/(dw_1 +dw_2)       )
!
!     dTr_1(t) / dt = 1/(dr_1*Cr_1) * (  Rn_r - H_r - LE_r
!                                      - 2*Kr_1*(Tr_1-Tr_2)/(dr_1 +dr_2)       )
!
!     dTr_k(t) / dt = 1/(dr_k*Cr_k) * (- 2*Kr_k-1*(Tr_k-Tr_k-1)/(dr_k-1 +dr_k)
!                                      - 2*Kr_k  *(Tr_k-Tr_k+1)/(dr_k+1 +dr_k) )
!
!       with
!
!   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!   Rn_w = abs_Rg_w
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * LWR
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
!  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * LWR
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
!  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
!
!  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
!
!  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
!
!  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
!
! with again
!                AC_can * Swall/Sroad * Twall + AC_can * Troad + AC_top * Ta + H_traffic/rho/Cp/Sroad
!   Ta_canyon = -------------------------------------------------------------------------------------
!                AC_can * Swall/Sroad         + AC_can         + AC_top
!
!
!                 AC_can * delt_road * qsat(Troad) + AC_top * qa + LE_traffic/rho/Lv/Sroad
!   qa_canyon = --------------------------------------------------------------------------
!                 AC_can * delt_road               + AC_top
!
!
! where H_traffic and LE_traffic are scaled to road area.
!
!
! The system is implicited (or semi-implicited).
!
! ZIMPL=1    ---> implicit system
! ZIMPL=0.5  ---> semi-implicit system
! ZIMPL=0    ---> explicit system
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
!!                  21/11/01 (V. Masson and A. Lemonsu) bug of latent flux
!!                           for very strong evaporation (all reservoir emptied
!!                           in one time-step)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT
!
USE MODE_THERMOS
!
USE MODI_URBAN_LW_COEF
USE MODI_TRIDIAG_GROUND
!
implicit none
!!!#include <arch_specific.hf>
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
REAL, DIMENSION(:), INTENT(OUT)   :: PLW_S_TO_R   ! LW contribution from
REAL, DIMENSION(:), INTENT(OUT)   :: PLW_S_TO_W   ! sky to road and wall
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=0.5      ! implicit coefficient
REAL :: ZEXPL=0.5      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROAD,2)+SIZE(PT_WALL,2)) ::  ZA,& ! lower diag.
                                                               ZB,& ! main  diag.
                                                               ZC,& ! upper diag.
                                                               ZY,& ! r.h.s.
                                                               ZX   ! solution

REAL, DIMENSION(SIZE(PQ_CANYON)) :: zqsat

!
REAL, DIMENSION(SIZE(PTA)) :: ZDN ! snow-covered surface fraction
REAL, DIMENSION(SIZE(PTA)) :: ZDF ! snow-free surface fraction
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_W  ! L.W. interactions
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_W  ! from first Temp.
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_R  ! on second Temp.
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_R  !
REAL, DIMENSION(SIZE(PTA)) :: ZLW_N_TO_W  ! idem. but from
REAL, DIMENSION(SIZE(PTA)) :: ZLW_N_TO_R  ! snow rad.
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT_ROAD ! dq_sat/dTs
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_AC_W   ! rho * conductance (for walls)
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_R  ! rho * conductance
!                                         !     * snow-free f.
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_R_WAT ! rho * conductance for water
!                                         !     * snow-free f.
REAL, DIMENSION(SIZE(PTA)) :: ZSAC_T      ! weighted sum
!                                         ! of conductances
!                                         ! for PT_CANYON
REAL, DIMENSION(SIZE(PTA)) :: ZSAC_Q      ! weighted sum
!                                         ! of conductances
!                                         ! for PQ_CANYON
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROAD,2)) :: ZMTC_O_D_ROAD
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PTA),SIZE(PT_WALL,2)) :: ZMTC_O_D_WALL
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROAD,2)) :: ZHC_D_ROAD
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA),SIZE(PT_WALL,2)) :: ZHC_D_WALL
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROAD
! road surface temperature
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WALL
! wall surface temperature
!
INTEGER :: IROAD_LAYER           ! number of road layers
INTEGER :: IWALL_LAYER           ! number of wall layers
INTEGER :: ILAYER                ! current layer
INTEGER :: JLAYER                ! loop counter
!-------------------------------------------------------------------------------
!
!
!*      1.     Layer thermal properties
!              ------------------------
!
!*      1.1    Roads
!              -----
!
IROAD_LAYER = SIZE(PT_ROAD,2)
ZMTC_O_D_ROAD(:,:) = 0.
!
DO JLAYER=1,IROAD_LAYER-1
  ZMTC_O_D_ROAD(:,JLAYER) = 2./(  PD_ROAD(:,JLAYER  )/PTC_ROAD(:,JLAYER  ) &
                                + PD_ROAD(:,JLAYER+1)/PTC_ROAD(:,JLAYER+1) )
  ZHC_D_ROAD   (:,JLAYER) = PHC_ROAD(:,JLAYER) * PD_ROAD (:,JLAYER)
END DO
!
ZMTC_O_D_ROAD(:,IROAD_LAYER) = 2. * PTC_ROAD(:,IROAD_LAYER) &
                                    / PD_ROAD (:,IROAD_LAYER)
ZHC_D_ROAD   (:,IROAD_LAYER) = PHC_ROAD(:,IROAD_LAYER) &
                               * PD_ROAD (:,IROAD_LAYER)
!
!*      1.2    Walls
!              -----
!
IWALL_LAYER = SIZE(PT_WALL,2)
ZMTC_O_D_WALL(:,:) = 0.
!
DO JLAYER=1,IWALL_LAYER-1
  ZMTC_O_D_WALL(:,JLAYER) = 2./(  PD_WALL(:,JLAYER  )/PTC_WALL(:,JLAYER  ) &
                                + PD_WALL(:,JLAYER+1)/PTC_WALL(:,JLAYER+1) )
  ZHC_D_WALL   (:,JLAYER) = PHC_WALL(:,JLAYER) * PD_WALL (:,JLAYER)
END DO
!
ZMTC_O_D_WALL(:,IWALL_LAYER) = 2. * PTC_WALL(:,IWALL_LAYER) &
                                  / PD_WALL (:,IWALL_LAYER)
ZMTC_O_D_WALL(:,IWALL_LAYER) = 1./(  1./ZMTC_O_D_WALL(:,IWALL_LAYER)    &
                                   + 1./(XCPD*PRHOA(:)*PAC_BLD(:))      )
!
ZHC_D_WALL   (:,IWALL_LAYER) = PHC_WALL(:,IWALL_LAYER) &
                             * PD_WALL (:,IWALL_LAYER)
!
!-------------------------------------------------------------------------------
!
!*      2.    Preliminaries
!             -------------
!
!*      2.1    snow-free surface fraction
!              --------------------------
!
ZDN=PDN_ROAD(:)
ZDF(:)=1.-ZDN
!
!*      2.2    flux properties
!              ---------------
!
ZSAC_T (:) =   ZDF(:)          * PAC_ROAD(:) &
             + PWALL_O_ROAD(:) * PAC_WALL(:) &
             +                   PAC_TOP (:)
ZSAC_Q (:) =   ZDF(:) * PDELT_ROAD(:)    * PAC_ROAD_WAT(:) &
             +                             PAC_TOP(:)
ZRHO_AC_W (:) = PRHOA(:) * PAC_WALL(:)
ZRHO_ACF_R(:) = PRHOA(:) * PAC_ROAD(:) * ZDF(:)
ZRHO_ACF_R_WAT(:) = PRHOA(:) * PAC_ROAD_WAT(:) * ZDF(:)
!
!*      2.3    Surface temperatures
!              --------------------
!
ZTS_ROAD(:) = PT_ROAD(:,1)
ZTS_WALL(:) = PT_WALL(:,1)
!
!
!*      2.4    qsat, dqsat/dTs, and humidity for roads
!              ---------------------------------------
!
ZDQSAT_ROAD(:) = DQSAT(ZTS_ROAD(:),PPS(:),PQSAT_ROAD(:))
!
!
!-------------------------------------------------------------------------------
!
!*      3.     LW properties
!              -------------
!
CALL URBAN_LW_COEF(PEMIS_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,   &
                   ZDN, ZDF, PESNOW_ROAD,                          &
                   ZLW_W_TO_W, ZLW_R_TO_W, ZLW_W_TO_R, ZLW_R_TO_R, &
                   PLW_S_TO_W, PLW_S_TO_R, ZLW_N_TO_W, ZLW_N_TO_R  )

!
!-------------------------------------------------------------------------------
!
!*      4.    Inside wall layer coefficients
!             ------------------------------
!
ILAYER=1
!
ZA(:,ILAYER) =   0.

ZB(:,ILAYER) =   ZHC_D_WALL(:,IWALL_LAYER) / PTSTEP                        &
     + ZIMPL * (   ZMTC_O_D_WALL(:,IWALL_LAYER  )                          &
                 + ZMTC_O_D_WALL(:,IWALL_LAYER-1)                          &
               )

ZC(:,ILAYER) =                                                             &
     + ZIMPL * ( - ZMTC_O_D_WALL(:,IWALL_LAYER-1)                          &
               )
!
ZY(:,ILAYER) =   ZHC_D_WALL(:,IWALL_LAYER) / PTSTEP                        &
                                             * PT_WALL(:,IWALL_LAYER)      &
                 + ZMTC_O_D_WALL(:,IWALL_LAYER) * PTI_BLD(:)               &
     + ZEXPL * ( - ZMTC_O_D_WALL(:,IWALL_LAYER  )                          &
                       * PT_WALL(:,IWALL_LAYER  )                          &
                 - ZMTC_O_D_WALL(:,IWALL_LAYER-1)                          &
                       * PT_WALL(:,IWALL_LAYER  )                          &
                 + ZMTC_O_D_WALL(:,IWALL_LAYER-1)                          &
                       * PT_WALL(:,IWALL_LAYER-1)                          &
               )
!
!-------------------------------------------------------------------------------
!
!*      5.     Other wall layers coefficients
!              ------------------------------
!
DO JLAYER=2,IWALL_LAYER-1

  ILAYER=IWALL_LAYER-JLAYER+1

  ZA(:,ILAYER) =                                                         &
         ZIMPL * ( - ZMTC_O_D_WALL(:,JLAYER  )                           &
                 )

  ZB(:,ILAYER) =   ZHC_D_WALL(:,JLAYER)/PTSTEP                           &
       + ZIMPL * (   ZMTC_O_D_WALL(:,JLAYER  )                           &
                   + ZMTC_O_D_WALL(:,JLAYER-1)                           &
                 )

  ZC(:,ILAYER) =                                                         &
         ZIMPL * ( - ZMTC_O_D_WALL(:,JLAYER-1)                           &
                 )
!
  ZY(:,ILAYER) =   ZHC_D_WALL(:,JLAYER)/PTSTEP * PT_WALL(:,JLAYER)     &
        + ZEXPL * (    ZMTC_O_D_WALL(:,JLAYER  ) * PT_WALL(:,JLAYER+1) &
                     - ZMTC_O_D_WALL(:,JLAYER  ) * PT_WALL(:,JLAYER  ) &
                     - ZMTC_O_D_WALL(:,JLAYER-1) * PT_WALL(:,JLAYER  ) &
                     + ZMTC_O_D_WALL(:,JLAYER-1) * PT_WALL(:,JLAYER-1) &
                  )
END DO
!
!-------------------------------------------------------------------------------
!
!*      6.     Surface wall layer coefficients
!              -------------------------------
!
ILAYER=IWALL_LAYER
!
!
ZA(:,ILAYER) =                                                             &
       ZIMPL * ( - ZMTC_O_D_WALL(:,1)                                      &
               )

ZB(:,ILAYER) =   ZHC_D_WALL(:,1)/PTSTEP                                    &
     + ZIMPL * ( -          4. *ZTS_WALL(:)**3 * ZLW_W_TO_W(:)             &
                 + ZRHO_AC_W(:) * XCPD                                     &
                           * (1.-PAC_WALL(:)*PWALL_O_ROAD(:)/ZSAC_T (:))   &
                 + ZMTC_O_D_WALL(:,1)                                      &
               )

ZC(:,ILAYER) =                                                             &
       ZIMPL * ( -          4.*ZTS_ROAD(:)**3 * ZLW_R_TO_W(:)              &
                 - ZRHO_AC_W(:) * XCPD * PAC_ROAD(:)*ZDF(:)/ZSAC_T (:)     &
               )


ZY(:,ILAYER) =   ZHC_D_WALL(:,1)/PTSTEP*ZTS_WALL(:)                        &
                 + PABS_SW_WALL(:)                                         &
                 + PLW_RAD    (:)    * PLW_S_TO_W(:)                       &
                 + ZTS_WALL   (:)**4 * ZLW_W_TO_W(:)                       &
                 + ZTS_ROAD   (:)**4 * ZLW_R_TO_W(:)                       &
                 + PTSNOW_ROAD(:)**4 * ZLW_N_TO_W(:)                       &
                 + ZRHO_AC_W(:) * XCPD * PTA(:)                            &
                                       * PAC_TOP(:) / ZSAC_T (:)           &
                 + PAC_WALL(:) / ZSAC_T (:)                                &
                                * (   PH_TRAFFIC (:) / (1.-PBLD(:))        &
                                    + PHSNOW_ROAD(:) * ZDN(:)       )      &
     + ZIMPL * (                                                           &
                 - 4.*ZTS_WALL(:)**4 * ZLW_W_TO_W(:)                       &
                 - 4.*ZTS_ROAD(:)**4 * ZLW_R_TO_W(:)                       &
               )                                                           &
     + ZEXPL * (   ZMTC_O_D_WALL(:,1) * PT_WALL(:,2)                       &
                 - ZMTC_O_D_WALL(:,1) * PT_WALL(:,1)                       &
                 - ZRHO_AC_W(:) * XCPD * (1.-PAC_WALL(:)*PWALL_O_ROAD(:)   &
                                                         /ZSAC_T (:)     ) &
                                * ZTS_WALL(:)                              &
                 + ZRHO_AC_W(:) * XCPD * PAC_ROAD(:)*ZDF(:)/ZSAC_T(:)      &
                                * ZTS_ROAD(:)                              &
               )
!
!-------------------------------------------------------------------------------
!
!*      7.     Surface road layer coefficients
!              -------------------------------
!
ILAYER=IWALL_LAYER+1
!
!
ZA(:,ILAYER) =                                                           &
       ZIMPL * ( - ZDF(:) * 4.*ZTS_WALL(:)**3 * ZLW_W_TO_R(:)            &
                 - ZRHO_ACF_R(:) * XCPD * PAC_WALL(:)                    &
                         * PWALL_O_ROAD(:) / ZSAC_T (:)                  &
               )

ZB(:,ILAYER) =   ZHC_D_ROAD(:,1)/PTSTEP                                         &
     + ZIMPL * ( - ZDF(:) * 4.*ZTS_ROAD(:)**3 * ZLW_R_TO_R(:)                   &
                 + ZRHO_ACF_R(:) * XCPD * (1.-ZDF(:)*PAC_ROAD(:)/ZSAC_T(:))     &
                 + ZRHO_ACF_R_WAT(:) * XLVTT * PDELT_ROAD(:) * ZDQSAT_ROAD(:)   &
                        * (1.-ZDF(:)*PAC_ROAD_WAT(:)*PDELT_ROAD(:)/ZSAC_Q(:))   &
                 + ZMTC_O_D_ROAD(:,1)                                           &
               )

ZC(:,ILAYER) =                                                           &
       ZIMPL * ( - ZMTC_O_D_ROAD(:,1)                                    &
                 )

ZY(:,ILAYER) =   ZHC_D_ROAD(:,1)/PTSTEP*PT_ROAD(:,1)                     &
                 + ZDF(:) * PABS_SW_ROAD(:)                              &
                 + ZDF(:) * PLW_RAD    (:)    * PLW_S_TO_R(:)            &
                 + ZDF(:) * ZTS_WALL   (:)**4 * ZLW_W_TO_R(:)            &
                 + ZDF(:) * ZTS_ROAD   (:)**4 * ZLW_R_TO_R(:)            &
                 + ZDF(:) * PTSNOW_ROAD(:)**4 * ZLW_N_TO_R(:)            &
                 + ZRHO_ACF_R(:) * XCPD * PTA(:)                         &
                                        * PAC_TOP(:) / ZSAC_T (:)        &
                 + ZDF(:) * PAC_ROAD(:) / ZSAC_T (:)                     &
                            * (   PH_TRAFFIC (:) / (1.-PBLD(:))          &
                                + PHSNOW_ROAD(:) * ZDN(:)       )        &
                 - ZRHO_ACF_R_WAT(:) * XLVTT * PDELT_ROAD(:)             &
                         * (  PQSAT_ROAD(:)                              &
                             -(   PQSAT_ROAD     (:) *PDELT_ROAD(:)      &
                                                  *ZDF(:)*PAC_ROAD_WAT(:)&
                                + PQA            (:) *     PAC_TOP(:)    &
                               ) / ZSAC_Q(:)                             &
                           )                                             &
                 + ZDF(:) * PAC_ROAD_WAT(:) * PDELT_ROAD(:) / ZSAC_Q (:) &
                            * (   PLE_TRAFFIC (:) / (1.-PBLD(:))         &
                                + PLESNOW_ROAD(:) * ZDN(:)       )       &
                 + ZDN(:) * PGSNOW_ROAD(:)                               &
     + ZIMPL * ( - ZDF(:) * 4.*ZTS_WALL(:)**4 * ZLW_W_TO_R(:)            &
                 - ZDF(:) * 4.*ZTS_ROAD(:)**4 * ZLW_R_TO_R(:)            &
                 + ZRHO_ACF_R_WAT(:) * XLVTT * PDELT_ROAD(:)             &
                               * (1.-PDELT_ROAD(:)*PAC_ROAD_WAT(:)*ZDF(:)&
                                                    /ZSAC_Q (:)        ) &
                               * ZDQSAT_ROAD(:) * ZTS_ROAD(:)            &
               )                                                         &
     + ZEXPL * (   ZRHO_ACF_R(:) * XCPD * PAC_WALL(:) * PWALL_O_ROAD(:)  &
                                 * ZTS_WALL(:) / ZSAC_T (:)              &
                 - ZRHO_ACF_R(:) * XCPD * ZTS_ROAD(:)                    &
                                 * ( 1. - PAC_ROAD(:)*ZDF(:)/ZSAC_T(:) ) &
                 - ZMTC_O_D_ROAD(:,1) * PT_ROAD(:,1)                     &
                 + ZMTC_O_D_ROAD(:,1) * PT_ROAD(:,2)                     &
               )
!
!-------------------------------------------------------------------------------
!
!*      8.     Other road layers coefficients
!              ------------------------------
!
DO JLAYER=2,IROAD_LAYER-1

  ILAYER=IWALL_LAYER+JLAYER

  ZA(:,ILAYER) =                                                       &
         ZIMPL * ( - ZMTC_O_D_ROAD(:,JLAYER-1)                         &
                 )

  ZB(:,ILAYER) =   ZHC_D_ROAD(:,JLAYER)/PTSTEP                         &
       + ZIMPL * (   ZMTC_O_D_ROAD(:,JLAYER-1)                         &
                   + ZMTC_O_D_ROAD(:,JLAYER  )                         &
                 )

  ZC(:,ILAYER) =                                                       &
         ZIMPL * ( - ZMTC_O_D_ROAD(:,JLAYER  )                         &
                 )
!
  ZY(:,ILAYER) =   ZHC_D_ROAD(:,JLAYER)/PTSTEP*PT_ROAD(:,JLAYER)       &
       + ZEXPL * (   ZMTC_O_D_ROAD(:,JLAYER-1) * PT_ROAD(:,JLAYER-1)   &
                   - ZMTC_O_D_ROAD(:,JLAYER-1) * PT_ROAD(:,JLAYER  )   &
                   - ZMTC_O_D_ROAD(:,JLAYER  ) * PT_ROAD(:,JLAYER  )   &
                   + ZMTC_O_D_ROAD(:,JLAYER  ) * PT_ROAD(:,JLAYER+1)   &
                 )
END DO
!
!-------------------------------------------------------------------------------
!
!*      9.     Inside road layer coefficients
!              ------------------------------
!
ILAYER=IWALL_LAYER+IROAD_LAYER
!
ZA(:,ILAYER) =                                                        &
       ZIMPL * ( - ZMTC_O_D_ROAD(:,IROAD_LAYER-1)                     &
                 )

ZB(:,ILAYER) =   ZHC_D_ROAD(:,IROAD_LAYER) / PTSTEP                   &
     + ZIMPL * (   ZMTC_O_D_ROAD(:,IROAD_LAYER-1)                     &
               )

ZC(:,ILAYER) =   0.
!
ZY(:,ILAYER) =   ZHC_D_ROAD(:,IROAD_LAYER) / PTSTEP                   &
                                             * PT_ROAD(:,IROAD_LAYER) &
     + ZEXPL * (   ZMTC_O_D_ROAD(:,IROAD_LAYER-1)                     &
                       * PT_ROAD(:,IROAD_LAYER-1)                     &
                 - ZMTC_O_D_ROAD(:,IROAD_LAYER-1)                     &
                       * PT_ROAD(:,IROAD_LAYER  )                     &
               )
!
!-------------------------------------------------------------------------------
!
!*     10.     Tri-diagonal system resolution
!              ------------------------------
!
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
DO JLAYER=1,IWALL_LAYER
  ILAYER=IWALL_LAYER-JLAYER+1
  PT_WALL(:,JLAYER) = ZX(:,ILAYER)
END DO
!
DO JLAYER=1,IROAD_LAYER
  ILAYER=IWALL_LAYER+JLAYER
  PT_ROAD(:,JLAYER) = ZX(:,ILAYER)
END DO
!
!-------------------------------------------------------------------------------
!
!*     12.    Road and wall absorbed infra-red radiation on snow-free surfaces
!             ----------------------------------------------------------------
!
PABS_LW_ROAD(:) =  PLW_S_TO_R(:)*PLW_RAD    (:)             &
                 + ZLW_R_TO_R(:)*PT_ROAD    (:,1)**4        &
                 + ZLW_W_TO_R(:)*PT_WALL    (:,1)**4        &
                 + ZLW_N_TO_R(:)*PTSNOW_ROAD(:)**4
!
PABS_LW_WALL(:) =  PLW_S_TO_W(:)*PLW_RAD    (:)             &
                 + ZLW_W_TO_W(:)*PT_WALL    (:,1)**4        &
                 + ZLW_R_TO_W(:)*PT_ROAD    (:,1)**4        &
                 + ZLW_N_TO_W(:)*PTSNOW_ROAD(:)**4
!
!-------------------------------------------------------------------------------
!
!*     13.    Air canyon temperature at time t+dt
!             -----------------------------------
!
PT_CANYON(:) = (  PT_ROAD    (:,1) * PAC_ROAD(:) * ZDF(:)             &
                + PT_WALL    (:,1) * PAC_WALL(:) * PWALL_O_ROAD(:)    &
                + PTA        (:)   * PAC_TOP (:)                      &
                + PH_TRAFFIC (:)   / (1.-PBLD(:)) / PRHOA(:) / XCPD   &
                + PHSNOW_ROAD(:)   * ZDN(:)       / PRHOA(:) / XCPD  )&
               / ZSAC_T(:)
!
!-------------------------------------------------------------------------------
!
!*     14.     canyon air specific humidities
!              ------------------------------
!
!*     14.1    New saturated specified humidity near the road surface
!              ------------------------------------------------------
!
!
PQSAT_ROAD(:) =  QSAT(PT_ROAD(:,1),PPS(:))
!
!*     14.2    Canyon air specific humidity
!              ----------------------------
!
PQ_CANYON(:) = (  PQSAT_ROAD     (:) * PAC_ROAD_WAT(:) * ZDF(:)      &
                                                  * PDELT_ROAD(:)    &
                + PQA            (:) * PAC_TOP(:)                    &
                + PLE_TRAFFIC (:) / (1.-PBLD(:)) / PRHOA(:) / XLVTT  &
                + PLESNOW_ROAD(:) * ZDN(:)       / PRHOA(:) / XLVTT )&
               / ZSAC_Q(:)


!
!*     14.3    Resolve instabilities developing for winter conditions
!              ----------------------------
!
! If canyon specific humidity gets too small or even negative
! set it to specific humidity at the lowest level (KW)
WHERE(PQ_CANYON(:) <= 1.e-6) PQ_CANYON(:) = PQA(:)

! If canyon specific humidity gets larger than saturation
! clip to saturation (KW)
zqsat(:) = QSAT(PT_CANYON(:),PPS(:))
WHERE(PQ_CANYON(:) > zqsat(:)) PQ_CANYON(:) = zqsat(:)

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ROAD_WALL_LAYER_E_BUDGET
