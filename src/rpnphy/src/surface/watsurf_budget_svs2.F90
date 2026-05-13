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
module watsurf_budget_svs2_mod
  implicit none
  public
contains
SUBROUTINE WATSURF_BUDGET_SVS2 ( DT, ESUBSNC, SUBSNC_CUM, &
     EG, EGV, ER_VL,ER_VH, ETR_VL,ETR_VH, RR, RR_VEG, RSNOW, RSNOWV, &
     WTA, WTG, ACROOT, WRMAX_VL, WRMAX_VH, WRMAX_FL, &
     SNM, SVM, SNCMA, WR_VL, WR_VH, WFL, WR_VLT, WR_VHT, WFLT, &
     PG, ETR_GRID, EG_GRID, N)
  !
  use sfc_options
  use svs_configs
  use svs2_tile_configs

  implicit none
!!!#include <arch_specific.hf>

  !
  INTEGER N,K

  REAL DT

  ! input
  real, dimension(n)        :: eg, egv, er_vl, etr_vl, er_vh, etr_vh, rr, rr_veg
  real, dimension(n)        :: esubsnc, subsnc_cum, sncma
  real, dimension(n,svs2_tilesp1) :: wta, wtg
  real, dimension(n,nl_svs) :: acroot
  real, dimension(n)        :: rsnow, rsnowv, wrmax_vl, wrmax_vh, wrmax_fl, snm, svm
  ! prognostic vars (I/0)
  real, dimension(n)        :: wr_vl, wr_vh, wfl, wr_vlt, wr_vht, wflt
  ! output
  real, dimension(n)        ::  PG, EG_GRID
  real, dimension(n,nl_svs) :: ETR_GRID

  !
  !Author
  !          N.Leroux, V.Vionnet (May 2024)
  !Revisions
  !
  ! 001      Includes forest litter layer - N.Leroux (Dec 2025)
  !Object
  !     Calculates the evolution of the intercepted liquid water 
  !     in the vegetation canopy (Wr).
  !     Also determines the water available to infiltrate soil and evapotranspiration 
  !     fluxes from the soil
  !
  !Arguments
  !
  !          - INPUT -
  !
  ! DT       timestep
  !
  !          --- Precipitation and Evaporation Rates ---
  !
  ! EG             evaporation rate (no fraction) over bare ground (1st soil layer) [kg/m2/s]
  ! EGV            evaporation rate (no fraction) over ground below high vegetation [kg/m2/s]
  ! ER_VL          direct evaporation rate (no fraction) from low veg. leaves [kg/m2/s]
  ! ETR_VL         evapotranspiration rate from low vegetation (no fraction) [kg/m2/s]
  ! ER_VH          direct evaporation rate (no fraction) from high veg. leaves [kg/m2/s]
  ! ETR_VH         evapotranspiration rate from high vegetation (no fraction) [kg/m2/s]
  ! RR             rain rate at the surfaceover bare ground/low veg [kg/m2/s]
  ! RR_VEG         rain rate at the surface below high veg, calculated in snow_interception_svs2 [kg/m2/s]
  ! RSNOW          flow of liquid water through the bare ground/low veg. snow pack - to the soil [kg/m2/s]
  ! RSNOWV         flow of liquid water through the snow-under-veg pack - to the soil [kg/m2/s]
  ! ESUBSNC      sublimation rate from incercepted snow calculated in ebudget [kg m-2 s-1]
  !
  !          --- (Surface) Cover Fraction  ---
  !
  ! WTA            weights for SVS2 surface tiles as seen from SPACE [0-1]
  ! WTG            weights for SVS2 surface tiles as seen from GROUND [0-1]
  !
  !          --- Vegetation characteristics ---
  !
  ! ACROOT(NL_SVS) active root fraction in each soil layer [0-1]
  !                (ETR is multiplied by ACROOT to get ETR for each layer)
  ! WRMAX_VL         max water content retained on low vegetation [kg/m2]
  ! WRMAX_VH         max water content retained on high vegetation [kg/m2]
  ! WRMAX_FL         max water content retained in forest litter layer [kg/m2]
  !
  !
  !          --- Water and snow in the different reservoirs and tiles ---
  !
  ! SNM           snow water equivalent (SWE) for bare ground and low veg snow [kg/m2]
  ! SVM           snow water equivalent (SWE) for snow under high veg [kg/m2]
  ! WR_VL          water content retained by low vegetation canopy [kg/m2]
  ! WR_VH          water content retained by high vegetation canopy [kg/m2]
  ! WFL          water content retained by forest litter layer [kg/m2]
  ! SUBSNC_CUM   Cumulated mass loss due to sublimation of incercepted snow [kg m-2]
  ! SNCMA       mass of intercepted snow in the canopy [kg/m2]
  !
  !          -  OUTPUT  -
  !
  !
  !          ---  Soil water fluxes and Runoff ---
  !
  ! ETR_GRID(NL_SVS) Evapotranspiration from each layer (grid box average) [m/s]
  ! EG_GRID(NL_SVS)  evaporation rate over bare ground and bare ground below high veg (grid box average) [kg/m2/s]
  ! PG         Water available for infiltration into soil [kg/m2/s]
  ! WR_VLT         new water content retained by low vegetation canopy [kg/m2]
  ! WR_VHT         new water content retained by high wvegetation canopy [kg/m2]
  ! WFLT           new water content retained by forest litter layer [kg/m2]
  !
  !          -  DIMENSIONS  -
  !
  ! N              number of grid cells
  !
  !
  INTEGER I
  real, dimension(n)        :: rveg_vh, rveg_vl, rveg_fl

  !
  !
  !-------------------------------------
  !   1.        EVOLUTION OF THE INTERCEPTED SNOW ON THE HIGH VEG
  !
  !                                  Remove sublimation from snow canopy reservoir
  !                                  Only applicable when the option ABV is applied
  IF (CANO_REF_FORCING .EQ.'ABV' .AND. LSNOW_INTERCEPTION_SVS2 ) THEN
      DO I=1,N
         SNCMA(I) = SNCMA(I) -  ESUBSNC(I)*DT
         SUBSNC_CUM(I) = SUBSNC_CUM(I) + ESUBSNC(I) * DT
      END DO
   ENDIF

  !
  !
  !-------------------------------------
  !   1.        EVOLUTION OF THE EQUIVALENT WATER CONTENT WR_VL, WR_VH, and WFL
  !
  !                                  Remove evaporation from canopy reservoir
  !                                  Then add precipitation. Careful: rain falls
  !                                  on canopy only if there is no snow under the
  !                                  vegetation. If there is snow on the ground,
  !                                  all rain goes directly to the snowpack and no
  !                                  liquid water remains in the vegetation
  !
   DO I=1,N
      !
      !   Mass balance of the high vegetation
      !
      IF( WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN

         !                                  Remove P-E from liquid water in vegetation
         !                                  Rain trapped on high vegetation is computed in snow_interception
         !                                  Sanity check: WR_VHT must be positive
         !                                  since ER is limited to WR/DT+RR in ebudget
         WR_VHT(I) = MAX(0., WR_VH(I) - DT * ER_VH(I))
         !                                  Compute canopy drip
         !                                  if Wr > Wrmax, there is runoff
         RVEG_VH(I) = MAX(0., (WR_VHT(I) - WRMAX_VH(I)) / DT )
         !
         !                                  We must be smaller than Wrmax
         !                                  after runoff is removed
         !
         WR_VHT(I) = MIN(WR_VHT(I), WRMAX_VH(I))

      ELSE
         WR_VHT(I) = 0.0
         RVEG_VH(I) = 0.
      ENDIF
      !
      !   Mass balance of the low vegetation
      !
      IF( WTG(I,indx_svs2_vl) .GE.EPSILON_SVS) THEN

         IF (SNM(I).GE.CRITSNOWMASS .OR. RSNOW(I).GE. EPSILON_SVS )THEN
            !                                  There is snow on the ground, remove evaporation
            !                                  then drain all liquid water from the vegetation
            WR_VLT(I) = WR_VL(I) - DT * ER_VL(I)
            !                                  Liquid water in vegetation becomes runoff
            RVEG_VL(I) = MAX(0.0,WR_VLT(I)/DT) ! Runoff mult. by dt later on, so water balance maintained..
            WR_VLT(I) = 0.0

         ELSE
            !                                  Remove P-E from liquid water in vegetation and rain can be trapped by low veg
            !                                  Sanity check: WR_VHT must be positive
            !                                  since ER is limited to WR/DT+RR in ebudget
            WR_VLT(I) = MAX(0., WR_VL(I) - DT * (ER_VL(I) -  RR(I)))
            !                                  Compute canopy drip
            !                                  if Wr > Wrmax, there is runoff
            RVEG_VL(I) = MAX(0., (WR_VLT(I) - WRMAX_VL(I)) / DT )
            !
            !                                  Wr must be smaller than Wrmax
            !                                  after runoff is removed
            !
            WR_VLT(I) = MIN(WR_VLT(I), WRMAX_VL(I))
         END IF
      ELSE
         WR_VLT(I) = 0.0
         RVEG_VL(I) = 0.
      ENDIF
      !
      !   Mass balance of the forest litter layer
      !
      IF (LFORLIT .AND. WTG(I,indx_svs2_vh) .GE. EPSILON_SVS) THEN 

         IF (SVM(I).GE.CRITSNOWMASS .OR. RSNOWV(I).GE. EPSILON_SVS ) THEN
            !                  There is snow on the ground or there was no during this time step (fraction of snow in vegh > 0)
            WFLT(I) = MAX(0., WFL(I) - DT * (WTG(I,indx_svs2_gv) / WTG(I,indx_svs2_vh) * EGV(I) - RSNOWV(I))) 
         ELSE
            !                  There is no snow, so rain becomes an input (fraction of snow in vegh = 0)
            WFLT(I) = MAX(0., WFL(I) - DT * (EGV(I) -  RR_VEG(I)))
         ENDIF
         !                                  Liquid water in forest litter becomes runoff
         RVEG_FL(I) =  MAX(0., (WFLT(I) - WRMAX_FL(I)) / DT )

         WFLT(I) = MIN(WFLT(I), WRMAX_FL(I))
         EGV(I) = 0. ! Prevents evaporation from the ground in ebudget_svs2
      ELSE
         RVEG_FL(I) = 0.
         WFLT(I) = 0.
      ENDIF

   END DO

  !
  !        2.     CALCULATE PRECIPITATION and VEGETATION+SNOW RUNOFF REACHING THE GROUND
  !               ------------------------------------------

  !
  !                                  Have vegetation runoff only if vegetation present
  !                                  Precipitation reaches ground directly only if no snow
  !                                  otherwise goes to snowpack, and reaches ground as snow runoff
  !
  DO I=1,N
      IF( (WTG(I,indx_svs2_vh)+ WTG(I,indx_svs2_vl)) .GE.EPSILON_SVS) THEN

         IF (.not. LFORLIT) THEN ! No forest litter
            ! There is vegetation -- first add snowpack runoff and runoff from vegetation
            PG(I) = (1. - WTG(I,indx_svs2_vh)) * RSNOW(I) + WTG(I,indx_svs2_vh) * RSNOWV(I) + &
                  WTG(I,indx_svs2_vl)*RVEG_VL(I) + WTG(I,indx_svs2_vh) * RVEG_VH(I)
         ELSE
            ! There is vegetation -- first add snowpack runoff and runoff from vegetation
            PG(I) = WTG(I,indx_svs2_sn) * RSNOW(I) + WTG(I,indx_svs2_vh) * RVEG_FL(I) + &
               WTG(I,indx_svs2_vl)*RVEG_VL(I)
         ENDIF

         IF( SNM(I).LT.CRITSNOWMASS .AND. RSNOW(I) .LT. EPSILON_SVS ) THEN
            ! No snow on bare ground and low veg, rain falls to bare ground
            PG(I) = PG(I) + WTG(I,indx_svs2_bg) * RR(I)
         ENDIF

         IF (.not. LFORLIT) THEN ! No forest litter
            IF ( SVM(I).LT.CRITSNOWMASS .AND. RSNOWV(I) .LT. EPSILON_SVS) THEN
               ! No snow at the end AND during the time step  below high veg, rain below canopy falls to ground below canopy
               PG(I) = PG(I) +  WTG(I,indx_svs2_vh) * RR_VEG(I)
            ENDIF
         ENDIF

       ELSE
         ! bare ground only
         PG(I) = RSNOW(I)
         IF ( SNM(I).LT.CRITSNOWMASS .AND. RSNOW(I) .LT. EPSILON_SVS) THEN
            ! have 100% bare ground and no snow at the end AND during the time step, all rain reaches surface
            PG(I) = PG(I) + RR(I)
         ENDIF
      ENDIF

  ENDDO




  !        3.      CALCULATE GRID AVERAGE EVAPOTRANSPIRATION
  !              -----------------------------------------
  !
  DO I=1,N
     DO K=1,NL_SVS
        ! Initial line in hydro_svs
        !ETR_GRID(I,K)=((VEGL(I)*(1.-PSN(I))+VEGH(I)*(1.-PSNVH(I)))*ACROOT(I,K)*ETR(I))/(1000.*DELZ(K))
        ! ETR_GRID(I,K)=((VEGL(I)*(1.-PSN(I))*ETR_VL(I)+VEGH(I)*(1.-PSNVH(I))*ETR_VH(I))*ACROOT(I,K))/(1000.*DELZ(K))
        ! weights need to be consistent with those in ebudget, otherwise, mass balance is not closed
        ETR_GRID(I,K) = ((WTG(I,indx_svs2_vl)*ETR_VL(I) +  WTG(I,indx_svs2_vh)*ETR_VH(I))*ACROOT(I,K))/(1000.*DELZ(K))
     END DO
     EG_GRID(I) = WTG(I,indx_svs2_bg) * EG(I)  + WTG(I,indx_svs2_gv) * EGV(I)
  END DO




  RETURN
END SUBROUTINE WATSURF_BUDGET_SVS2
end module watsurf_budget_svs2_mod
