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
SUBROUTINE HYDRO_SVS ( DT, &
     EG, ER, ETR, RR, RSNOW, RSNOWV, &
     IMPERVU, VEGL, VEGH, PSN, PSNVH, ACROOT, WRMAX, &
     WSAT, KSAT, PSISAT, BCOEF, FBCOF, WFCINT, GRKEF, &
     SNM, SVM, WR, WRT, WD, WDT, WF, WFT, &
     KSATC, KHC, PSI, GRKSAT, WFCDP, &
     F, LATFLW, RUNOFF, N,  WATPND, MAXPND)
  !
  use sfc_options
  use svs_configs
  implicit none
!!!#include <arch_specific.hf>

  !     
  INTEGER N,K 
  ! CONSTANTS for horizontal decay of GRKSAT
  ! 
  REAL, PARAMETER :: GRKSAT_C1=10.0
  REAL, PARAMETER :: GRKSAT_C2=5.0


  REAL DT, W
  
  INTEGER KFICE ! Option for the correction factor for hydraulic conductivity
  !    0: Correction factor taken from Zhang and Gray (1997). Same as CLASS 3.6
  !    1: Impedance factor taken from SURFEX (Boone et al., 2000)
  !    3: No modification of hydraulic conductivity in presence of ice 

  INTEGER WAT_REDIS ! Option for the redistribution of water in case of over-saturation after the soil_fluxes solver
  !    0: Default param:
  !    1: New param:  - in case of vertical redistribution, liquid water that flows to the next layer is limited by 
  !                       the effective porosity. Any water in excess is added to the lateral flow
  !                   - harmonic mean of the hydraulic conductivity  

  ! input
  real, dimension(n)        :: eg, er, etr, rr, impervu
  real, dimension(n)        :: psn, psnvh, vegh, vegl
  real, dimension(n,nl_svs) :: bcoef, fbcof, acroot, ksat
  real, dimension(n,nl_svs) :: psisat, wfcint, wsat
  real, dimension(n)        :: grkef, rsnow, rsnowv, wrmax, snm, svm 
  ! prognostic vars (I/0)
  real, dimension(n)        :: wr, wrt
  real, dimension(n,nl_svs) :: wd , wdt, wf, wft
  ! output
  real, dimension(n,nl_svs) :: grksat, ksatc
  real, dimension(n,nl_svs-1):: khc, psi
  real, dimension(n)        :: wfcdp
  real, dimension(n,nl_svs+1):: f
  real, dimension(n,nl_svs) :: latflw
  real, dimension(n)        :: runoff
  real, dimension(n)        :: watpnd, maxpnd

  !
  !Author
  !          N.Alavi, S.Zhang, E. Gaborit, V. Fortin, S. Belair, V. Vionnet et al.  (June 2015) 
  !Revisions
  !
  ! 001  
  !Object
  !     Calculates the evolution of the soil water contents
  !     liquid water retained in the vegetation canopy (Wr).  
  !     Also determine the runoff, lateral flow and drainage
  !     using interflow parametrization for sloping train
  !     (based on Soulis et al. 2000 and 2011)
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
  ! ER             direct evaporation rate (no fraction) from veg. leaves [kg/m2/s]
  ! ETR            evapotranspiration rate (no fraction) [kg/m2/s]
  ! RR             rain rate at the surface [kg/m2/s]
  ! RSNOW          flow of liquid water through the bare ground/low veg. snow pack - to the soil [kg/m2/s]
  ! RSNOWV         flow of liquid water through the snow-under-veg pack - to the soil [kg/m2/s]
  !
  !          --- (Surface) Cover Fraction  ---
  !
  ! IMPERVU        fraction of land surface that is considered impervious (related to urban areas) [0-1]
  ! VEGL           fraction of LOW vegetation [0-1]
  ! VEGH           fraction of HIGH vegetation [0-1]
  ! PSN            fraction of bare ground or low veg. covered by snow [0-1]
  ! PSNVH          fraction of HIGH vegetation covered by snow [0-1]
  !
  !          --- Vegetation characteristics ---
  !
  ! ACROOT(NL_SVS) active root fraction in each soil layer [0-1]
  !                (ETR is multiplied by ACROOT to get ETR for each layer)
  ! WRMAX          max water content retained on vegetation [kg/m2]
  !
  !          --- Surface characteristics ---
  ! MAXPND       maximum ponding depth [m] (used only if ponding is activated)
  !
  !          --- Soil characteristics ---
  !
  ! WSAT  (NL_SVS) volumetric water content at soil saturation per layer [m3/m3]
  ! KSAT  (NL_SVS) vertical hydraulic conductivity at saturation per layer [m/s]
  ! PSISAT(NL_SVS) value of soil water suction at air-entry (near saturation) per layer [m]
  ! BCOEF (NL_SVS) slope of the retention curve per layer
  ! FBCOF (NL_SVS) parameter derived from BCOEF per layer to determine field capacity (Soulis et al. 2011)
  ! WFCINT(NL_SVS) volumetric water content at field capacity for interflow (per layer) [m3/m3]
  ! GRKEF          ratio of tile slope to tile length (needed for watdrain to calculate lateral flow) [m-1]
  !
  !          --- Prognostic variables of SVS not modified by HYDRO_SVS ---
  !
  ! SNM           snow water equivalent (SWE) for bare ground and low veg snow [kg/m2]
  ! SVM           snow water equivalent (SWE) for snow under high veg [kg/m2]
  !
  !          - INPUT/OUTPUT  -
  !
  !          --- Prognostic variables of SVS modified by HYDRO_SVS ---
  !
  ! WR             water content retained by vegrtation canopy [kg/m2]
  ! WRT            new water content retained by vegetation canopy [kg/m2]
  ! WD (NL_SVS)    soil volumetric water content (per layer) [m3/m3]
  ! WDT(NL_SVS)    new soil volumetric water content (per layer) [m3/m3]
  ! WF (NL_SVS)    frozen soil volumetric water (per layer) [m3/m3]
  ! WFT(NL_SVS)    new frozen soil volumetric water (per layer) [m3/m3]
  ! WATPND         depth of ponded water [m] (updated only if ponding is activated)
  !
  !          -  OUTPUT  -
  !
  !          ---  Diagnostic Soil Parameters (used only for testing SVS) ---
  !
  ! KSATC (NL_SVS) KSAT adjusted for presence of ice (per layer) [m/s]
  ! KHC (NL_SVS-1) soil hydraulic conductivity at soil boundaries (of layer) [m/s]
  ! PSI (NL_SVS-1) soil water potential at the soil boundaries (of layer) [m]
  ! GRKSAT(NL_SVS) horizontal hydraulic conductivity at saturation (per layer) [m/s]
  ! WFCDP          volumetric water content at field capacity at the bottom of the soil column [m3/m3]
  !
  !          ---  Soil water fluxes and Runoff ---
  !
  ! F(NL_SVS+1)    water flux between soil layers and through bottom (in meters, converted to mm at the end ) [kg/m2/dt]
  ! LATFLW(NL_SVS) lateral flow (interflow) (in meters in calculation, converted to mm at the end) per layer [kg/m2/dt]
  ! RUNOFF         surface runoff [kg/m2/dt]
  !
  !          -  DIMENSIONS  -
  !
  ! N              number of grid cells
  !
  !
  INTEGER I

  ! local arrays

  real                        :: fice

  real, dimension(n,nl_svs)   :: delzvec
  real, dimension(n,nl_svs)   :: asatfc, wsatc, asat0, grkefl, ksatmean
  real, dimension(n,nl_svs)   :: etr_grid
  real, dimension(n)          :: asat1, basflw, satsfc, subflw , pg

  real, dimension(n)          :: wrt_vl,wrt_vh,rveg_vl,rveg_vh
  real                        :: wat_down

  real, dimension(n)          :: pond_ret ! Amount of surface runoff retained in the pond
  
  ! For the Runge-Kutta method
  real, dimension(n,nl_svs)   :: wd_rk, dwd_rk1, dwd_rk2, dwd_rk3, dwd_rk4
  real, dimension(n,nl_svs)   :: over_rk1, over_rk2, over_rk3, over_rk4
  real, dimension(n,nl_svs+1) :: f_rk

  !***********************************************************************
  !
  !   0.        INITIALIZE TO ZERO THE VERTICAL AND LATER FLUXES
  !

  DO I=1,N
     DO K=1,NL_SVS
        LATFLW(I,K)=0.
        F(I,K)=0.
     ENDDO
     F(I,NL_SVS+1)=0.
  ENDDO

  ! Option for soil freezing
  !    0: Correction factor taken from Zhang and Gray (1997). Same as CLASS 3.6
  !    1: Impedance factor taken from SURFEX (Boone et al., 2000)
  !    3: No modification of hydraulic conductivity in presence of ice 
  KFICE = 0

  !
  ! Option for the redistribution of water in case of over-saturation after the soil_fluxes solver
  !    0: Default param:
  !    1: New param:  - in case of vertical redistribution, liquid water that flows to the next layer is limited by 
  !                       the effective porosity. Any water in excess is added to the lateral flow
  !                   - harmonic mean of the hydraulic conductivity  
  IF(LSOIL_FREEZING_SVS1) THEN
       WAT_REDIS = 1
  ELSE
       WAT_REDIS = 0
  ENDIF
  !
  !-------------------------------------
  !   1.        EVOLUTION OF THE EQUVALENT WATER CONTENT Wr
  !      
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
    IF(VEGH(I).GE.EPSILON_SVS) THEN

     IF (SVM(I).GE.CRITSNOWMASS)THEN
        !                                  There is snow on the ground, remove evaporation
        !                                  then drain all liquid water from the vegetation
        WRT_VH(I) = WR(I) - DT * ER(I)
        !                                  Liquid water in vegetation becomes runoff
        RVEG_VH(I) = MAX(0.0,WRT_VH(I)/DT) ! Runoff mult. by dt later on, so water balance maintained..
        WRT_VH(I) = 0.0

     ELSE
        !                                  Remove P-E from liquid water in vegetation
        !                                  Sanity check: WRT must be positive
        !                                  since ER is limited to WR/DT+RR in ebudget
        WRT_VH(I) = MAX(0., WR(I) - DT * (ER(I) -  RR(I)))
        !                                  Compute canopy drip
        !                                  if Wr > Wrmax, there is runoff
        RVEG_VH(I) = MAX(0., (WRT_VH(I) - WRMAX(I)) / DT ) 
        !
        !                                  Wr must be smaller than Wrmax
        !                                  after runoff is removed
        !
        WRT_VH(I) = MIN(WRT_VH(I), WRMAX(I))
     END IF

    ELSE
        WRT_VH(I) = 0.0
        RVEG_VH(I) = 0.
    ENDIF
     !

    IF(VEGL(I)*(1-PSN(I)).GE.EPSILON_SVS) THEN

     IF (SNM(I).GE.CRITSNOWMASS)THEN
        !                                  There is snow on the ground, remove evaporation
        !                                  then drain all liquid water from the vegetation
        WRT_VL(I) = WR(I) - DT * ER(I)
        !                                  Liquid water in vegetation becomes runoff
        RVEG_VL(I) = MAX(0.0,WRT_VL(I)/DT) ! Runoff mult. by dt later on, so water balance maintained..
        WRT_VL(I) = 0.0

     ELSE
        !                                  Remove P-E from liquid water in vegetation
        !                                  Sanity check: WRT must be positive
        !                                  since ER is limited to WR/DT+RR in ebudget
        WRT_VL(I) = MAX(0., WR(I) - DT * (ER(I) -  RR(I)))
        !                                  Compute canopy drip
        !                                  if Wr > Wrmax, there is runoff
        RVEG_VL(I) = MAX(0., (WRT_VL(I) - WRMAX(I)) / DT ) 
        !
        !                                  Wr must be smaller than Wrmax
        !                                  after runoff is removed
        !
        WRT_VL(I) = MIN(WRT_VL(I), WRMAX(I))
     END IF
    ELSE
        WRT_VL(I) = 0.0
        RVEG_VL(I) = 0.
    ENDIF

    IF ( (VEGH(I) + VEGL(I)*(1-PSN(I)) ) .GE. EPSILON_SVS ) THEN
       WRT(I) = (VEGL(I)*(1-PSN(I))*WRT_VL(I)+VEGH(I)*WRT_VH(I))/(VEGL(I)*(1-PSN(I))+VEGH(I))
    ELSE
       WRT(I) = 0.0
    ENDIF
     !
  END DO

  !
  !        2.     CALCULATE PRECIPITATIO and VEGETATION+SNOW RUNOFF REACHING THE GROUND 
  !               ------------------------------------------

  !
  !                                  Have vegetation runoff only if vegetation present
  !                                  Precipitation reaches ground directly only if no snow
  !                                  otherwise goes to snowpack, and reaches ground as snow runoff
  !
  DO I=1,N
      IF(VEGH(I)+VEGL(I)*(1.-PSN(I)).ge.EPSILON_SVS) THEN
         ! have vegetation -- have vegetation runoff 

         IF( SNM(I).GE.CRITSNOWMASS .AND. SVM(I).GE. CRITSNOWMASS ) THEN
            !both snow packs exists, rain falls directly to snow, consider runoff from vegetation also.
            PG(I) = (1.-VEGH(I)) * RSNOW(I) + VEGH(I) * RSNOWV(I) &
                 + VEGL(I)*(1.-PSN(I))*RVEG_VL(I) + VEGH(I) * RVEG_VH(I)


         ELSE IF ( SNM(I).GE.CRITSNOWMASS .AND. SVM(I).LT.CRITSNOWMASS ) THEN
            ! only low vegetation snow pack present, rain for low vegetation portion
            ! falls through to snow. No rain reaches bare ground because snowpack on low veg. and bare ground exists
            PG(I) = (1.-VEGH(I)) * RSNOW(I) &
                 + VEGL(I)*(1-PSN(I))*RVEG_VL(I) + VEGH(I) * RVEG_VH(I)


         ELSE IF ( SNM(I).LT.CRITSNOWMASS .AND. SVM(I).GE.CRITSNOWMASS ) THEN
            !only high vegetation snow pack present, rain for high vegetation portion
            ! falls through to snow. No snow on low veg. and bare ground, so rain can reach bare ground directly
            PG(I) = VEGH(I) * RSNOWV(I) &
               + VEGL(I)*(1.-PSN(I))*RVEG_VL(I) + VEGH(I) * RVEG_VH(I) &
               + (1. - PSN(I)) * (1. -VEGL(I) - VEGH(I) ) * RR(I)
         ELSE 
            ! no snow present, rain can reach bare ground directly
            PG(I) = VEGL(I)*(1.-PSN(I))*RVEG_VL(I) + VEGH(I) * RVEG_VH(I) &
                 + (1 - PSN(I)) * (1. -VEGL(I) - VEGH(I) ) * RR(I)           
         ENDIF

      ELSE
         ! bare ground only
         IF ( SNM(I).GE.CRITSNOWMASS) THEN
            ! have snow on bare ground, only consider snow runoff
            PG(I) = RSNOW(I)
         ELSE
            ! have 100% bare ground and now snow, all rain reaches surface
            PG(I) = RR(I)
         ENDIF
      ENDIF    

  ENDDO

  !
  !
  !        3.     CALCULATE THE REQUIRED PARAMETERS FOR NEW HYDRLOGY ROUTINE 
  !               ------------------------------------------

  DO I=1,N
     DO K=1,NL_SVS
        !Vectorize layer thicknesses over space for watdrn module
        DELZVEC(I,K)=DELZ(K)

        !Adjust ksat and wsat for presence of ice
        IF(KFICE==0) THEN
            FICE = (1.0-MAX(0.0,MIN((WSAT(I,K)-CRITWATER)/WSAT(I,K),WF(I,K)/WSAT(I,K))))**2.
        ELSE IF (KFICE ==1) THEN    
            FICE =  EXP(LOG(10.0)*(-6*WF(I,K)/(WF(I,K)+WD(I,K))))
        ELSE IF (KFICE ==3) THEN    
            FICE = 1.   
        ELSE
            FICE = 1.   
        ENDIF 
        KSATC(I,K) = KSAT(I,K)*FICE
        WSATC(I,K)= MAX((WSAT(I,K)-WF(I,K)-0.00001), CRITWATER)

        ! Calculate parameters needed for WATDRAIN
        ! ASAT0  bulk saturation at the begining of time step (WD/WDSAT)

        ASAT0 (I,K) = WD(I,K) / WSATC(I,K)
        ASATFC(I,K) = WFCINT(I,K)/ WSAT(I,K)

        ! Horizontal soil hydraulic conductivity decays with depth
        GRKSAT (I,K) = GRKSAT_C1 * EXP( GRKSAT_C2 *(DL_SVS(NL_SVS)-DL_SVS(K))/DL_SVS(NL_SVS))*KSATC(I,K)
        GRKEFL (I,K) = GRKEF(I)*GRKSAT(I,K)
     END DO
  END DO

  !        4.      CALCULATE SURFACE RUNOFF
  !              ------------------------------------------
  !Call watdrain to calculate runoff


  CALL WATDRN(DELZVEC(:,1),BCOEF(:,1),WSATC(:,1),GRKSAT(:,1),GRKEFL(:,1),ASATFC(:,1),ASAT0(:,1),ASAT1,SUBFLW,BASFLW,SATSFC,N,1,N,DT)


  DO I=1,N
     !use ksat to calculate runoff (mm/s)
     RUNOFF(I) = MAX( (SATSFC(I)*PG(I)+(1-SATSFC(I))*MAX(PG(I)-KSATC(I,1)*1000.,0.0)) , 0.0 )             

! EG_code related to ponding of water
     IF (lwater_ponding_svs1) THEN
        pond_ret(I) = min( RUNOFF(I)*DT , (maxpnd(I)-watpnd(I))*1000.0 )
        RUNOFF(I) = RUNOFF(I) - pond_ret(I)/DT
     ELSE
        pond_ret(I) = 0.0
     END IF

     ! assuming that 33% of urban cover is totally impervious (RUNOFF = PG over impervious surface)
     RUNOFF(I) = RUNOFF(I) * (1. - IMPERVU(I) ) + PG(I) * IMPERVU(I)        
     ! remove runoff from the ampount of water reaching the ground
     PG(I) = PG(I) - RUNOFF(I) - ( (1. - IMPERVU(I)) * pond_ret(I)/DT )           ! (mm/s)
     RUNOFF(I) = RUNOFF(I)*DT

     IF (lwater_ponding_svs1) THEN
        ! update amount of ponding water
        watpnd(I) = watpnd(I) + (1. - IMPERVU(I) ) * pond_ret(I) / 1000.0

     ENDIF 

  END DO


  !        5.      CALCULATE BASEFLOW AT THE BOTTOM OF THE TOTAL SOIL PROFILE
  !              -----------------------------------------

  ! BASE FLOW for last soil layer
  ! Compute wfc as a function of thickness of the total depth DL_SVS(NL_SVS) depth based on Soulis et al. (2012)

  DO I=1,N
     WFCDP(I) = WSATC(I,NL_SVS)*FBCOF(I,NL_SVS)*(PSISAT(I,NL_SVS)/DL_SVS(NL_SVS))**(1./BCOEF(I,NL_SVS))
  END DO

  !Call WATDRAIN to calculate baseflow

  CALL WATDRN(DELZVEC(:,NL_SVS),BCOEF(:,NL_SVS),WSATC(:,NL_SVS),KSATC(:,NL_SVS), &
       GRKEFL(:,NL_SVS),ASATFC(:,NL_SVS),ASAT0(:,NL_SVS),ASAT1,SUBFLW,BASFLW,SATSFC,N,1,N,DT)

  DO I=1,N

     IF(WD(I,NL_SVS).GT.WFCDP(I)) THEN  
        !baseflow only happens when soil water contant of last layer exceeds field capacity
        F(I,NL_SVS+1)=BASFLW(I)
     ELSE 
        F(I,NL_SVS+1)=0.0
     END IF

  END DO


  !        6.      CALCULATE NEW SOIL WATER CONTENT FOR EACH LAYER
  !              -----------------------------------------
  !
  ! Compute water removed by evapotranspiration from each layer (rate in m/s, grid box average)
  DO I=1,N
     DO K=1,NL_SVS
        ETR_GRID(I,K)=((VEGL(I)*(1.-PSN(I))+VEGH(I)*(1.-PSNVH(I)))*ACROOT(I,K)*ETR(I))/(1000.*DELZ(K))
     END DO
  END DO
  !
  !Compute water fluxes between soil layers based on Richard's equation 

  ! Boundary condition at the top is net precipitation      
  DO I=1,N
     F(I,1)=(PG(I)-(1.-VEGL(I)-VEGH(I))*(1-PSN(I))* EG(I))*DT/1000.
  END DO

  !Compute water fluxes between soil layers, find K AND PSI at the boundaries     
  ! do it for all layers except NL_SVS (water flux is computed above from watdrain).       

  IF (hydro_svs_method.EQ.0) THEN
     ! First-order forward method       
     CALL SOIL_FLUXES( DT, &
          WSATC, KSATC, PSISAT, BCOEF, ETR_GRID, WD, &
          F, WDT, DWD_RK1, OVER_RK1, KHC, PSI, N)
     DO I=1,N
        DO K=1,NL_SVS
           WDT(I,K)=WDT(I,K)-OVER_RK1(I,K)
        END DO
     END DO

  ELSE !hydro_svs_method=1
     ! Runge-Kutta 4th order method
     ! see for example http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
     ! Divide the source/sink terms by half to estimate WD at midpoint
     DO I=1,N
        DO K=1,NL_SVS+1

           F_RK(I,K) = F(I,K)/2.
        END DO
     END DO
     ! k1=f(y*(t0),t0)
     CALL SOIL_FLUXES( DT/2., &
          WSATC, KSATC, PSISAT, BCOEF, ETR_GRID, WD, &
          F_RK, WD_RK, DWD_RK1, OVER_RK1, KHC, PSI, N)      
     ! k2=(f(y*(t0)+k1*h/2,t0+h/2)
     CALL SOIL_FLUXES( DT/2., &
          WSATC, KSATC, PSISAT, BCOEF, ETR_GRID, WD_RK, &
          F_RK, WD_RK, DWD_RK2, OVER_RK2, KHC, PSI, N)
     ! k3=f(y*(t0)+k2*h/2,t0+h/2)
     CALL SOIL_FLUXES( DT, &
          WSATC, KSATC, PSISAT, BCOEF, ETR_GRID, WD_RK, &
          F, WD_RK, DWD_RK3, OVER_RK3, KHC, PSI, N)
     ! k4=f(y*(t0)+k3*h,t0+h)
     CALL SOIL_FLUXES( DT, &
          WSATC, KSATC, PSISAT, BCOEF, ETR_GRID, WD_RK, &
          F, WD_RK, DWD_RK4, OVER_RK4, KHC, PSI, N)
     ! y*(t0+h)=y*(t0)+(k1*h/6+k2*h/3+k3*h/3+k4*h/6)
     ! careful: DWD_RK1 and DWD_RK2 were calculated with a timestep h=DT/2
     ! whereas DWD_RK3 and DWD_RK4 was calculated with h=DT
     ! so k1*h/6=1/3*DWD_RK1, k2*h/3=2/3*DWD_RK2, k3*h/3=1/3*DWD_RK3 and k4*h/6=1/6*DWD_RK4
     DO I=1,N
        DO K=1,NL_SVS
           WDT(I,K)=WD(I,K)+( &
                DWD_RK1(I,K)/3.+ &
                DWD_RK2(I,K)*2./3.+ &
                DWD_RK3(I,K)/3.+ &
                DWD_RK4(I,K)/6.)
           WDT(I,K)=WDT(I,K)-( &
                OVER_RK1(I,K)/3.+ &
                OVER_RK2(I,K)*2./3.+ &
                OVER_RK3(I,K)/3.+ &
                OVER_RK4(I,K)/6.)
        END DO
     END DO
  END IF

  ! Check if the calculated WDT is between critical water (0.01) and wsat
  !if WDT is less than CRITWATER get water from the next layer
  !if WD is more than saturation, add eccess water to lateral flow or baseflow (F) 
  !update water flow (f) of each layer based on the new soil moisture

  DO I=1,N
     DO K=1,NL_SVS
        IF (WDT(I,K).LT.CRITWATER)  THEN

           ! if we are in the last soil layer, soil water content of the layer below cannot be updated
           IF(K.NE.NL_SVS) WDT(I,K+1)=WDT(I,K+1)- &
                (CRITWATER-WDT(I,K))*DELZ(K)/DELZ(K+1)
           F(I,K+1)=F(I,K+1)-(CRITWATER-WDT(I,K))*DELZ(K)
           WDT(I,K)=CRITWATER

        ELSEIF (WDT(I,K).GT.WSATC(I,K))  THEN

           IF (K.NE.KHYD) THEN

              IF(WAT_REDIS==1) THEN
                 IF(K.NE.NL_SVS) THEN 
                    KSATMEAN(I,K)=KSATC(I,K)*KSATC(I,K+1)*(DELZ(K)+&
                         DELZ(K+1))/(KSATC(I,K)*DELZ(K+1)+KSATC(I,K+1)*DELZ(K))
                    
                 ELSE
                    KSATMEAN(I,K) = KSATC(I,K)
                 ENDIF
                 
                 ! excess water removal via a combination of a downward and a lateral flux
                 ! fraction of excess water going into layer below
                 W = KSATMEAN(I,K)/(KSATMEAN(I,K)+GRKEFL(I,K)*DELZ(K))

                 ! Downward liquid water flux is limited to avoid the saturation of the layer below  
                 IF(K.NE.NL_SVS) THEN   
                    WAT_DOWN = MIN(W*(WDT(I,K)-WSATC(I,K))*DELZ(K),  &
                         MAX(0.,(WSATC(I,K+1)-WDT(I,K+1))*DELZ(K+1)) )
                 ELSE
                    WAT_DOWN = W*(WDT(I,K)-WSATC(I,K))*DELZ(K)
                 ENDIF
                 ! adding this fraction of excess water to the layer below
                 F(I,K+1)=F(I,K+1)+ WAT_DOWN
                 ! the rest of the excess water goes to lateral flow
                 LATFLW(I,K) = (WDT(I,K)-WSATC(I,K))*DELZ(K) - WAT_DOWN
                 ! if we are in the last soil layer, soil water content of the layer below cannot be updated
                 IF(K.NE.NL_SVS)  WDT(I,K+1)=WDT(I,K+1)+WAT_DOWN/DELZ(K+1)
                  
              ELSE

                 ! excess water removal via a combination of a downward and a lateral flux
                 ! fraction of excess water going into layer below
                 W = KSATC(I,K)/(KSATC(I,K)+GRKEFL(I,K)*DELZ(K))
                 ! adding this fraction of excess water to the layer below
                 F(I,K+1)=F(I,K+1)+W*(WDT(I,K)-WSATC(I,K))*DELZ(K)
                 ! the rest of the excess water goes to lateral flow
                 LATFLW(I,K)=(1.0-W)*(WDT(I,K)-WSATC(I,K))*DELZ(K)
                 ! if we are in the last soil layer, soil water content of the layer below cannot be updated
                 IF(K.NE.NL_SVS) WDT(I,K+1)=WDT(I,K+1)+W*(WDT(I,K)-WSATC(I,K))*DELZ(K)/DELZ(K+1)

              ENDIF
              
             
           ELSEIF (K.EQ.KHYD) THEN
              ! excess water removal via lateral flow
              LATFLW(I,K)=(WDT(I,K)-WSATC(I,K))*DELZ(K)
           END IF
           WDT(I,K)=WSATC(I,K)         

        END IF
     END DO  
  END DO

  !        7.      CALCULATE LATERAL FLOW
  !              -----------------------------------------
  !calculate parameter needed for WATDRAIN based on new WDT 

  DO I=1,N
     DO K=1,NL_SVS
        ASAT0 (I,K) = WDT(I,K) / WSATC(I,K)     
     END DO
  END DO

  !Call WATDRAIN to calculate SUBFLW from each layer 

  DO K=1,NL_SVS   

     CALL WATDRN(DELZVEC(:,K),BCOEF(:,K),WSATC(:,K),GRKSAT(:,K),GRKEFL(:,K), &
          ASATFC(:,K),ASAT0(:,K),ASAT1,SUBFLW,BASFLW,SATSFC,N,1,N,DT)

     DO I=1,N

        IF(WDT(I,K).GT.WFCINT(I,K))THEN
           !check SUBFLW does not exceed the availabe water in the layer
           SUBFLW(I) = MAX(0.,MIN(WDT(I,K)*DELZ(K),SUBFLW(I)))
           LATFLW (I,K) = LATFLW(I,K)+SUBFLW(I)
           !remove LATFLW from the layer
           WDT(I,K)  = WDT(I,K)-SUBFLW(I)/DELZ(K)
        ENDIF

     END DO

  END DO

  !---------------------------------------------------

  !        8.     TIME-EVOLUTION OF THE FROZEN SOIL WATER
  !               ---------------------------------------
  !

  DO I=1,N
     !        TODO
     do k=1,nl_svs
        WFT(I,K) = WF(I,K)   
     enddo
  END DO
  !
  !---------------------------------------------------
  !
  !--------------------------------------------------
  ! SET MINIMUM VALUES OF SOIL WATER AND MAKE SURE FROZEN SOIL WATER NON-NEGATIVE
  DO I=1,N
     DO K=1,NL_SVS
        WDT(I,K) = MAX(WDT(I,K),CRITWATER)
        !
        IF (WDT(I,K)+WFT(I,K).LT.CRITWATER) THEN
           WDT(I,K) = CRITWATER
           WFT(I,K) = 0.0
        END IF
        WFT(I,K) = MAX( WFT(I,K), 0.0   ) 
     ENDDO
     !

  END DO
  !--------------------------------------------------------
  !--------------------------------------------------
  ! CALCULATE DIAGNOSTICS, CONVERT UNITS, AND ACCUMULATORS 

  DO I=1,N

     DO K=1,NL_SVS
        ! CONVERT TO MM
        LATFLW (I,K) = LATFLW(I,K) * M_TO_MM
     ENDDO

     DO K=1,NL_SVS+1
        ! CONVERT TO MM
        F(I,K) = F(I,K) * M_TO_MM           
     ENDDO

  END DO

  RETURN
END SUBROUTINE HYDRO_SVS
