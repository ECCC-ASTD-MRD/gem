!copyright (C) 2001  MSC-RPN COMM  %%%RPNPHY%%%

      SUBROUTINE SNOW_ALONE (TSNS,TSND,RHOSL,ALPHAS,WL, & 
                           SNODP,SM, &   
                           PS,VMOD, VDIR, RHOA,THETAA,RG,RAT, &  
                           HU,RR,SR,T,T2M, &  
                           U10M,V10M,TAVG, &  
                           MELTS_TOT,MELTS_RN, &   
                           RNET,HFLUX,LE,EFLUX,RSNOW, &  
                           RHOSNO, RESA, &  
                           DT,Z0,Z0LOC,FCOR, ZUSL,ZTSL, LAT, PSN,N) 
!
      use sfclayer, only: sl_sfclayer,SL_OK
      use tdpack
      use sfc_options
      use svs_configs
      implicit none
!!!#include <arch_specific.hf>

!
      INTEGER N
!
      REAL TSNS(N), TSND(N), RHOSL(N), ALPHAS(N), WL(N)
      REAL SNODP(N), SM(N)

      REAL PS(N), VMOD(N), VDIR(N), RHOA(N), THETAA(N), RG(N), RAT(N)
      REAL HU(N), RR(N), SR(N), T(N), T2M(N)
      REAL U10M(N), V10M(N)
      REAL TAVG(N),MELTS_TOT(N), MELTS_RN(N)
      REAL RNET(N), HFLUX(N), LE(N), EFLUX(N), RSNOW(N)
      REAL RHOSNO(N), LAT(N), RESA(N)
      REAL DT,Z0(N), Z0LOC(N), FCOR(N), ZUSL(N), ZTSL(N),PSN(N)

!
!
!Author
!             S. Belair & M. Abrahamowicz (May 2009)
!
!Revisions
!
! 001         Bug fixes: M. Abrahamowicz, S.Z. Husain, S. Zhang, N. Gauthier,
!                        E. Gaborit, V. Vionnet          
!Object
!             Stand-alone snow model 
!
!Arguments
!
!
!
!             - Input/Output (Prognostic variables of the snow scheme) -
!
! TSNS        snow temperature -- S for "skin"
! TSND        mean snow temperature -- D for "deep"
! RHOSL       density of snow (relative to that of liquid water)
! ALPHAS      albedo of snow
! WL          liquid water in the snow pack
! SNODP       snow depth
! SM          snow mass (equivalent water content of the snow reservoir)
!
!             - Input (Forcing) -
! PS          Surface pressure
! VMOD        wind speed at the model lowest level (ZUSL)
! VDIR        wind direction at the model lowest level (ZUSL)
! RHOA        density of air at the model lowest level
! THETAA      Potential temperature at the model lowest level
! RG          solar radiation incident at the surface (downward solar)
! RAT         longwave radiation incident at the surface (NIR)
! HU          Specific humidity of air at the model lowest level
! RR          liquid precipitation rate at the surface in [mm/s]
! SR          solid  precipitation rate at the surface in [mm/s]
! T           surface air temperature
! T2M         2-m air temperature
! U10M        U-component of wind at 10 m
! V10M        V-component of wind at 10 m
!
!             - Other inputs -
! DT          time step
! Z0          momentum roughness length (no snow)
! Z0LOC      local momentum roughness (no snow)      
! FCOR        Coriolis factor
! LAT       latitude
! ZUSL        height of wind input(measured from model base at topo height + Z0)
! ZTSL        height of temperature and humidity input
!
!             - Output (Energy and water budgets of the snow pack) -
! TAVG        Average snow pack temperature use in melt/freez calc.
! MELTS_TOT   total snow melting (accumulator)
! MELTS_RN    snow melting due to rain (accumulator)
! RNET        net radiation at the snow surface
! HFLUX       sensible heat flux from the snow surface
! EFLUX       water vapor flux from the snow surface 
! LE          latent heat flux from the snow surface 
! RSNOW       liquid water out of the snow pack
! RHOSNO      density of snow (kg/m3) for output only
! RESA        aerodynamical surface resistance for snow
!
include "isbapar.cdk"

!
      INTEGER I
!
!
      REAL LAMI, CICE, DAY, CWAT, DT_12MIN
      
      REAL EMISSN, RAIN1, RAIN2, MLTRAIN 
      REAL CRMIN, CRMAX, TAUHOUR, RHOE,  MYOMEGA
      REAL RSNOW_DAY, RHOICE, ANSMIN, MAX_EFLUX

      real, dimension(n) :: lams, zcs, zqs, ctu, zqsat, zdqsat, zqsatt, &
           rora, a, b, c, tsnst, tsndt, rhomax, fmltrain, &
           smt, wlt, alphast, rhosfall, rhoslt, smx, kdiffu, &
           dampd, z0h, bcoef, dmelt, dsnowdt, ftemp, wlmax,  &
           freez_l1, work_l1, work_l2, melt_l1, melt_l2, melt_rain, rsnow_critmass
      
!
!
!
!***********************************************************************
!
!
!
!                                THE FOLLOWING SHOULD BE PUT IN 
!                                A COMMON COMDECK
!
!
      LAMI    = 2.22 
      CICE    = 2.106E3  ! specific heat of ice 
      DAY     = 86400.
      EMISSN  = 0.97
      CRMIN   = 0.03
      CRMAX   = 0.10
      RHOE    = 0.20
      TAUHOUR = 3600.
      RHOICE  = 0.9
      ANSMIN  = 0.5
      MYOMEGA   = ( 2*PI )/ DAY
      CWAT    = 4.187E3  ! Specific heat of water
      DT_12MIN = 12.*60. ! 12 min in sec     
      
!
      RAIN1   = 2.8e-5 ! mm/s
      RAIN2   = 2.8e-4 ! mm/s
!
!
!
!*            REFRESH ALL INPUT VARIABLES IF THERE IS NEGLIGIBLE SNOW
!                -----------------------------------------------------
!         CHECK THAT THIS INITIALIZATION MAKES SENSE
!
      DO I=1,N
         IF (SM(I).LT.CRITSNOWMASS) THEN
            ALPHAS(I)   = ANSMAX
            RHOSL(I)    = RHOSDEF
            !                             For snow temperature, set it to AIR temperature
!                             capped at the triple point of water
            TSNS(I)     = MIN(T(I),TRPL)
            TSND(I)     = MIN(T(I),TRPL)

            !                             To conserve water, send trace snow mass and liquid water in snow pack to runoff before setting mass to zero

            RSNOW_CRITMASS(I) = ( SM(I) + WL(I) ) / DT
            !                             Reset snow mass, liquid water and snow depth to zero
            SM(I)       = 0.  
            WL(I)       = 0.
            SNODP(I)    = 0.          
         ELSE
            RSNOW_CRITMASS(I) = 0.0
         END IF
      END DO
!
!
!
!
!
!*       1.     THE HEAT CAPACITY AND THERMAL DIFFUSIVITY OF THE SNOW
!               -----------------------------
!
!                          
!
      DO I=1,N
!                          First calculate the conductivity of snow (Yen,1981)
        LAMS(I) = LAMI * RHOSL(I)**1.88
!
!                          Heat capacity
!
        ZCS(I) = 2.0 * SQRT( PI/(LAMS(I)*RAUW*RHOSL(I)*CICE*DAY))
!                          
!                          Thermal diffusivity (You et al., 2014)
        KDIFFU(I) =  LAMS(I) / ( CICE * RAUW*RHOSL(I)) 
!
      END DO
!
!        2.     SPECIFIC HUMIDITY AT SNOW SURFACE
!               ----------------------------------------
!
!
!                       Calculate specific humidity at snow surface
!                       (For snow, specific & saturation humdity are
!                       the same as snow is always saturated)
!
      DO I=1,N
         ZQS(I) = FOQST( TSNS(I), PS(I) )
      END DO


!
!
!***     3.     SURFACE TRANSFER COEFFICIENTS FOR HEAT AND MOMENTUM (CH and CD)
!*             ---------------------------------------------------------------
!
!  
      DO I=1,N
        Z0H(I)      = Z0HSNOW
      END DO

      if ( svs_dynamic_z0h ) then
         ! here do not use dynamic z0h calculation
         ! however, to be consistent with bare ground and veg.
         ! use local z0m

         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
                          TSNS, ZQS, Z0LOC, Z0H, LAT, FCOR, &
                          L_min=sl_Lmin_soil, &
                          coeft=CTU )
         
         if (i /= SL_OK) then
            call physeterror('snow_alone', 'error returned by sl_sfclayer()')
            return
         endif

      else

         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
                          TSNS, ZQS, Z0, Z0H, LAT, FCOR, &
                          L_min=sl_Lmin_soil, &
                          coeft=CTU )

         if (i /= SL_OK) then
            call physeterror('snow_alone', 'error returned by sl_sfclayer()')
            return
         endif
      endif

      DO I=1,N
!
        RESA(I) = 1. / CTU(I)
!
      END DO


!
!
!
!       4.     TIME INTEGRATION of the SNOW SKIN TEMPERATURE (TSNS)
!               -------------------------------------------------
!
!                            Thermodynamic functions
!
      DO I=1,N
        ZQSAT(I)  = FOQST( TSNS(I),PS(I) )
        ZDQSAT(I) = FODQS( ZQSAT(I),TSNS(I) )
      END DO
!
!
!                            function zrsra
!
      DO I=1,N
        RORA(I) = RHOA(I) / RESA(I)
      END DO
!
!     
!
!                             terms za, zb, and zc for the
!                                    calculation of tsns(t)
!
      DO I=1,N
!
        A(I) = 1. / DT + ZCS(I) * &  
                (4. * EMISSN * STEFAN * (TSNS(I)**3) &   
                + RORA(I) * ZDQSAT(I) * (CHLC+CHLF) &  
                + RORA(I) * CPD) &  
                + 2. * PI / DAY
!
        B(I) = 1. / DT + ZCS(I) * &   
                (3. * EMISSN * STEFAN * (TSNS(I)** 3) &   
                + RORA(I) * ZDQSAT(I) * (CHLC+CHLF) )
!
        C(I) = 2. * PI * TSND(I) / DAY &   
                + ZCS(I) *  & 
                ( RORA(I) * CPD * THETAA(I)  &  
                + RG(I) * (1. - ALPHAS(I)) + EMISSN * RAT(I) &   
                - RORA(I) &  
                * (CHLC+CHLF) * (ZQSAT(I)-HU(I)) )
!
      END DO
!
!
      DO I=1,N
         TSNST(I) = ( TSNS(I)*B(I) + C(I) ) / A(I)
      ENDDO

!
!!       5.     DEEP SNOW TEMPERATURE (TSND) AT TIME 'T+DT'
!               ------------------------------------------
!
      DO I=1,N
        TSNDT(I) = (TSND(I) + DT*TSNST(I)/DAY) / (1.+DT/DAY)
      END DO


!
!
!
!
!*       6.     MELTING AND FREEZING TENDENCIES OF SNOW
!               ---------------------------------------
!
!                             
!                             Calculate tendencies using T+ temperatures
!                             Apply energy related to melt/freez to T+ afterwards
!
!
      DO I=1,N
!                             Calculate average snow pack temperature
!                             assuming no time-depency or shift in phase
!                             of the sinusoidal thermal wave (of the force-restore) 
!                             penetrating the snow pack (roughly following You et al., 2014)
!
!                             
!                             Damping depth in m , assuming diurnal forcing dominates
!                           
         DAMPD(I) = SQRT(  2 * KDIFFU(I) / MYOMEGA  )
!                             
!                             Average snow pack temp. for "superficial/damping" layer 
!        
         IF(SNODP(I).GT.0.0) THEN       
            DMELT(I) = MIN(DAMPD(I),SNODP(I))
            BCOEF(I) = ( 1. - exp( - DMELT(I)/DAMPD(I) ) ) * DAMPD(I) / DMELT(I)
         ELSE
            DMELT(I) = 0.0
            BCOEF(I) = 0.0
         ENDIF

         TAVG(I) = BCOEF(I) * TSNST(I) + (1-BCOEF(I)) * TSNDT(I)
         
!      
!                         Common portion of the MELT and FREEZ TERMS
!  
         ! "layer 1"
         WORK_L1(I) = (TAVG(I)-TRPL) / ( ZCS(I)*CHLF*DT )
         ! "layer 2"
         WORK_L2(I) = (TSNDT(I)-TRPL) / ( ZCS(I)*CHLF*DT )
      END DO
!
!
!
!                             MELT tendencies  -- NO FREEZING FOR NOW
!                             Also calculate the maximum snow density
!


      DO I=1,N
         
         ! layer 1 
         if( work_l1(i).ge.0.0 .and. SM(I).ge.CRITSNOWMASS ) then
            ! have melting
            MELT_L1(I)  = MIN( WORK_L1(I) , SM(I)*(DMELT(I)/SNODP(I))/DT )
            !RHOMAX(I)   = 0.6
            FREEZ_L1(I) = 0.0
            RHOMAX(I) = 600. - 204.7 / (SNODP(I)+EPSILON_SVS) *  & 
            ( 1.-EXP(-SNODP(I)/0.673))
            RHOMAX(I) = 0.001 * RHOMAX(I)
         else if( work_l1(i).lt.0.0 .and. SM(I).ge.CRITSNOWMASS ) then
            ! have freezing --- don't allow it for now
            MELT_L1(I)   = 0.0
            FREEZ_L1(I)  = 0.0
            !RHOMAX(I) = 0.3
                        
            !FREEZ_L1(I)  = MIN( -WORK_L1(I) , WL(I)/DT )
            RHOMAX(I) = 450. - 204.7 / (SNODP(I)+EPSILON_SVS) *  & 
              ( 1.-EXP(-SNODP(I)/0.673))
            RHOMAX(I) = 0.001 * RHOMAX(I)
         else
            freez_l1(i)=0.0
            melt_l1(I)=0.0
            rhomax(I)=RHOSDEF
         endif
         !
         ! layer 2
         !
         if(work_l2(i).gt.0.0  .and.  SM(I).ge.CRITSNOWMASS .and.dmelt(i).lt.snodp(i)) then
            ! melting
            MELT_L2(I)  = MIN( WORK_L2(I) , SM(I)*(((SNODP(I)-DMELT(I))/SNODP(I))/DT ))
         else
            MELT_L2(I)  = 0.0
         endif
        

      END DO
!
!
!
!        7.     EFFECT OF RAIN ON TEMPERATURE OF SNOW PACK
!               ------------------------------------------
!
!                                When rain is falling on snow,
!                                melting is accelerated due to heat
!                                transfers between the incident rain
!                                and the snow pack (since rain is
!                                usually warmer then the snow pack).
!
!                                It is hypothesized that the temperature
!                                of falling water is the same as that
!                                of air at diag. level.
!
    IF(SVS_SNOW_RAIN=='BELAIR03' .OR. SVS_SNOW_RAIN=='BELAIR03_DTGEM') THEN
      DO I=1,N
        IF (RR(I).LT.RAIN1) THEN
          FMLTRAIN(I) = 0.
        ELSE IF (RR(I).GT.RAIN2) THEN
          FMLTRAIN(I) = 1.
        ELSE
          FMLTRAIN(I) = ( RR(I) - RAIN1 ) / ( RAIN2 - RAIN1 )
        END IF
      END DO
!
      DO I=1,N
       
         IF (T2M(I).GT.TRPL.AND.SM(I).GT.0.0.AND.RR(I).GT.0.) THEN

            IF(SVS_SNOW_RAIN=='BELAIR03') THEN 
                ! Original formulation from Belair et al. (2003)
            MLTRAIN = ( T2M(I)-TRPL ) / ( 2.*ZCS(I)*CHLF*DT )
            ELSE IF(SVS_SNOW_RAIN=='BELAIR03_DTGEM') THEN 
                ! Revised formulation to remove the dependency of melt to the time step
                ! Use the time step of GEM at the time when the parameterisation has been developped (12 min)
                MLTRAIN = ( T2M(I)-TRPL ) / ( 2.*ZCS(I)*CHLF*DT_12MIN )
            ENDIF
           
            MELT_RAIN(I) = FMLTRAIN(I) * MLTRAIN
            
            MELT_RAIN(I) = MIN( MELT_RAIN(I), MAX( (SM(I)*(DMELT(I)/SNODP(I))/DT) - MELT_L1(I),0.0) )
            
         ELSE
            
            MLTRAIN = 0.0
            MELT_RAIN(I) = 0.0
        
        END IF
      END DO  

   ELSE IF (SVS_SNOW_RAIN =='NONE') THEN
      DO I=1,N
            MELT_RAIN(I) = 0.0
      END DO
 
   ELSE IF (SVS_SNOW_RAIN =='L21' ) THEN
 
     !  Assume all the energy brought by rain is used ot melt the snowpack.
     ! Parameterization used in Vionnet et al. (2020) and Leonardini et al. (2021)
      DO I=1,N
         IF (T2M(I).GT.TRPL.AND.SM(I).GT.0.0.AND.RR(I).GT.0.) THEN
            MELT_RAIN(I)  = RR(I)*  CWAT* ( T2M(I)-TRPL ) /CHLF
            MELT_RAIN(I) = MIN( MELT_RAIN(I), MAX( (SM(I)*(DMELT(I)/SNODP(I))/DT) - MELT_L1(I),0.0) )
         ELSE
            MELT_RAIN(I) = 0.0
        END IF
      END DO
 
    ENDIF
!
!
!
!
!                              Melting-Freezing tendency for the
!                              SM and WL reservoirs
!
      DO I=1,N
        DSNOWDT(I) = ( FREEZ_L1(I)  -MELT_L1(I)-MELT_L2(I)-MELT_RAIN(I) ) * DT
!
        MELTS_TOT(I) = MELTS_TOT(I) + (MELT_L1(I)+MELT_L2(I)+MELT_RAIN(I))*DT
        MELTS_RN(I) = MELTS_RN(I) + MELT_RAIN(I)*DT
      END DO

!
!
!        8.     EFFECT OF MELT/FREEZE ON SNOWPACK TEMP.
!               ------------------------------------------
!
!
!                              new temperatures Tsns(t) & Tsnd(t) after melting/freezing
!                            
!
      DO I=1,N
         TSNST(I) = TSNST(I) +   ZCS(I) * CHLF * (FREEZ_L1(I)-MELT_L1(I) )  * DT
         TSNDT(I) = TSNDT(I) +   ZCS(I) * CHLF * (           -MELT_L2(I) ) * DT
         
!       ALLOW TEMP. TO GO ABOVE ZERO TO TRY TO CONSERVE ENERGY WITH FORCE-RESTORE !!!
!       Make sure don't exceed triple pt.                             
!     
       ! TSNST(I) = MIN( TSNST(I) , TRPL )
        !TSNDT(I) = MIN( TSNDT(I) , TRPL )
!
      END DO      

!
!
!
!
!        9.     FLUX CALCULATIONS FOR SNOW COVERED SURFACE ONLY
!               ------------------------------------------------
!
!
      DO I=1,N
!                                            recalculate the qsat function
!
        ZQSATT(I) =  FOQST(  TSNST(I)  ,  PS(I)   )

!
!                                            net radiation
!
        RNET(I)  = (1. - ALPHAS(I)) * RG(I) + EMISSN *  & 
                 (RAT(I) - STEFAN * (TSNST(I)** 4))

!
!                                            sensible heat flux
!
        HFLUX(I) = RHOA(I) * CPD * (TSNST(I) - THETAA(I)) / RESA(I)


        FTEMP(I) = ( TSNST(I) - THETAA(I) ) / RESA(I)
!
!                                            latent heat of evaporation from
!                                            the snow canopy
!                                            
!
        EFLUX(I) = (ZQSATT(I) - HU(I)) / RESA(I)

!       IMPOSE MAXIMUM on EFLUX, based on available liquid water after melt/freez of reservoir... make sure latent heat consistent with this...
        if(PSN(I).gt.0.0) then
           MAX_EFLUX= ( (SM(I) + DSNOWDT(I)) / DT + SR(I) )  / RHOA(I) / PSN(I) 
        else
           MAX_EFLUX=0.0
        endif

        EFLUX(I)=  MIN( EFLUX(I), MAX_EFLUX )
              
!                                   
        LE(I)    = RHOA(I) * (CHLC+CHLF) * EFLUX(I)

      END DO


!
!
!
!!       10.     EVOLUTION OF THE SNOW EQUIVALENT WATER CONTENT (Wst)
!               --------------------------------------------
!
!
      DO I=1,N
!
!
!                               evolution of Ws
!
!            ! DANGER HERE: We use SM as a check of snow presence,
!            ! In the case of no-snow this variable could be updated
!            ! anyhow due to the construct of the subroutine, so
!            ! we add an additional check here i.e. 
! 
!!           CALCULATE SMT VARIABLE IF AND ONLY IF:
!            A) there is snow i.e. ws >= critsnow
!         or B) the snow rate is non-zero
!              
        IF(SM(I).ge.CRITSNOWMASS.or.SR(I).gt.0.0) THEN
           SMT(I) = SM(I) - DT * (RHOA(I)*PSN(I)*EFLUX(I) - SR(I)) + DSNOWDT(I) 
           SMT(I) = MAX( SMT(I), 0.0)
        ELSE
           SMT(I) = 0.0
        ENDIF
        
      END DO

!
!
!       11.     EVOLUTION OF LIQUID WATER IN THE SNOW PACK
!               ------------------------------------------
!
!                               Calculate the maximum liquid water
!                               that can be retained in the snow pack
!
      DO I=1,N
        IF (RHOSL(I).LT.RHOE) THEN
          WLMAX(I) = ( CRMIN + (CRMAX-CRMIN)*(RHOE-RHOSL(I))/ RHOE)* SMT(I)
        ELSE
          WLMAX(I) = CRMIN * SMT(I)
        END IF
      END DO


!
!
!                               Calculate runoff of liquid water from 
!                               the snow pack
!
      DO I=1,N
        IF (WL(I).LE.WLMAX(I)) THEN
          RSNOW(I) = ( WL(I) / TAUHOUR ) * EXP( WL(I)-WLMAX(I) )
        ELSE
          RSNOW(I) = WLMAX(I) / TAUHOUR + (WL(I)-WLMAX(I)) / DT
        END IF
        RSNOW(I) = MAX( 0.      , RSNOW(I) )
        RSNOW(I) = MIN( (WL(I)-DSNOWDT(I))/DT+RR(I) ,  RSNOW(I) )
      END DO
!
!
!                               Calculate the new values for WL and
!                               for the liquid water reaching the ground
!
      DO I=1,N
         IF(SM(I).ge.CRITSNOWMASS.or.SR(I).gt.0.0) THEN
            !if existing snow, or fresh snow fall
            WLT(I) = WL(I) +  RR(I) * DT - RSNOW(I)* DT - DSNOWDT(I)
            WLT(I) = MAX( 0., WLT(I) )
         ELSE
            WLT(I) = 0.0
         ENDIF
      END DO
!
!
!
!!      12.     EVOLUTION OF SNOW ALBEDO
!               ------------------------
!
!
!
      DO I=1,N
        RSNOW_DAY = MIN( RSNOW(I)*DAY, 10. )
      END DO
!
!                                       the evolution of the snow albedo differs
!                                       if there is melting or not
!
      DO I=1,N
!
!
        IF (SMT(I).GT.0.0.AND.DSNOWDT(I).LT.0.0) THEN
!
!                                       when there is freezing
!
           ALPHAST(I) = (ALPHAS(I)-ANSMIN)*EXP(-0.01*DT/3600.) &  
                +  ANSMIN & 
                +  SR(I)*DT/WCRN*(ANSMAX-ANSMIN)
!
!
        ELSE IF (SMT(I).GT.0.0.AND.DSNOWDT(I).GE.0.0) THEN
!
!                                       when there is melting
!
           ALPHAST(I) = ALPHAS(I) - TODRY*DT/DAY  &   
                + SR(I)*DT/WCRN*(ANSMAX-ANSMIN)
!
!
        ELSE
           ALPHAST(I) = ANSMAX
        ENDIF
!                                       limits of the albedo
!
        ALPHAST(I) = MAX( ANSMIN, ALPHAST(I) )       
        ALPHAST(I) = MIN( ANSMAX, ALPHAST(I) )
!
      END DO
!
!
!
!*       13.     EVOLUTION OF SNOW DENSITY
!                -------------------------
!
!
!                           Density of falling snow
!
      DO I=1,N
        IF (SMT(I).GT.0.0) THEN
           RHOSFALL(I) = 109. + 6.*(T2M(I)-TRPL) +   &
                              26.*(U10M(I)**2+V10M(I)**2)**0.25
           RHOSFALL(I) = MIN(MAX((RHOSFALL(I)*0.001),RHOMIN), 0.250)
        ELSE
           RHOSFALL(I) = RHOSDEF
        END IF
      END DO


!
!                           Evolution of the snow density depends
!                           on 3 factors:  
!    
!                           - decrease of density due to fresh new snow
!                           - increase of density due to aging
!                           - increase of density due to freezing of 
!                             liquid water in the snow pack
!
!
!                           A) decrease due to fresh new snow
!
      DO I=1,N
         SMX(I)   = MAX( SMT(I),SR(I)*DT)
         IF (SMT(I).GT.0.0) THEN
!            
            RHOSLT(I) = ( (SMX(I)-SR(I)*DT) * RHOSL(I)  & 
                        + (SR(I)*DT) * RHOSFALL(I)) / SMX(I)
         ELSE
            ! default
            RHOSLT(I) = RHOSDEF
         END IF
      END DO
!
!
!                           B) increase due to aging
!
      DO I=1,N
        IF (SMT(I).GT.0.0.AND.RHOSLT(I).LT.RHOMAX(I)) THEN
           RHOSLT(I) = (RHOSLT(I)-RHOMAX(I))*EXP(-0.01*DT/3600.)  &  
                          + RHOMAX(I)
        END IF
      END DO
!
!
!                           C) increase due to freezing
!
      DO I=1,N
        IF (SMT(I).GT.0.0) THEN
          RHOSLT(I) =  ( SMT(I)*RHOSLT(I) + FREEZ_L1(I)*DT*RHOICE ) &  
                        / ( SMT(I) + FREEZ_L1(I) * DT )
!                          Make sure within bounds 
          RHOSLT(I) = MIN( RHOICE, RHOSLT(I) )
          RHOSLT(I) = MAX( RHOMIN, RHOSLT(I) )

        END IF
      END DO
!
!                          Calculate snow depth based on snow density
!
      DO I=1,N
         SNODP(I) = SMT(I)/(RHOSLT(I)*RAUW)
!
!
      END DO
!
!
!
!*       14.     UPDATE the PROGNOSTIC VARIABLES
!                -------------------------------
!
      DO I=1,N
        TSNS(I)     = TSNST(I)
        TSND(I)     = TSNDT(I)
        WL(I)       = WLT(I)
        SM(I)       = SMT(I)
        ALPHAS(I)   = ALPHAST(I)
        RHOSL(I)    = RHOSLT(I)
        RHOSNO(I)   = RHOSLT(I)*RAUW
      END DO
!
!
!
!                  If negligible amount of snow,
!                  set ALL OUTPUT variables 
!                  to default values and/or zero.
!                 
!
      DO I=1,N
        IF (SM(I).LT.CRITSNOWMASS) THEN
          ALPHAS(I)   = ANSMAX
          RHOSL(I)    = RHOSDEF
          RHOSNO(I)   = RHOSLT(I)*RAUW
          TSNS(I)     = 300.0
          TSND(I)     = 300.0
          TAVG(I)     = 300.0
          RNET(I)     = 0.0
          HFLUX(I)    = 0.0
          LE(I)       = 0.0
          EFLUX(I)    = 0.0
          ! To conserve water, send trace snow mass and
          ! liquid water in snow pack to runoff before setting mass to zero
          RSNOW_CRITMASS(I) = RSNOW_CRITMASS(I) + ( SM(I) + WL(I) ) / DT
          !                             Reset snow mass, liquid water and snow depth to zero
          SM(I)       = 0.  
          WL(I)       = 0.
          SNODP(I)    = 0.            
        END IF
        ! Add trace runoff to runoff variable
        RSNOW(I) = RSNOW(I) + RSNOW_CRITMASS(I)
      END DO
!
      RETURN
      END
