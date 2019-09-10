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

subroutine snow_veg1(TSNS,TSND,RHOSL,ALPHAS,WL, &
     SNODP,SM, &
     PS,VMOD,VDIR,RHOA,THETAA,RG,RAT, &
     HU,RR,SR,T,T2M, &
     U10M,V10M,VEGH,SKYVIEW,VTR, &
     TVGD, EMISVH,TAVG, &
     RNET,HFLUX,LE,EFLUX,RSNOW, &
     RHOSNO, RESA, &
     DT,Z0,Z0HSNOW,FCOR,LAT,ZUSL,ZTSL,PSNVHA, N)
   use tdpack
   use sfclayer_mod, only: sl_sfclayer,SL_OK
   use sfc_options
   use svs_configs
   implicit none
!!!#include <arch_specific.hf>

   integer N

   real TSNS(N), TSND(N), RHOSL(N), ALPHAS(N), WL(N)
   real SNODP(N), SM(N), SKYVIEW(N), VTR(N), EMISVH(N), VEGH(N)
   real TVGD(N)

   real PS(N), VMOD(N), VDIR(N), RHOA(N), THETAA(N), RG(N), RAT(N)
   real HU(N), RR(N), SR(N), T(N), T2M(N)
   real U10M(N), V10M(N)
   real RNET(N), HFLUX(N), LE(N), EFLUX(N), RSNOW(N)
   real RHOSNO(N), TAVG(N), RESA(N)
   real DT, Z0HSNOW, Z0(N), FCOR(N), LAT(N), ZUSL(N), ZTSL(N),PSNVHA(N)

   !@Author S. Belair & M. Abrahamowicz (May 2009)

   !@Revisions
   ! 001         M. Abrahamowicz, S.Z. Husain, S. Zhang, N. Gauthier,
   !                        E. Gaborit, V. Vionnet

   !@Object Stand-alone snow model for snow under vegetation

   !@Arguments
   !       - Input/Output (Prognostic variables of the snow-under-veg scheme) -
   ! TSNS        snow-under-veg temperature -- S for "skin"
   ! TSND        mean snow-under-veg temperature -- D for "deep"
   ! RHOSL       density of snow-under-veg (relative to that of liquid water)
   ! ALPHAS      albedo of snow-under-veg
   ! WL          liquid water in the snow-under-veg pack
   ! SNODP       snow-under-veg depth
   ! SM          snow-under-veg mass (equivalent water content of the snow reservoir)
   ! SKYVIEW     Sky view factor for tall/high vegetation
   ! VTR        (HIGH) Vegetation transmissivity
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
   ! VEGH        fraction of HIGH vegetation [0-1]
   ! SKYVIEW     Sky view factor for tall/high vegetation
   ! VTR         (HIGH) Vegetation transmissivity
   ! TVGD        Mean/Deep Vegetation Temperature
   ! EMISVH      High Vegetation Emissivity
   !             - Other inputs -
   ! DT          time step
   ! Z0          momentum roughness length (no snow)
   ! Z0HSNOW     Constant for thermal roughness of snow
   ! FCOR        Coriolis factor
   ! LAT         latitude
   ! ZUSL          height of wind input(measured from model base at topo height + Z0)
   ! ZTSL          height of temperature and humidity input

   !             - Output (Energy and water budgets of the snow pack) -
   ! TAVG        Average snow pack temperature use in melt/freez calc.
   ! RNET        net radiation at the snow surface
   ! HFLUX       sensible heat flux from the snow surface
   ! EFLUX       water vapor flux from the snow-unde-veg surface
   ! LE          latent heat flux from the snow-unde-veg surface
   ! RSNOW       liquid water out of the snow pack
   ! RHOSNO      density of snow-under-veg (kg/m3) for output only
   ! RESA        aerodynamical surface resistance for snow

   include "isbapar.cdk"

   integer I


   real LAMI, CICE, DAY

   real EMISSN, RAIN1, RAIN2, MLTRAIN
   real CRMIN, CRMAX, TAUHOUR, RHOE, MYOMEGA
   real RSNOW_DAY, RHOICE, ANSMIN
   real MAX_EFLUX

   real, dimension(n) :: lams, zcs, zqs, ctu, zqsat, zdqsat, zqsatt, &
        rora, a, b, c, tsnst, tsndt, rhomax, fmltrain, &
        smt, wlt, alphast, rhosfall, rhoslt, smx, rvrun, kdiffu, &
        dampd, z0h, bcoef, dmelt, dsnowdt, ftemp, wlmax, rsnow_critmass, &
        freez_l1, work_l1, work_l2, melt_l1, melt_l2, melt_rain




   !************************************************************************



   !                                THE FOLLOWING SHOULD BE PUT IN
   !                                A COMMON COMDECK


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

   RAIN1   = 2.8e-5 ! mm/s
   RAIN2   = 2.8e-4 ! mm/s



   !*            REFRESH ALL INPUT VARIABLES IF THERE IS NEGLIGIBLE SNOW
   !                -----------------------------------------------------
   !         CHECK THAT THIS INITIALIZATION MAKES SENSE

   !         If no high vegetation present, reset all snow variables because they do not exist

   !         IF snow mass less than critsnowmass, but high vegetation present, reset snow variables, but
   !         send trace snow mass and liquid water content to runoff
   do I=1,N
      if (VEGH(I).ge.EPSILON_SVS) then

         if (SM(I).lt.CRITSNOWMASS) then
            ALPHAS(I)   = ANSMAX
            RHOSL(I)    = RHOSDEF
            !                             For snow temperature, set it to AIR temperature
            !                             capped at the triple point of water
            TSNS(I)     = min(T(I),TRPL)
            TSND(I)     = min(T(I),TRPL)

            !                             To conserve water, send trace snow mass and liquid water in snow pack to runoff before setting mass to zero
            RSNOW_CRITMASS(I) = ( SM(I) + WL(I) ) / DT
            !                             Reset snow mass, liquid water and snow depth to zero
            SM(I)       = 0.
            WL(I)       = 0.
            SNODP(I)    = 0.
         else
            RSNOW_CRITMASS(I) = 0.0
         end if


      else
         ! no high veg, so snow under high veg does not exist... reset all snow variables
         ALPHAS(I)   = ANSMAX
         RHOSL(I)    = RHOSDEF
         TSNS(I)     = min(T(I),TRPL)
         TSND(I)     = min(T(I),TRPL)
         SNODP(I)    = 0.
         WL(I)       = 0.
         SNODP(I)    = 0.
         SM(I)       = 0.
         RSNOW_CRITMASS(I) = 0.
      endif

   end do





   !*       1.     THE HEAT CAPACITY AND THERMAL DIFFUSIVITY OF THE SNOW
   !               -----------------------------



   do I=1,N
      !                          First calculate the conductivity of snow (Yen,1981)
      LAMS(I) = LAMI * RHOSL(I)**1.88

      !                          Heat capacity

      ZCS(I) = 2.0 * sqrt( PI/(LAMS(I)*RAUW*RHOSL(I)*CICE*DAY))

      !                          Thermal diffusivity (You et al., 2014)
      KDIFFU(I) =  LAMS(I) / ( CICE * RAUW*RHOSL(I))

   end do

   !        2.     SPECIFIC HUMIDITY AT SNOW SURFACE
   !               ----------------------------------------


   !                       Calculate specific humidity at snow surface
   !                       (For snow, specific & saturation humdity are
   !                       the same as snow is always saturated)

   do I=1,N
      ZQS(I) = FOQST( TSNS(I), PS(I) )
   end do




   !***     3.     SURFACE TRANSFER COEFFICIENTS FOR HEAT AND MOMENTUM (CH and CD)
   !*             ---------------------------------------------------------------


   do I=1,N
      Z0H(I)      = Z0HSNOW
   end do



   i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
        TSNS, ZQS, Z0, Z0H, LAT, FCOR, &
        L_min=sl_Lmin_soil, &
        coeft=CTU )

   if (i /= SL_OK) then
      call physeterror('snow_veg', 'error returned by sl_sfclayer()')
      return
   endif

   do I=1,N

      RESA(I) = 1. / CTU(I)

   end do


   !       4.     TIME INTEGRATION of the SNOW SKIN TEMPERATURE (TSNS)
   !               -------------------------------------------------

   !                            Thermodynamic functions

   do I=1,N
      ZQSAT(I)  = FOQST( TSNS(I),PS(I) )
      ZDQSAT(I) = FODQS( ZQSAT(I),TSNS(I) )
   end do


   !                            function zrsra

   do I=1,N
      RORA(I) = RHOA(I) / RESA(I)
   end do



   !                             terms za, zb, and zc for the
   !                                    calculation of tsns(t)
   !                             Note: Am using emissivity of one
   !                             for vegetation

   do I=1,N

      A(I) = 1. / DT + ZCS(I) * &
           (4. * EMISSN * STEFAN * (TSNS(I)**3) &
           + RORA(I) * ZDQSAT(I) * (CHLC+CHLF) &
           + RORA(I) * CPD) &
           + 2. * PI / DAY

      B(I) = 1. / DT + ZCS(I) * &
           (3. * EMISSN * STEFAN * (TSNS(I)** 3) &
           + RORA(I) * ZDQSAT(I) * (CHLC+CHLF) )

      C(I) = 2. * PI * TSND(I) / DAY &
           + ZCS(I) * &
           ( RORA(I) * CPD * THETAA(I)&
           + VTR(I) * RG(I) * (1. - ALPHAS(I))&
           + EMISSN *       SKYVIEW(I)  * RAT(I) &
           + EMISSN * (1. - SKYVIEW(I)) * EMISVH(I) * STEFAN * (TVGD(I)**4)&
           - RORA(I)  &
           * (CHLC+CHLF) * (ZQSAT(I)-HU(I)) )

   end do


   do I=1,N
      TSNST(I) = ( TSNS(I)*B(I) + C(I) ) / A(I)
   enddo


   !!       5.     DEEP SNOW TEMPERATURE (TSND) AT TIME 'T+DT'
   !               ------------------------------------------

   do I=1,N
      TSNDT(I) = (TSND(I) + DT*TSNST(I)/DAY) / (1.+DT/DAY)
   end do






   !*       6.     MELTING AND FREEZING TENDENCIES OF SNOW
   !               ---------------------------------------


   !                             Calculate tendencies using T+ temperatures
   !                             Apply energy related to melt/freez to T+ afterwards


   do I=1,N
      !                             Calculate average snow pack temperature
      !                             assuming no time-depency or shift in phase
      !                             of the sinusoidal thermal wave (of the force-restore)
      !                             penetrating the snow pack (roughly following You et al., 2014)


      !                             Damping depth in m , assuming diurnal forcing dominates

      DAMPD(I) = sqrt(  2 * KDIFFU(I) / MYOMEGA  )

      !                             Average snow pack temp. for "superficial/damping" layer

      if(SNODP(I).gt.0.0) then
         DMELT(I) = min(DAMPD(I),SNODP(I))
         BCOEF(I) = ( 1. - exp( - DMELT(I)/DAMPD(I) ) ) * DAMPD(I) / DMELT(I)
      else
         DMELT(I) = 0.0
         BCOEF(I) = 0.0
      endif

      TAVG(I) = BCOEF(I) * TSNST(I) + (1-BCOEF(I)) * TSNDT(I)


      !                         Common portion of the MELT and FREEZ TERMS

      ! "layer 1"
      WORK_L1(I) = (TAVG(I)-TRPL) / ( ZCS(I)*CHLF*DT )
      ! "layer 2"
      WORK_L2(I) = (TSNDT(I)-TRPL) / ( ZCS(I)*CHLF*DT )
   end do



   !                             MELT tendencies  -- NO FREEZING FOR NOW
   !                             Also calculate the maximum snow density



   do I=1,N

      ! layer 1
      if( work_l1(i).gt.0.0 .and. SM(I).ge.CRITSNOWMASS ) then
         ! have melting
         MELT_L1(I)  = min( WORK_L1(I) , SM(I)*(DMELT(I)/SNODP(I))/DT )
         RHOMAX(I)   = 0.6
         FREEZ_L1(I) = 0.0
         ! RHOMAX(I) = 600. - 20.47 / (SNODP(I)+EPSILON_SVS) *  &
         !( 1.-EXP(-SNODP(I)/0.0673))
         !RHOMAX(I) = 0.001 * RHOMAX(I)
      else if( work_l1(i).lt.0.0 .and. SM(I).ge.CRITSNOWMASS ) then
         ! have freezing --- don't allow it for now
         MELT_L1(I)   = 0.0
         FREEZ_L1(I)  = 0.0
         RHOMAX(I) = 0.3

         !FREEZ_L1(I)  = MIN( -WORK_L1(I) , WL(I)/DT )
         !RHOMAX(I) = 450. - 20.47 / (SNODP(I)+EPSILON_SVS) *  &
         !  ( 1.-EXP(-SNODP(I)/0.0673))
         !RHOMAX(I) = 0.001 * RHOMAX(I)
      else
         freez_l1(i)=0.0
         melt_l1(I)=0.0
         rhomax(I)=RHOSDEF
      endif

      ! layer 2

      if(work_l2(i).gt.0.0  .and.  SM(I).ge.CRITSNOWMASS .and.dmelt(i).lt.snodp(i)) then
         ! melting
         MELT_L2(I)  = min( WORK_L2(I) , SM(I)*(((SNODP(I)-DMELT(I))/SNODP(I))/DT ))
      else
         MELT_L2(I)  = 0.0
      endif


   end do



   !        7.     EFFECT OF RAIN ON TEMPERATURE OF SNOW PACK
   !               ------------------------------------------

   !                                When rain is falling on snow,
   !                                melting is accelerated due to heat
   !                                transfers between the incident rain
   !                                and the snow pack (since rain is
   !                                usually warmer then the snow pack).

   !                                It is hypothesized that the temperature
   !                                of falling water is the same as that
   !                                of air at diag. level.

   do I=1,N
      if (RR(I).lt.RAIN1) then
         FMLTRAIN(I) = 0.
      else if (RR(I).gt.RAIN2) then
         FMLTRAIN(I) = 1.
      else
         FMLTRAIN(I) = ( RR(I) - RAIN1 ) / ( RAIN2 - RAIN1 )
      end if
   end do

   do I=1,N

      if (T2M(I).gt.TRPL.and.SM(I).gt.0.0.and.RR(I).gt.0.) then
         MLTRAIN = ( T2M(I)-TRPL ) / ( 2.*ZCS(I)*CHLF*DT )

         MELT_RAIN(I) = FMLTRAIN(I) * MLTRAIN

         MELT_RAIN(I) = min( MELT_RAIN(I), max( (SM(I)*(DMELT(I)/SNODP(I))/DT) - MELT_L1(I),0.0) )

      else

         MLTRAIN = 0.0
         MELT_RAIN(I) = 0.0

      end if
   end do




   !                              Melting-Freezing tendency for the
   !                              SM and WL reservoirs

   do I=1,N
      DSNOWDT(I) = ( FREEZ_L1(I)  -MELT_L1(I)-MELT_L2(I)-MELT_RAIN(I) ) * DT

   end do



   !        8.     EFFECT OF MELT/FREEZE ON SNOWPACK TEMP.
   !               ------------------------------------------


   !                              new temperatures Tsns(t) & Tsnd(t) after melting/freezing


   do I=1,N
      TSNST(I) = TSNST(I) +   ZCS(I) * CHLF * (FREEZ_L1(I)-MELT_L1(I) )  * DT
      TSNDT(I) = TSNDT(I) +   ZCS(I) * CHLF * (           -MELT_L2(I) ) * DT

      !       ALLOW TEMP. TO GO ABOVE ZERO TO TRY TO CONSERVE ENERGY WITH FORCE-RESTORE !!
      !       Make sure don't exceed triple pt.

      ! TSNST(I) = MIN( TSNST(I) , TRPL )
      !TSNDT(I) = MIN( TSNDT(I) , TRPL )

   end do





   !        9.     FLUX CALCULATIONS FOR SNOW COVERED SURFACE ONLY
   !               ------------------------------------------------


   do I=1,N
      !                                            recalculate the qsat function

      ZQSATT(I) =  FOQST(  TSNST(I)  ,  PS(I)   )


      !                                            net radiation

      RNET(I)  = VTR(I) * (1. - ALPHAS(I)) * RG(I)&
           + EMISSN * (SKYVIEW(I)* RAT(I) - STEFAN * (TSNST(I)** 4))&
           + EMISSN * (1. - SKYVIEW(I))* EMISVH(I) * STEFAN * (TVGD(I)**4)


      !                                            sensible heat flux

      HFLUX(I) = RHOA(I) * CPD * (TSNST(I) - THETAA(I)) / RESA(I)


      FTEMP(I) = ( TSNST(I) - THETAA(I) ) / RESA(I)

      !                                            latent heat of evaporation from
      !                                            the snow canopy


      EFLUX(I) = (ZQSATT(I) - HU(I)) / RESA(I)

      !       IMPOSE MAXIMUM on EFLUX, based on available liquid water after melt/freez of reservoir... make sure latent heat consistent with this...

      if(PSNVHA(I).gt.0.0) then
         MAX_EFLUX=(  (SM(I) + DSNOWDT(I)) / DT + SR(I) )  / RHOA(I) / PSNVHA(I)
      else
         MAX_EFLUX=0.0
      endif

      EFLUX(I)=  min( EFLUX(I), MAX_EFLUX )


      LE(I)    = RHOA(I) * (CHLC+CHLF) * EFLUX(I)

   end do





   !!       10.     EVOLUTION OF THE SNOW EQUIVALENT WATER CONTENT (Wst)
   !               --------------------------------------------

   !                               evaporation over snow

   do I=1,N



      !                               evolution of Ws

      !            ! DANGER HERE: We use SM as a check of snow presence,
      !            ! In the case of no-snow this variable could be updated
      !            ! anyhow due to the construct of the subroutine, so
      !            ! we add an additional check here i.e.

      !!           CALCULATE SMT VARIABLE IF AND ONLY IF:
      !            A) there is snow i.e. ws >= critsnow
      !         or B) the snow rate is non-zero

      if(SM(I).ge.CRITSNOWMASS.or.SR(I).gt.0.0) then
         SMT(I) = SM(I) - DT * (RHOA(I)*PSNVHA(I)*EFLUX(I) - SR(I)) + DSNOWDT(I)
         SMT(I) = max(SMT(I), 0.0)
      else
         SMT(I) = 0.0
      endif

   end do



   !       11.     EVOLUTION OF LIQUID WATER IN THE SNOW PACK
   !               ------------------------------------------

   !                               Calculate the contribution of rain
   !                               and vegetation runnoff to WL
   !                               WATCH OUT: Assume all rain goes directly
   !                               to snowpack ... essential to close the water
   !                               budget

   do I=1,N
      RVRUN(I)= RR(I)
   end do




   !                               Calculate the maximum liquid water
   !                               that can be retained in the snow pack

   do I=1,N
      if (RHOSL(I).lt.RHOE) then
         WLMAX(I) = ( CRMIN + (CRMAX-CRMIN)*(RHOE-RHOSL(I))/ RHOE)* SMT(I)
      else
         WLMAX(I) = CRMIN * SMT(I)
      end if
   end do




   !                               Calculate runoff of liquid water from
   !                               the snow pack

   do I=1,N
      if (WL(I).le.WLMAX(I)) then
         RSNOW(I) = ( WL(I) / TAUHOUR ) * exp( WL(I)-WLMAX(I) )
      else
         RSNOW(I) = WLMAX(I) / TAUHOUR + (WL(I)-WLMAX(I)) / DT
      end if
      RSNOW(I) = max( 0.      , RSNOW(I) )
      RSNOW(I) = min( (WL(I)-DSNOWDT(I))/DT+RVRUN(I) ,  RSNOW(I) )
   end do


   !                               Calculate the new values for WL and
   !                               for the liquid water reaching the ground

   do I=1,N

      if(SM(I).ge.CRITSNOWMASS.or.SR(I).gt.0.0) then
         !if existing snow, or fresh snow fall
         WLT(I) = WL(I) + RVRUN(I) * DT - RSNOW(I)* DT - DSNOWDT(I)
         WLT(I) = max( 0.0, WLT(I) )
      else
         WLT(I) = 0.0
      endif
   end do



   !!      12.     EVOLUTION OF SNOW ALBEDO
   !               ------------------------



   do I=1,N
      RSNOW_DAY = min( RSNOW(I)*DAY, 10. )
   end do

   !                                       the evolution of the snow albedo differs
   !                                       if there is melting or not

   do I=1,N


      if (SMT(I).gt.0.0.and.DSNOWDT(I).lt.0.0) then

         !                                       when there is freezing

         ALPHAST(I) = (ALPHAS(I)-ANSMIN)*exp(-0.01*DT/3600.) &
              +  ANSMIN &
              +  SR(I)*DT/WCRN*(ANSMAX-ANSMIN)


      else if (SMT(I).gt.0.0.and.DSNOWDT(I).ge.0.0) then

         !                                       when there is melting

         ALPHAST(I) = ALPHAS(I) - TODRY*DT/DAY  &
              + SR(I)*DT/WCRN*(ANSMAX-ANSMIN)


      else
         ALPHAST(I) = ANSMAX
      endif
      !                                       limits of the albedo

      ALPHAST(I) = max( ANSMIN, ALPHAST(I) )
      ALPHAST(I) = min( ANSMAX, ALPHAST(I) )


   end do



   !*       13.     EVOLUTION OF SNOW DENSITY
   !                -------------------------


   !                           Density of falling snow

   do I=1,N
      if (SMT(I).gt.0.0) then
         RHOSFALL(I) = 109. + 6.*(T2M(I)-TRPL) +   &
              26.*(U10M(I)**2+V10M(I)**2)**0.25
         RHOSFALL(I) = min(max((RHOSFALL(I)*0.001),RHOMIN), 0.250)
      else
         RHOSFALL(I) = RHOSDEF
      end if
   end do



   !                           Evolution of the snow density depends
   !                           on 3 factors:

   !                           - decrease of density due to fresh new snow
   !                           - increase of density due to aging
   !                           - increase of density due to freezing of
   !                             liquid water in the snow pack


   !                           A) decrease due to fresh new snow

   do I=1,N
      SMX(I)   = max( SMT(I),SR(I)*DT)
      if (SMT(I).gt.0.0) then

         RHOSLT(I) = ( (SMX(I)-SR(I)*DT) * RHOSL(I)  &
              + (SR(I)*DT) * RHOSFALL(I)) / SMX(I)
      else
         ! default
         RHOSLT(I) = RHOSDEF
      end if
   end do


   !                           B) increase due to aging

   do I=1,N
      if (SMT(I).gt.0.0.and.RHOSLT(I).lt.RHOMAX(I)) then
         RHOSLT(I) = (RHOSLT(I)-RHOMAX(I))*exp(-0.01*DT/3600.)  &
              + RHOMAX(I)
      end if
   end do


   !                           C) increase due to freezing

   do I=1,N
      if (SMT(I).gt.0.0) then
         RHOSLT(I) =  ( SMT(I)*RHOSLT(I) + FREEZ_L1(I)*DT*RHOICE ) &
              / ( SMT(I) + FREEZ_L1(I) * DT )
         !                          Make sure within bounds
         RHOSLT(I) = min( RHOICE, RHOSLT(I) )
         RHOSLT(I) = max( RHOMIN, RHOSLT(I) )

      end if
   end do

   !                          Calculate snow depth based on snow density

   do I=1,N
      SNODP(I) = SMT(I)/(RHOSLT(I)*RAUW)


   end do



   !*       14.     UPDATE the PROGNOSTIC VARIABLES
   !                -------------------------------

   do I=1,N
      TSNS(I)     = TSNST(I)
      TSND(I)     = TSNDT(I)
      WL(I)       = WLT(I)
      SM(I)       = SMT(I)
      ALPHAS(I)   = ALPHAST(I)
      RHOSL(I)    = RHOSLT(I)
      RHOSNO(I)   = RHOSLT(I)*RAUW
   end do



   !                  If negligible amount of snow,
   !                  set ALL OUTPUT variables
   !                  to default values and/or zero.


   do I=1,N

      if (VEGH(I).ge.EPSILON_SVS) then

         if (SM(I).lt.CRITSNOWMASS) then
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
         endif

         ! Add trace runoff to runoff variable
         RSNOW(I) = RSNOW(I) + RSNOW_CRITMASS(I)

      else
         ! no high veg, so snow under high veg does not exist... reset all snow variables
         ALPHAS(I)   = ANSMAX
         RHOSL(I)    = RHOSDEF
         RHOSNO(I)   = RHOSLT(I)*RAUW
         TSNS(I)     = 300.0
         TSND(I)     = 300.0
         TAVG(I)     = 300.0
         WL(I)       = 0.0
         SM(I)       = 0.0
         SNODP(I)    = 0.0
         RNET(I)     = 0.0
         HFLUX(I)    = 0.0
         LE(I)       = 0.0
         EFLUX(I)    = 0.0
         RSNOW(I) = 0.0
      endif
   enddo

   return
end subroutine snow_veg1
