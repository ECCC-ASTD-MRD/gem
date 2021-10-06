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
!-------------------------------------- LICENCE END --------------------------

subroutine ebudget4(T, TS, T2, W2, WF, WL, &
     WS, DT, ALPHAS, RAINRATE, &
     RG, ALVG, ALBT, EMIST, EMSVC, RAT, THETAA, HU, PS, RHOA, &
     VEG, HRSURF, HV, DEL, RESA, RS, CT, &
     CG, ZCS, PSN, PSNV, PSNG, WSAT, D2, SNODP, &
     TST, T2T, RNET, HFLUX, LE, LEG, LEV, &
     LES, LER, LETR, GFLUX, EFLUX, &
     LEFF, DWATERDT, DSNOWDT, FREEZS, RHOMAX, &
     MELTS_TOT, MELTS_RN, FTEMP, FVAP, N)
   use tdpack
   use sfc_options, only: rad_off, atm_external, isba_melting_fix, &
        snow_emiss, &
        snow_emiss_const, isba_soil_emiss, isba_soil_emiss_const
   implicit none
!!!#include <arch_specific.hf>
   !@Object
   !
   !     Calculates the evolution of the surface and deep-soil temperature
   !     (i.e., Ts and T2), as well as all the surface fluxes.
   !@Arguments
   !          - Input -
   ! T         surface air temperature
   ! TS        surface temperature
   ! T2        mean surface temperature
   ! W2        soil water
   ! WF        frozen soil water
   ! WL        liquid water in the snow pack
   ! WS        equivalent water content of the snow reservoir
   ! DT        timestep
   ! ALPHAS    albedo of snow
   ! RAINRATE  rainrate
   ! RG        global radiation (downward solar)
   ! ALVG      surface albedo associated with vegetation type
   ! ALBT      total surface albedo (snow + vegetation)
   ! RAT       atmospheric radiation incident on the ground (NIR)
   ! THETAA    air potential temperature at the lowest level
   ! HU        specific humidity of air at the lowest level
   ! PS        surface pressure
   ! RHOA      air density near the surface
   ! VEG       fraction of the grid covered by vegetation
   ! HRSURF    relative humidity of the surface
   ! HV        Halstead coefficient (relative humidity of veg. canopy)
   ! DEL       portion of the leaves covered by water
   ! RESA      aerodynamic surface resistance
   ! RS        stomatal resistance
   ! CT        total heat capacity coefficient
   ! ZCS       heat capacity coefficient for the snow canopy
   ! PSN       fraction of the grid covered by snow
   ! PSNV      fraction of vegetation covered by snow
   ! PSNG      fraction of bare ground covered by snow
   ! WSAT      volumetric water content at saturation
   ! D2        soil depth
   !           - Output -
   ! TST       new surface temperature
   ! T2T       new mean surface temperature
   ! RNET      net radiation
   ! HFLUX     sensible heat flux
   ! LE        latent heat flux
   ! LEG       latent heat flux over bare ground
   ! LEV       latent heat flux over vegetation
   ! LES       latent heat flux over snow
   ! LER       direct latent heat flux from vegetation leaves
   ! LETR      evapotranspiration latent heat flux
   ! GFLUX     ground flux
   ! EFLUX     water vapor flux
   ! DWATERDT  net tendency of melting-freezing of soil water
   ! DSNOWDT   net tendency of melting-freezing of snow water
   ! FREEZS    tendency of freezing water in the snow pack
   ! MELTS_TOT total snow melting (accumulator)
   ! MELTS_RN  snow melting due to rain (accumulator)

   integer N
   real T(N), TS(N), T2(N), WS(N), DT, ALPHAS(N)
   real RAINRATE(N)
   real W2(N), WF(N), WL(N), WSAT(N), D2(N), SNODP(N)
   real RG(N), ALVG(N), ALBT(N), RAT(N), THETAA(N)
   real HU(N), PS(N), RHOA(N), VEG(N)
   real HRSURF(N), HV(N), DEL(N), RESA(N), RS(N), CT(N)
   real CG(N), ZCS(N), PSN(N), PSNV(N), PSNG(N)
   real TST(N), T2T(N), RNET(N), HFLUX(N), LE(N)
   real LEG(N), LEV(N), LES(N), LER(N), LETR(N), GFLUX(N)
   real EFLUX(N), EMIST(N), EMSVC(N)
   real FTEMP(N), FVAP(N)
   real LEFF(N), DWATERDT(N), DSNOWDT(N), FREEZS(N)
   real RHOMAX(N), MELTS_TOT(N), MELTS_RN(N)

   !@Author S. Belair (January 1997)
   !@Revisions
   ! 001      S. Belair (November 1998)
   !             Use physics package thermodynamic functions

   ! 002      S. Belair (December 1998)
   !             Correction to the latent heat coefficients due to
   !             soil water freezing.
   !             Calculation of the FREEZG and MELTG tendencies for
   !             soil water freezing and melting.
   !             Tendencies for the surface temperature TS.

   ! 003      S. Belair (January 1999)
   !             Tendencies for melting and freezing of snow

   ! 004      B. Bilodeau (December 1999)
   !             real downward radiation as an argument
   !             different albedos for vegetation and snow as arguments

   ! 005      S. Belair (January 2000)
   !             diagnostic equation for the maximum density of snow
   !             effect of incident rain on the snowpack

   ! 006      B. Bilodeau (January 2001)
   !             Automatic arrays

   ! 007      S. Belair (Fall 2003)
   !             include the effect of soil freeze/thaw

   ! 008      B. Bilodeau (April 2007)
   !             add snow melting diagnostics

   ! 009      Olympics Team (2010)
   !             impose condition on snow melt impact on sfc temp.

   ! 010      N. Gauthier (April 2013)
   !             correction to units in computation of max snow density

   ! 011      M. Abrahamowicz (August 2013)
   !             Do not consider RADIA=NIL key when FLUVERT='SURFACE'
   !             For FLUVERT='SURFACE only:
   !                - Rewrite check on TST adjustment related to snow melt
   !                - Bug fix in snow melt term related to TST adjustment
   !                - Bug fix in snow melt term related to rain melt
   !@Notes
   !***  METHOD
   !*    ------

   !     1- find the grid-averaged albedo, emissivity, and roughness length
   !     2- compute the za, zb, and zc terms involved in the numerical
   !        resolution of the equations for Ts and T2.
   !     3- find Ts(t) and T2(t).
   !     4- derive the surface fluxes.*

   real,parameter :: PETIT = 1.E-7

   integer :: i
   real :: KCOEF, RHOW
   real :: MLTRAIN, RAIN1, RAIN2
   real :: TEMPO
   real, dimension(n) :: ZQSAT, ZDQSAT, ZQSATT, RORA, A, B, C, TN, &
        ZHV, FREEZFRAC, FREEZG, MELTG, TNMT0, MELTS, WORK, EG, ES, EV, &
        FMLTRAIN, EMISSS, EMISSN

   !***********************************************************************

   ! THE FOLLOWING SHOULD BE PUT IN A COMMON COMDECK

   RHOW   = 1000.
   KCOEF  = 1.E-6

   if (isba_melting_fix) KCOEF = 1.55E-7

   RAIN1  = 2.8e-8
   RAIN2  = 2.8e-7

   ! 1. GRID-AVERAGED ALBEDO, EMISSIVITY, AND ROUGHNESS LENGTH
   !    -----------------------------------------------------
   !  (considering snow surfaces)

   select case (snow_emiss)
   case DEFAULT
      emissn(:) = snow_emiss_const
   end select
   select case (isba_soil_emiss)
   case ('CLIMATO')
      emisss(:) = emsvc(:)
   case DEFAULT
      emisss(:) = isba_soil_emiss_const
   end select
   if (rad_off .and. .not. atm_external) then
      !Do not include radiative forcings
      ALBT = 1.
      EMIST = 0.
   else
      do I=1,N

         ! For the albedo
         ALBT(I) = ( 1.-PSN(I) )*ALVG(I) + PSN(I)*ALPHAS(I)

         ! For the emissivity

         ! ATTENTION ... ATTENTION ... ATTENTION ...
         ! The surface emissivity for bare soil and
         ! vegetation is now fixed at 0.95 ... this
         ! could be changed in the future.  We could
         ! use an analysis for example.  Or determine
         ! the emissivity from the vegetation or soil
         ! types.

         EMIST(I) = ( 1.-PSN(I) ) * EMISSS(I) + &
              PSN(I)  * EMISSN(I)
      end do
   endif


   ! 2. LATENT HEAT COEFFICIENTS - CORRECTION DUE TO FREEZING
   !    AND MELTING OF SOIL WATER
   !    -----------------------------------------------------
   ! Using the fraction of frozen water
   ! in the soil, calculate the "effective"
   ! latent heat of evaporation/sublimation

   do I=1,N
      FREEZFRAC(I) = WF(I) / (W2(I)+WF(I)+PETIT)
      LEFF(I)      = FREEZFRAC(I)      * (CHLC+CHLF) &
           + (1.-FREEZFRAC(I)) *  CHLC
   end do

   ! 3. COEFFICIENTS FOR THE TIME INTEGRATION OF  TS
   !     --------------------------------------------

   ! Thermodynamic functions
   do I=1,N

      ZQSAT(I)  = FOQST( TS(I),PS(I) )
      ZDQSAT(I) = FODQS( ZQSAT(I),TS(I) )

   end do

   ! function zrsra
   do I=1,N
      RORA(I) = RHOA(I) / RESA(I)
   end do

   ! terms za, zb, and zc for the calculation of ts(t)

   do I=1,N

      A(I) = 1. / DT + CT(I) * (4. * EMIST(I) * STEFAN &
           * (TS(I)**3) +  RORA(I) * ZDQSAT(I) * &
           ( CHLC*VEG(I)*(1-PSNV(I))*HV(I) &
           + LEFF(I)*(1.-VEG(I))*(1.-PSNG(I))*HRSURF(I) &
           + (CHLC+CHLF)*PSN(I) )+ RORA(I) * CPD) &
           + 2. * PI / 86400.

      B(I) = 1. / DT + CT(I) * (3. * EMIST(I) * STEFAN &
           * (TS(I)** 3) + RORA(I) * ZDQSAT(I) * &
           ( CHLC*VEG(I)*(1-PSNV(I))*HV(I) &
           + LEFF(I)*(1.-VEG(I))*(1.-PSNG(I))*HRSURF(I) &
           + (CHLC+CHLF)*PSN(I) ) )

      C(I) = 2. * PI * T2(I) / 86400. + CT(I) &
           * (RORA(I) * CPD * THETAA(I) + RG(I) &
           * (1. - ALBT(I)) + EMIST(I)*RAT(I) - RORA(I) &
           * ( CHLC*VEG(I)*(1.-PSNV(I))*HV(I) &
           * (ZQSAT(I)-HU(I)) &
           + LEFF(I)*(1.-VEG(I))*(1.-PSNG(I)) &
           * (HRSURF(I)*ZQSAT(I)-HU(I)) &
           + (CHLC+CHLF)*PSN(I) * (ZQSAT(I)-HU(I)) ) )

   end do

   do I=1,N
      TST(I) = ( TS(I)*B(I) + C(I) ) / A(I)
   end do

   ! 4. MELTING AND FREEZING TENDENCIES OF SNOW
   !    ---------------------------------------

   ! Calculate the temperature TN (include cover effect of vegetation)
   do I=1,N
      TN(I) = (1.-VEG(I))*TST(I) + VEG(I)*T2(I)
   end do


   ! Common portion of the MELTS and FREEZS equations
   do I=1,N
      WORK(I) = PSN(I) * (TN(I)-TRPL) / ( ZCS(I)*CHLF*DT )
   end do

   ! MELTS and FREEZS tendencies
   ! Also calculate the maximum snow density
   do I=1,N
      if (WORK(I).lt.0.) then
         MELTS(I)  = 0.0
         FREEZS(I) = min( -WORK(I), WL(I)/DT )
         !   RHOMAX(I) = 450. - 20470. / (SNODP(I)+PETIT) * &
         !  ( 1.-EXP(-SNODP(I)/67.3))
         RHOMAX(I) = 450. - 204.7 / (SNODP(I)+PETIT) * &
              ( 1.-exp(-SNODP(I)/0.673))
         RHOMAX(I) = 0.001 * RHOMAX(I)
      else
         MELTS(I)  = min( WORK(I) , WS(I)/DT )
         FREEZS(I) = 0.0
         !   RHOMAX(I) = 600. - 20470. / (SNODP(I)+PETIT) * &
         !  ( 1.-EXP(-SNODP(I)/67.3))
         RHOMAX(I) = 600. - 204.7 / (SNODP(I)+PETIT) * &
              ( 1.-exp(-SNODP(I)/0.673))
         RHOMAX(I) = 0.001 * RHOMAX(I)
      end if
   end do

   ! new temperature Ts(t) after melting
   do I=1,N
      TEMPO=CT(I) * CHLF * (FREEZS(I)-MELTS(I)) * DT
      if ( (tempo < 0.) .and. (TST(I) > TRPL ) ) then

         ! Trying to COOL DOWN  surface temp. that is
         ! initially above zero, below zero by melting snow ...
         ! don't let it happen

         TST(I) = max(TRPL,TST(I)+TEMPO)

      else if ( (tempo > 0.) .and. (TST(I) < TRPL ) ) then

         ! Trying to WARM UP surface temp. that is initially
         ! below zero, above zero by freezing liquid water
         ! in  snow ... don't let it happen

         TST(I) = min(TRPL,TST(I)+TEMPO)

      else if ( (tempo > 0.) .and. (TST(I) >= TRPL ) ) then

         ! Trying to WARM UP surface temp. that is initially
         ! above zero by freezing liquid water
         ! in  snow ... don't let it happen if  isba_no_warm_sn_freez=.TRUE.

!!$         if ( .not.  isba_no_warm_sn_freez ) then
!!$
!!$            TST(I) = TST(I)+TEMPO
!!$
!!$         endif
         continue

      else

            TST(I) = TST(I)+TEMPO

      endif
   end do

   !  4. EFFECT OF RAIN ON TEMPERATURE OF SNOW PACK
   !     ------------------------------------------

   !  When rain is falling on snow,
   !  melting is accelerated due to heat transfers between the
   !  incident rain and the snow pack (since rain is
   !  usually warmer then the snow pack).

   !  It is hypothesized that the temperature of falling water is the
   !  same as that of air at the lowest atmospheric level.

   do I=1,N
      if (RAINRATE(I).lt.RAIN1) then
         FMLTRAIN(I) = 0.
      else if (RAINRATE(I).gt.RAIN2) then
         FMLTRAIN(I) = 1.
      else
         FMLTRAIN(I) = ( RAINRATE(I) - RAIN1 ) / &
              ( RAIN2       - RAIN1 )
      end if
   end do

   do I=1,N
      if (T(I).gt.TRPL.and.WS(I).gt.0.0.and. &
           RAINRATE(I).gt.0.) then

         MLTRAIN = max(  ( T(I)-max(TST(I),TRPL) ) / ( 2.*ZCS(I)*CHLF*DT )  , 0. )
         MELTS(I) = MELTS(I) + FMLTRAIN(I) * MLTRAIN
         MELTS_RN (I) = MELTS_RN (I) + MLTRAIN*DT

      end if
   end do

   ! Melting-Freezing tendency for the WS and WL reservoirs

   do I=1,N
      DSNOWDT(I) = ( FREEZS(I)-MELTS(I) ) * DT
      ! DIAGNOSTIC : accumulated melting (in kg/m2 or mm)
      MELTS_TOT(I) = MELTS_TOT(I) + MELTS(I)*DT
   end do

   !  5. FREEZING AND MELTING TENDENCIES FOR SOIL WATER
   !     ----------------------------------------------
   !  Recalculate TN with the new surfacetemperature
   !  (include cover effect of vegetation canopy AND snow pack)

   do I=1,N
      TN(I) = (1-PSN(I)) * ( (1-VEG(I))*TST(I)+VEG(I)*T2(I) ) &
           +    PSN(I)  *    T2(I)

      !    Calculate the K coefficient in the
      !    FREEZG and MELTG flux terms

      TNMT0(I) = TN(I)  - TRPL
   end do

   ! Calculate the freezing and melting tendencies

   if (isba_melting_fix) then

      do I=1,N
         if (TNMT0(I).ge.0.) then
            FREEZG(I) = 0.
            MELTG(I)  = (KCOEF/D2(I)) * TNMT0(I)
         else
            FREEZG(I) = -(KCOEF/D2(I)) * (W2(I)/WSAT(I)) * TNMT0(I)
            MELTG(I)  = 0.
         end if
         DWATERDT(I) = (FREEZG(I)-MELTG(I))*DT
      end do

   else

      do I=1,N
         if (TNMT0(I).ge.0.) then
            FREEZG(I) = 0.
            MELTG(I)  = KCOEF * TNMT0(I)
         else
            FREEZG(I) = - KCOEF * (W2(I)/WSAT(I)) * TNMT0(I)
            MELTG(I)  = 0.
         end if
         DWATERDT(I) = (FREEZG(I)-MELTG(I))*DT
      end do

   endif

   ! Make sure we don't remove more liquid
   ! or solid water than there is in the soil

   do I=1,N
      if (DWATERDT(I).ge.0.) then
         DWATERDT(I) = min( W2(I), DWATERDT(I) )
      else
         DWATERDT(I) = max( -WF(I), DWATERDT(I) )
      end if
   end do



   ! Update the surface temperature TS

   !     DO I=1,N
   ! TST(I) = TST(I) + CG(I)*CHLF*DWATERDT(I)
   !     END DO

   ! 6. T2 AT TIME 'T+DT'
   !    -----------------

   do I=1,N
      T2T(I) = (T2(I) + DT*TST(I)/86400.) / &
           (1.+DT/86400.)

      ! Include the effect of soil freeze/thaw
      T2T(I) = T2T(I) &
           + RHOW*CG(I)*CHLF*D2(I)*DWATERDT(I)*DT/86400.
   end do


   ! 7. FLUX CALCULATIONS
   !    -----------------


   do I=1,N
      ! recalculate the qsat function
      ZQSATT(I) = FOQST(  TST(I),  PS(I)   )

      ! net radiation
      RNET(I) = (1. - ALBT(I)) * RG(I) + EMIST(I) * &
           (RAT(I) - STEFAN * (TST(I)** 4))

      ! sensible heat flux
      HFLUX(I) = RHOA(I) * CPD * (TST(I) - THETAA(I)) / RESA(I)
      FTEMP(I) = (TST(I) - THETAA(I)) / RESA(I)

      ! latent heat of evaporation from the ground
      LEG(I) = RHOA(I) * LEFF(I) * (1.-VEG(I))*(1.-PSNG(I)) &
           * (HRSURF(I)* ZQSATT(I) - HU(I)) / RESA(I)

      EG(I) = (1.-VEG(I))*(1.-PSNG(I)) &
           * (HRSURF(I)* ZQSATT(I) - HU(I)) / RESA(I)

      ! latent heat of evaporation from the snow canopy
      LES(I) = RHOA(I) * (CHLC+CHLF) * PSN(I) * &
           (ZQSATT(I) - HU(I)) / RESA(I)

      ES(I) =  PSN(I) * (ZQSATT(I) - HU(I)) / RESA(I)

      ! latent heat of evaporation from vegetation
      LEV(I) = RHOA(I) * CHLC * VEG(I)*(1-PSNV(I)) &
           * HV(I) * (ZQSATT(I) - HU(I)) / RESA(I)

      EV(I) =  VEG(I)*(1-PSNV(I)) &
           * HV(I) * (ZQSATT(I) - HU(I)) / RESA(I)

      ! latent heat of evapotranspiration
      ZHV(I) = max(0., sign(1.,ZQSAT(I) - HU(I)))
      LETR(I) = ZHV(I) * (1. - DEL(I)) * RHOA(I) &
           * CHLC * VEG(I)*(1-PSNV(I)) &
           * (ZQSATT(I)  - HU(I)) / (RESA(I) + RS(I))


      ! latent heat of direct evaporation
      LER(I)  = LEV(I) - LETR(I)

      ! total latent heat of evaporation
      LE(I) = LEG(I) + LEV(I) + LES(I)

      ! total water vapor flux
      EFLUX(I) = EG(I) + EV(I) + ES(I)
      FVAP(I)  = EFLUX(I)

      ! heat flux into the ground
      GFLUX(I) = RNET(I) - HFLUX(I) - LE(I)
   end do

   return
end subroutine ebudget4
