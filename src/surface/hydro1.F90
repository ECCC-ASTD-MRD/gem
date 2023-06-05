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

subroutine hydro3(DT, T, TS, WG, W2, WF, WL, WR, WS, ALPHAS, RHOS, SNODP, &
     RAINRATE, SNOWRATE, T2M, U10M, V10M, &
     LEV, LETR, LEG, LES, ZC1, ZC2, C3,  WGEQ, CT, LAI, VEG, D2, &
     WSAT, WFC, PSN, PSNG, PSNV, &
     LEFF, DWATERDT, DSNOWDT, FREEZS, RHOMAX, TST, &
     WGT, W2T, WFT, WLT, WRT, WST, ALPHAST, RHOST, DRAIN, RUNOFF, RSNOW, N)
   use tdpack
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>
   !@Object
   !     Calculates the evolution of the water variables, i.e., the superficial
   !     and deep-soil volumetric water content (wg and w2), the equivalent
   !     liquid water retained in the vegetation canopy (Wr), the equivalent
   !     water of the snow canopy (Ws), and also of the albedo and density of
   !     the snow (i.e., ALBS and RHOS).  Also determine the runoff and drainage
   !     into the soil.
   !@Arguments
   !          - Input -
   ! DT       timestep
   ! TS       surface temperature
   ! WG       superficial volumetric water content
   ! W2       soil volumetric water content
   ! WF       frozen soil water
   ! WL       liquid water in the snow pack
   ! WR       water content retained by the vegetation canopy
   ! WS       equivalent water content of the snow reservoir
   ! ALPHAS   albedo of the snow
   ! RHOS     density of the snow (relative to that of liquid water)
   ! RAINRATE rain rate at the surface
   ! SNOWRATE snow rate at the surface
   ! LEV      latent heat of evaporation over vegetation
   ! LETR     latent heat of evapotranspiration
   ! LEG      latent heat over bare ground
   ! LES      latent heat over snow
   ! ZC1      coefficient for the moisture calculations
   ! ZC2      coefficient for the moisture calculations
   ! C3       coefficient for the drainage
   ! WGEQ     equilibrium volumetric water content
   ! CT       heat capacity coefficient for the total grid
   ! LAI      leaf area index
   ! VEG      fraction of the grid covered by vegetation
   ! D2       rooting depth of the soil
   ! WSAT     volumetric water content at soil saturation
   ! WFC      volumetric water content at field capacity
   ! PSN      fraction of the grid covered by snow
   ! PSNG     fraction of bare ground covered by snow
   ! PSNV     fraction of vegetation covered by snow
   ! LEFF     effective latent heat
   ! DWATERDT exchanged soil water between the frozen and liquid reservoirs
   ! DSNOWDT  exchanged snow water between the frozen and liquid reservoirs
   ! FREEZS   freezing tendency of liquid water in the snow pack
   !          - Input / Output -
   ! TST      new surface temperature
   !           - Output -
   ! WGT      new superficial volumetric water content
   ! W2T      new soil volumetric water content
   ! WFT      new frozen soil water
   ! WLT      new liquid water in the snow pack
   ! WRT      new water content retained by the vegetation canopy
   ! WST      new equivalent water content of the snow reservoir
   ! ALPHAST  new albedo of snow
   ! RHOST    new relative density of snow
   ! DRAIN    drainage
   ! RUNOFF   runoff

   integer N
   real DT, T(N), TS(N), WG(N), W2(N), WR(N), WS(N)
   real WF(N), WL(N)
   real ALPHAS(N), RHOS(N), RAINRATE(N), SNOWRATE(N)
   real T2M(N), U10M(N), V10M(N)
   real SNODP(N)
   real LEV(N), LETR(N), LEG(N), LES(N), ZC1(N), ZC2(N)
   real C3(N), WGEQ(N), CT(N), LAI(N), VEG(N), D2(N)
   real WSAT(N), WFC(N), PSN(N), PSNG(N), PSNV(N)
   real LEFF(N), DWATERDT(N), DSNOWDT(N), FREEZS(N)
   real RHOMAX(N)
   real TST(N)
   real WGT(N), W2T(N), WRT(N), WST(N), ALPHAST(N), RHOST(N)
   real WFT(N), WLT(N)
   real DRAIN(N), RUNOFF(N), RSNOW(N)

   !@Author S. Belair (January 1997)
   !@Revisions
   ! 001      S. Belair (December 1997)
   !                Include the subgrid-scale runoff according
   !                to the Nanjing model (Dumenil and Todini 1992,
   !                Habets and Noilhan 1996).
   !
   ! 002      S. Belair (October 1998)
   !                Partition of precipitation at the surface in
   !                liquid and solid part (that's a temporary
   !                patch).
   !
   ! 003      S. Belair (November 1998)
   !                Small bug for the partition of precip. AND
   !                Refresh properties for snow when there is no snow
   !
   ! 004      S. Belair (December 1998)
   !                Time evolution of the frozen soil water variable
   !
   ! 005      S. Belair (January 1999)
   !                Include new prognostic variable for liquid water
   !                retained in the snow pack
   !
   ! 006      S. Belair (September 1999)
   !                New tendency term for the density of snow related
   !                to freezing of liquid water in the snow pack
   !
   ! 007      S. Belair (January 2000)
   !                Effect of rain on the temperature of the snowpack is
   !                now calculated in "ebudget.ftn"
   !                Calculate minimum albedo of snow "ansmin" and
   !                the maximum density of snow "rhomax"
   !
   ! 008      S. Belair (July 2000)
   !                Bug correction for runoff under the snowpack
   !
   ! 009      B. Bilodeau (January 2001)
   !                Automatic arrays
   !
   ! 010      S. Belair and B. Bilodeau (May 2001)
   !                New density for fresh snow
   !
   !

   include "isbapar.cdk"

   integer :: I

   real :: CRMIN, CRMAX, RHOE, TAUHOUR, RHOICE

   real :: RSNOW_DAY

   real, dimension(n) :: RR, SR, EV, ETR, EG, RUIR, RUIR2, WRMAX, PG, &
        ES, RUIG, RUIP, RUIT, WSX, BEXP, IM, I0IM, PGTMP, RVEG, WLMAX, &
        ANSMIN, RHOSFALL

   !***********************************************************************

   ! THE FOLLOWING SHOULD BE IN A COMMON

   CRMIN    = 0.03
   CRMAX    = 0.10
   RHOE     = 0.20
   TAUHOUR  = 3600.
   RHOICE   = 0.9

   !                   Rainrate and snowrate (in mm/s).

   do I=1,N

      RR(I) = 1000. * RAINRATE(I)
      SR(I) = 1000. * SNOWRATE(I)
 
      ! Freezing rain correction
      if (isba_zr_freeze) then
         if (T(I) <= TRPL .and. TS(I) <= TRPL) then
            SR(I) = SR(I) + RR(I)
            RR(I) = 0.
         endif
      endif

   end do



   ! 1. EVOLUTION OF THE EQUIVALENT WATER CONTENT Wr
   !    --------------------------------------------
   !    evaporation rates

   do I=1,N
      EV(I)  = LEV(I)  / CHLC
      ETR(I) = LETR(I) / CHLC
      EG(I)  = LEG(I)  / LEFF(I)
   end do

   ! evolution of the intercepted water
   ! (if we don t consider the runoff)

   do I=1,N
      WRT(I) = WR(I) - DT &
           * (EV(I) - ETR(I) &
           - VEG(I)*(1-PSNV(I)) * RR(I))
   end do

   ! when Wr < 0, the direct evaporation (i.e., EV-ETR) removes too much
   ! liquid water from the vegetation reservoir.  This is considered as
   ! negative runoff, and it is stocked in RUIR2.
   !
   do I=1,N

      RUIR2(I) = min(0.,WRT(I)/DT)

      ! Wr must be positive

      WRT(I) = max(0., WRT(I))

      ! if Wr > Wrmax, there is runoff

      WRMAX(I) = 0.2*VEG(I)*LAI(I)
      RUIR(I) = max(0., (WRT(I) - WRMAX(I)) / DT )

      !  Wr must be smaller then Wrmax

      WRT(I) = min(WRT(I), WRMAX(I))


      ! Runoff from the vegetation canopy is thus given by RVEG

      RVEG(I) = RUIR(I) + RUIR2(I)

      ! The rate of liquid water reaching the ground is the
      ! precipitation plus the vegetation runoff (we also consider the
      ! negative runoff).

      PG(I) = (1.-VEG(I)) * (1.-PSNG(I)) * RR(I) &
           +               (1.-PSN(I))  * RVEG(I)

   end do

   ! Finally, if Wr < CRITWATER, we force it to 0.0

   do I=1,N
      if (WRT(I).lt.CRITWATER) WRT(I) = 0.0
   end do

   !
   ! 2. EVOLUTION OF THE EQUIVALENT WATER CONTENT Ws
   !    --------------------------------------------
   !    evaporation over snow

   do I=1,N

      ES(I) = LES(I) / (CHLC+CHLF)

      ! evolution of Ws

      WST(I) = WS(I) - DT * (ES(I) - SR(I)) + DSNOWDT(I)

   end do

   !
   ! Check if too much snow has been melted

   do I=1,N
      if (WS(I).lt.0.0) then
         TST(I) = TST(I) - CT(I)*CHLF*DT*WS(I)
         WS(I)  = 0.0
      end if
   end do

   ! if Ws < CRITSNOW, force it to 0.0 and refresh properties
   ! (albedo and density) of the snow pack

   do I=1,N
      if (WST(I).lt.CRITSNOW) then
         WST(I) = 0.0
         ALPHAST(I) = ANSMAX
         ALPHAS(I)  = ANSMAX
         RHOST(I)   = RHOSDEF
         RHOS(I)    = RHOSDEF
      end if
   end do

   ! 3. EVOLUTION OF LIQUID WATER IN THE SNOW PACK
   !    ------------------------------------------
   !    Calculate the maximum liquid water
   !    that can be retained in the snow pack

   do I=1,N
      if (RHOS(I).lt.RHOE) then
         WLMAX(I) = ( &
              CRMIN &
              + (CRMAX-CRMIN) * (RHOE-RHOS(I)) / RHOE &
              ) &
              * WST(I)
      else
         WLMAX(I) = CRMIN * WST(I)
      end if
   end do

   !  Calculate runoff of liquid water from the snow pack

   do I=1,N
      if (WL(I).le.WLMAX(I)) then
         RSNOW(I) = WL(I) / TAUHOUR * &
              exp( WL(I)-WLMAX(I) )
      else
         RSNOW(I) = WLMAX(I) / TAUHOUR &
              + (WL(I)-WLMAX(I)) / DT
      end if
      RSNOW(I) = max( 0.      , RSNOW(I) )
      RSNOW(I) = min( WL(I)/DT, RSNOW(I) )
   end do

   ! Calculate the new values for WL and for the liquid water 
   ! reaching the ground

   do I=1,N
      WLT(I) = WL(I) + PSN(I)*(RR(I)+RVEG(I))*DT &
           - RSNOW(I)*DT - DSNOWDT(I)
      WLT(I) = max( 0., WLT(I) )
   end do

   ! Include the effect of runoff from the snowpack as a source for 
   ! water reaching the ground

   do I=1,N
      PG(I) = PG(I) + RSNOW(I)
   end do

   !  4. SUBGRID-SCALE RUNOFF
   !     --------------------
   !     Preliminary calculations for:
   !     BEXP = exponent b in the formulation
   !     IM   = maximum infiltration capacity
   !     I0IM = fraction i0 / im where i0 is infiltration
   !            capacity of the non-saturated subgrid elements
   !     PGTMP = PG (liquid water fluxes at the surface) with the 
   !             right units (i.e., meters).
   !
   !     Note that the sum of W2+WF is used.
   !     When soil is frozen, greater chance of runoff.
   !
   !     First remove excedent of soil water
   !     (greater than the soil saturation)

   do I=1,N
      RUIT(I) = max( 0., W2(I)+WF(I)-WSAT(I) )
      if (RUIT(I).gt.W2(I)) &
           WF(I) = WF(I)+W2(I)-RUIT(I)
      W2(I) = max( 0., W2(I)-RUIT(I) )
   end do

   ! Subgrid-scale runoff parameterization

   do I=1,N
      BEXP(I)  = 1.
      IM(I)    = (1.+BEXP(I))*WSAT(I)*D2(I)
      I0IM(I)  = 1. &
           - max( 0., ( 1.- (W2(I)+WF(I))/WSAT(I)) ) &
           **( 1./(BEXP(I)+1.) )
      PGTMP(I) = PG(I) * DT / 1000.
      PGTMP(I) = max( PGTMP(I), 0.0 )

      ! Calculate the runoff
      RUNOFF(I) = PGTMP(I) + IM(I)/(BEXP(I)+1) * &
           ( &
           max( 0., ( 1.-I0IM(I)-PGTMP(I)/IM(I) ) )**(BEXP(I)+1.) &
           -  max( 0., ( 1.-I0IM(I) ) )**(BEXP(I)+1.) &
           )
      RUNOFF(I) = max( RUNOFF(I), 0.0 )
      PG(I) = PG(I) - RUNOFF(I)*1000./DT
   end do

   do I=1,N
      RUNOFF(I) = RUNOFF(I)*1000. + RUIT(I)
   end do

   ! 5. EVOLUTION OF THE SUPERFICIAL WATER CONTENT WG
   !     ---------------------------------------------

   do I=1,N
      ! updated values for wg

      WGT(I) = (WG(I) - DT * (ZC1(I) * (EG(I) - PG(I)) &
           / 1000. -  ZC2(I) * WGEQ(I) / 86400.) ) &
           / (1. + DT * ZC2(I) / 86400.)

      !  horizontal runoff when wg > wsat

      RUIG(I) = max(0., WG(I) - WSAT(I) )

   end do


   ! 6. EVOLUTION OF THE DEEP-SOIL WATER CONTENT W2
   !    -------------------------------------------

   do I=1,N

      ! updated values for w2

      W2T(I) = W2(I) - DT*(EG(I) + ETR(I) - PG(I)) &
           / (D2(I) * 1000.)

      ! deep-soil horizontal runoff

      RUIP(I) = RUIG(I) * 10./(D2(I)*1000.)

   end do


   ! 7. DRAINAGE FROM THE DEEP SOIL
   !    ---------------------------
   !    when w2 > wfc, there is drainage
   do I=1,N

      DRAIN(I) = -DT*max(0.,W2T(I)-WFC(I))*C3(I) &
           / (D2(I)*86400.)

      ! the deep-soil volumetric water content w2 is modified consequently

      W2T(I) = W2T(I) + DRAIN(I)

   end do

   ! 8. TIME-EVOLUTION OF THE FROZEN SOIL WATER
   !    ---------------------------------------

   do I=1,N
      WFT(I) = WF(I)  + DWATERDT(I)
      W2T(I) = W2T(I) - DWATERDT(I)
   end do

   ! 9. EVOLUTION OF SNOW ALBEDO
   !    ------------------------
   !    First calculate the minimum value for the snow albedo ANSMIN

   do I=1,N
      RSNOW_DAY = min( RSNOW(I)*24.*3600., 10. )
      ANSMIN(I) = 0.5 - 0.2 * &
           (  1. - cos(RSNOW_DAY*PI / &
           ( 2.*10. ) )               )
   end do

   ! the evolution of the snow albedo differs if there is melting or not

   do I=1,N

      ! when there is melting

      if (WST(I).gt.0.0.and.DSNOWDT(I).lt.0.0) &
           ALPHAST(I) = (ALPHAS(I)-ANSMIN(I))*exp(-0.01*DT/3600.) &
           +  ANSMIN(I) &
           +  SR(I)*DT/WCRN*(ANSMAX-ANSMIN(I))

      if (WST(I).gt.0.0.and.DSNOWDT(I).ge.0.0) &
           ALPHAST(I) = ALPHAS(I) - TODRY*DT/86400. &
           + SR(I)*DT/WCRN*(ANSMAX-ANSMIN(I))

   end do

   ! limits of the albedo

   do I=1,N

      if (WST(I).gt.0.0) then
         ALPHAST(I) = min( ANSMAX, ALPHAST(I) )
         ALPHAST(I) = max( ANSMIN(I), ALPHAST(I) )
      else
         ALPHAST(I) = ANSMAX
      end if

   end do

   ! 10. EVOLUTION OF SNOW DENSITY
   !     -------------------------

   ! Density of falling snow

   do I=1,N
      if (WST(I).gt.0.0) then
         RHOSFALL(I) = 109. + 6.*(T2M(I)-TRPL) + &
              26.*(U10M(I)**2+V10M(I)**2)**0.25
         RHOSFALL(I) = min(max((RHOSFALL(I)*0.001),RHOMIN), 0.250)
      else
         RHOSFALL(I) = RHOSDEF
      end if
   end do

   ! Evolution of the snow density depends on 3 factors:
   !  - decrease of density due to fresh new snow
   !  - increase of density due to aging
   !  - increase of density due to freezing of
   !    liquid water in the snow pack
   !
   !  A) decrease due to fresh new snow

   do I=1,N
      if (WST(I).gt.0.0) then
         WSX(I)   = max( WST(I),SR(I)*DT)
         RHOST(I) = ( (WSX(I)-SR(I)*DT) * RHOS(I) &
              + (SR(I)*DT) * RHOSFALL(I)) / WSX(I)
      end if
   end do

   !  B) increase due to aging

   do I=1,N
      if (WST(I).gt.0.0.and.RHOST(I).lt.RHOMAX(I)) then
         RHOST(I) = (RHOST(I)-RHOMAX(I))*exp(-0.01*DT/3600.) &
              + RHOMAX(I)
      end if
   end do

   !  C) increase due to freezing

   do I=1,N
      if (WST(I).gt.0.0) then
         RHOST(I) = &
              ( WST(I)*RHOST(I) + FREEZS(I)*DT*RHOICE ) &
              / ( WST(I)+FREEZS(I)*DT )
      end if
   end do

   ! lower limit for the snow density

   do I=1,N
      if (WST(I).gt.0.0) then
         RHOST(I) = min( RHOICE, RHOST(I) )
         RHOST(I) = max( RHOMIN, RHOST(I) )
      else
         RHOST(I) = RHOSDEF
      end if
   end do

   do I=1,N
      SNODP(I) = WST(I)/(RHOST(I)*RAUW)
   end do

   ! 11. RUN-OFF
   !     -------

   !  when w2 > wsat, there is runoff

   do I=1,N

      RUIT(I) = max( 0., W2T(I)+WFT(I)-WSAT(I) )
      if (RUIT(I).gt.W2T(I)) &
           WFT(I) = WFT(I)+W2T(I)-RUIT(I)
      W2T(I) = max( 0., W2T(I)-RUIT(I) )

      RUNOFF(I) = RUNOFF(I) + RUIT(I)

   end do

   ! 12. LIMITS FOR WG AND W2
   !     --------------------

   do I=1,N

      WGT(I) = min( WGT(I), WSAT(I) )
      WGT(I) = max( WGT(I), CRITWATER   )

      W2T(I) = min( W2T(I), WSAT(I) )
      WFT(I) = min( WFT(I), WSAT(I) )

      if (W2T(I)+WFT(I).lt.CRITWATER) then
         W2T(I) = CRITWATER
         WFT(I) = 0.0
      end if

      !  we ensure coherence between the frozen and liquid soil water

      if (W2T(I).lt.CRITWATER) then
         WFT(I) = WFT(I) - ( CRITWATER-W2T(I) )
         W2T(I) = CRITWATER
      endif

      WFT(I) = max( WFT(I), 0.0   )

   end do

   return
end subroutine hydro3
