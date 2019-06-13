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
!-------------------------------------- LICENCE END --------------------------

subroutine SOILI2(TS, WG, W2, WF, WS, RHOS, VEG, &
     CGSAT, WSAT, WWILT, BCOEF, C1SAT, C2REF, ACOEF, PCOEF, CV, Z0, Z0VEG, &
     CG, ZC1, ZC2, WGEQ, CT, ZCS, PSN, PSNG, PSNV, PSNZ0, Z0TOT, Z0H, Z0HVEG, N)
   use tdpack_const, only: GRAV, PI
   use sfc_options, only: isba_snow_z0veg
   implicit none
!!!#include <arch_specific.hf>
   !@Object
   !     Calculates the coefficients related to the soil (i.e., CG, CT,
   !     ZC1, ZC2, WGEQ) and to the snow canopy (i.e., ZCS, ps, psng,
   !     psnv, and psnz0)
   !@Arguments
   !          - Input -
   ! TS       surface temperature
   ! WG       superficial volumetric water content
   ! W2       soil volumetric water content
   ! WF       volumetric content of frozen water in the ground
   ! WS       equivalent water content of the snow reservoir
   ! RHOS     density of snow
   ! VEG      fraction of vegetation
   ! CGSAT    soil thermal coefficient at saturation
   ! WSAT     saturated volumetric moisture content
   ! WWILT    wilting point volumetric water content
   ! BCOEF    slope of the retention curve
   ! C1SAT    value of ZC1 at saturation
   ! C2REF    reference value of ZC2
   ! ACOEF    a coefficient for the wgeq formulation
   ! PCOEF    p coefficient for the wgeq formulation
   ! CV       heat capacity of the vegetation
   ! Z0       roughness length
   ! Z0VEG    vegetation roughness length
   !           - Output -
   ! CG       heat capacity of the bare soil
   ! ZC1, ZC2 coefficients for the moisture calculations
   ! WGEQ     equilibrium surface volumetric moisture content
   ! CT       averaged heat capacity of the grid
   ! ZCS      heat capacity of the snow
   ! PSN      fraction of the grid covered by snow
   ! PSNG     fraction of the bare soil covered by snow
   ! PSNV     fraction of the vegetation covered by snow
   ! PSNZ0    fraction of the grid covered by snow for the
   !          calculation of Z0
   ! Z0TOT    total roughness length (including the effect of snow)
   ! Z0H      roughness length for the heat transfers
   ! Z0HVEG   roughness length for the heat transfers (vegetation only)

   integer N
   real TS(N), WG(N), W2(N), WF(N), WS(N), RHOS(N), VEG(N)
   real CGSAT(N), WSAT(N), WWILT(N), BCOEF(N), C1SAT(N)
   real C2REF(N), ACOEF(N), PCOEF(N), CV(N), Z0(N), Z0VEG(N)
   real CG(N), ZC1(N), ZC2(N), WGEQ(N), CT(N), ZCS(N)
   real PSN(N), PSNG(N), PSNV(N), PSNZ0(N), Z0TOT(N)
   real Z0H(N), Z0HVEG(N)

   !@Author S. Belair (January 1997)
   !@Revisions
   ! 001      S. Belair (November 1998)
   !             ZOH calculated from Z0TOT instead of Z0
   ! 002      S. Belair (December 1998)
   !             Correction to the thermal coefficient for bare
   !             soil in order to include the effect of frozen water
   ! 003      S. Belair (January 2000)
   !             change the roughness length of snow (now 0.005)
   !             and change the formulation for "psng", i.e., the
   !             fraction of bare soil covered by snow
   ! 004      B. Bilodeau and S. Belair (February 2000)
   !             Formulation for WGEQ modified to avoid underflows
   ! 005      B. Bilodeau (January 2001)
   !             Automatic arrays
   ! 006      S. Belair and L.-P. Crevier (February 2001)
   !             Put roughness length of snow equal to Z0
   !             Cap Z0H at 0.2 m

   integer i

   real LAMI, CI, DAY, WCRN, BETAS, Z0SNOW
   real LAMTYP, LAMW, CW
   real, dimension(n) :: LAMS, CW1MAX, ETA, WMAX, SIGSQUARE, BETA, &
        X, Z0_PSNV

   !***********************************************************************
   !  Define some constants for the ice
   !  NOTE:  these definitions should be put in a COMMON

   LAMI   = 2.22
   CI     = 2.106E3
   LAMTYP = 1.5
   LAMW   = 0.56
   CW     = 4.2E3
   DAY    = 86400.
   WCRN   = 10.0
   BETAS  = 0.408
   Z0SNOW = 0.005

   ! 1.     THE HEAT CAPACITY OF BARE-GROUND
   !         --------------------------------
   !  Actually, the 'C' coefficients in ISBA do not represent heat capacities,
   !  but rather their inverse.  So in the following formulation, CG is large
   !  when W2 is small, thus leading to small values for the heat capacity.
   !  In other words, a small amount of energy will result in great
   !  temperature variations (for drier soils).

   do I=1,N
      CG(I) = CGSAT(I) * ( WSAT(I)/ max(W2(I)+WF(I),0.001))** &
           ( 0.5*BCOEF(I)/log(10.) )

   end do

   ! 2. CORRECTION OF THE HEAT CAPACITY OF GROUND DUETO FROZEN WATER
   !    ---------------------------------------------
   !    When soil water is completely or partially frozen,
   !    the thermal coefficient CG has to be
   !    modified to include the impact of the
   !    thermal conductivity and capacity of ice.
   !
   !************************************************
   !   ATTENTION: cette boucle est temporairement hors d'usage
   !              puisqu'il y a des points ou des valeurs < 0
   !              sont mis a la puissance 0.5
   !***********************************************
   !***      DO I=1,N
   !***        LAMC(I)  = (4.*PI) / (CG(I)*CG(I)*DAY)
   !***        CTYP(I)  = LAMC(I) / LAMTYP
   !***        CGCOR(I) = 1.
   !***     1           + WF(I)*WF(I)*(LAMI-LAMW)*(CI-CW)/LAMC(I)
   !***     1           + WF(I)*( (LAMI-LAMW)/LAMTYP + (CI-CW)/CTYP(I) )
   !***        CG(I)    = CG(I) / CGCOR(I)**0.5
   !***      END DO

   !# Cg can not be greater than 2.E-5
   do I=1,N
      CG(I) = min( CG(I), 2.0E-5 )
   end do

   ! 3. THE HEAT CAPACITY OF THE SNOW
   !    -----------------------------
   !    First calculate the conductivity of snow
   do I=1,N
      LAMS(I) = LAMI * RHOS(I)**1.885
      ZCS(I) = 2.0 * sqrt( PI/(1000.*LAMS(I)* &
           RHOS(I)*CI*DAY) )
   end do


   !4. FRACTIONS OF SNOW COVERAGE
   !    -------------------------
   if (isba_snow_z0veg) then
      Z0_PSNV(:) = Z0VEG(:)
   else
      Z0_PSNV(:) = Z0(:)
   endif
   do I=1,N

      !# fraction of ground covered by snow
      if (WS(I).lt.WCRN) then
         PSNG(I) = WS(I) / WCRN
      else
         PSNG(I) = 1.
      end if

      !# fraction of vegetation covered by snow
      PSNV(I) = WS(I) / (WS(I)+RHOS(I)*5000.*Z0_PSNV(I))

      !# fraction of the grid covered by snow
      PSN(I) = ( 1-VEG(I) )*PSNG(I) + VEG(I)*PSNV(I)

      !# fraction of the grid covered by snow for the calculation of z0
      PSNZ0(I) = WS(I) / ( WS(I)+WCRN+BETAS*GRAV*Z0(I) )

   end do


   !5. GRID-AVERAGED HEAT CAPACITY
   !   ---------------------------
   !  (with contribution from the ground, vegetation, and snow areas)
   !
   do I=1,N
      CT(I) = 1. / ( (1.-VEG(I))*(1.-PSNG(I)) / CG(I) &
           +         VEG(I)   *(1.-PSNV(I)) / CV(I) &
           +                PSN(I)          / ZCS(I)  )
   end do

   ! 6. COEFFICIENT ZC1
   !    --------------
   !    The coefficient ZC1 is calculated two
   !    different ways depending on the humidity of the soil

   do I=1,N
      ! First situation:  humid soil
      ! Then the calculation follows eq. (19) of Noilhan and Planton(1989)

      if (WG(I).gt.WWILT(I)) then

         ZC1(I) = C1SAT(I)*(WSAT(I)/WG(I))**(0.5*BCOEF(I)+1)

      else

         ! Second situation: dry soil
         ! We use the Gaussian formulation of
         ! Giordanni (1993) and Braud et al. (1993)

         CW1MAX(I)    = ( 1.19*WWILT(I)-5.09 )*0.01*TS(I) &
              + (-1.464*WWILT(I)+17.86)
         CW1MAX(I)    = max( CW1MAX(I)  ,  0.02 )
         ETA(I)       =   (-1.815E-2*TS(I)+6.41)*WWILT(I) &
              + (6.5E-3*TS(I)-1.4)
         WMAX(I)      = ETA(I)*WWILT(I)
         SIGSQUARE(I) = - ( WMAX(I)*WMAX(I) )  / &
              ( 2.0*log( 0.01/CW1MAX(I) ) )
         BETA(I)      = 1. - 2.*VEG(I)*( 1.-VEG(I) )

         ZC1(I) = CW1MAX(I)*BETA(I) &
              *exp( -(WG(I)-WMAX(I))*(WG(I)-WMAX(I)) / &
              (2.*SIGSQUARE(I)) )

         ! C1 must be divided by an arbitrarynormalization depth of 0.1 m
         ZC1(I) = 10.0*ZC1(I)

      end if

   end do

   !* 7. COEFFICIENT ZC2
   !     ---------------
   do I=1,N
      ZC2(I) = C2REF(I)*W2(I) / &
           max( (WSAT(I)-W2(I)), 0.001)
   end do

   ! 8. EQUILIBRIUM VOLUMETRIC WATER CONTENT WGEQ
   !    -----------------------------------------

   do I=1,N
      X(I) = (W2(I)+WF(I)) / WSAT(I)
      WGEQ(I) = X(I) - ACOEF(I)*X(I)**PCOEF(I)
      WGEQ(I) = WGEQ(I)*WSAT(I)
   end do

   ! 9. ROUGHNESS LENGTH THAT CONSIDERS THE SNOW
   !    ----------------------------------------
   !    CAREFUL:  Here we calculate the roughness length for heat 
   !              (it is a factor times the roughness length for momentum).
   !
   !    ATTENTION ... ATTENTION ... ATTENTION ...
   !    The roughness lenght for heat transfers
   !    is now calculated here, in ISBA, using Z0TOT
   !    which includes the contributions from the
   !    vegetation, orography, and snow.  In the future,
   !    ZOH should be determine another way, in PHYDEBU.

   do I=1,N
      !        Z0TOT(I) = (1.-PSNZ0(I)) * LOG(Z0(I))
      !     1           +     PSNZ0(I)  * LOG(Z0SNOW)
      !        Z0TOT(I) = EXP( Z0TOT(I) )
      Z0TOT(I) = Z0(I)
   end do

   do I=1,N
      ! Z0H(I) = MAX(   0.1 * Z0TOT(I)  ,  1.E-5   )
      Z0H(I) = min(   0.2 * Z0TOT(I)  ,  0.2   )
      Z0HVEG(I) = min(   0.2 * Z0VEG(I)  ,  0.2   )
   end do

   return
end subroutine SOILI2
