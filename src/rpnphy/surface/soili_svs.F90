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

subroutine SOILI_SVS (WD, &
     WF, SNM, SVM, RHOS, RHOSV, &
     VEGH, VEGL, &
     CGSAT, WSAT, WWILT, BCOEF, &
     CVH, CVL, ALVH, ALVL,  &
     EMISVH, EMISVL, ETG, RGLVH , RGLVL, STOMRVH, STOMRVL,  &
     GAMVH,GAMVL, &
     LAIVH, LAIVL, &
     Z0MVH, Z0MVL, Z0, CLAY, SAND, DECI, EVER,LAID,  &
     Z0HSNOW,Z0MBG, Z0M_TO_Z0H, &
     WTA, CG, PSNGRVL,  &
     Z0H, ALGR, EMGR, PSNVH, PSNVHA,   &
     ALVA, LAIVA, CVPA, EVA, Z0HA, RGLA, STOMRA ,&
     GAMVA, N )
   use tdpack_const, only: PI
   use svs_configs
   use sfc_options, only: read_emis
   implicit none
!!!#include <arch_specific.hf>


   integer N
   real Z0HSNOW, Z0M_TO_Z0H, Z0MBG
   real WD(N,NL_SVS),WF(N,NL_SVS)

   real SNM(N), RHOS(N)
   real RHOSV(N), Z0MVH(N), VEGH(N), VEGL(N), SVM(N)
   real CGSAT(N), WSAT(N,NL_SVS), WWILT(N,NL_SVS), BCOEF(N,NL_SVS)
   real Z0(N)
   real CG(N), WTA(N,indx_svs_ag)
   real PSNGRVL(N)
   real Z0H(N), ALGR(N), CLAY(N), SAND(N)
   real DECI(N), EVER(N), LAID(N)
   real EMGR(N), PSNVH(N), PSNVHA(N),  LAIVH(N)
   real ALVA(N), LAIVA(N), CVPA(N), EVA(N)
   real LAIVL(N), CVH(N), CVL(N), ALVL(N), ALVH(N)
   real EMISVH(N), EMISVL(N), ETG(N), Z0MVL(N), RGLVH(N), RGLVL(N)
   real Z0HA(N), RGLA(N), STOMRA(N), STOMRVH(N), STOMRVL(N)
   real GAMVL(N), GAMVH(N), GAMVA(N)


!Author
!          S. Belair et al. (January 2009)
!Revisions
! 001      Bug fixes: M. Abrahamowicz, S.Z. Husain. N. Gauthier,
!          M. Bani Shahabadi
!
!Object
!     Calculates the coefficients related to the soil (i.e., CG, CT,
!     ) and to the snow canopy
!     (i.e., ZCS, psn,psnvh)
!
!Arguments
!
!
!          - Input -
! WD(NL)   soil volumetric water content in soil layer (NL soil layers)
! WF(NL)   volumetric content of frozen water in the ground
! SNM      equivalent water content of the snow reservoir
! SVM      equivalent water content of the snow-under-vegetation reservoir
! RHOS     relative density of snow
! RHOSV    relative density of snow-under-vegetation
! VEGH     fraction of HIGH vegetation
! VEGL     fraction of LOW vegetation
! CGSAT    soil thermal coefficient at saturation
! WSAT     saturated volumetric moisture content
! WWILT    wilting point volumetric water content
! BCOEF    slope of the retention curve
! LAIVH    Vegetation leaf area index for HIGH vegetation only
! LAIVH    Vegetation leaf area index for LOW  vegetation only
! Z0MVH    Local roughness associated with HIGH vegetation only (no
!          orography)
! Z0MVL    Local roughness associated with LOW vegetation only (no
!          orography)
! CV       heat capacity of the vegetation
! CVH      heat capacity of HIGH vegetation
! CVL      heat capacity of LOW  vegetation
! EMISVH   emissivity of HIGH vegetation
! EMISVL   emissivity of LOW  vegetation
! ETG      emissivity of land surface (no snow) READ AT ENTRY,
!          USED TO STAMP ALL NO SNOW EMISSIVITIES. IF NOT READ,
!          VARIABLE SET TO ZERO !!!!
! RGLVH    param. stomatal resistance for HIGH vegetation
! RGLVL    param. stomatal resistance for LOW  vegetation
! STOMRVH  minim. stomatal resistance for HIGH vegetation
! STOMRVL  minim. stomatal resistance for LOW  vegetation
! GAMVH    stomatal resistance param. for HIGH vegetation
! GAMVL    stomatal resistance param. for LOW  vegetation
! Z0       momentum roughness length (no snow)
! CLAY     percentage of clay in soil (mean of 3 layers)
! SAND     percentage of sand in soil (mean of 3 layers)
! DECI     fraction of high vegetation that is deciduous
! EVER     fraction of high vegetation that is evergreen
! LAID     LAI of deciduous trees
! Z0HSNOW  Constant for thermal roughness of snow
! Z0MBG     Constant momentum roughness for bare ground
! Z0M_TO_Z0H   Conversion factor to convert from momemtum roughness
!          to thermal roughness
!
!           - Output -
! WTA      Weights for SVS surface types as seen from SPACE
! CG       heat capacity of the bare soil
! PSNGRVL  fraction of the bare soil or low veg. covered by snow
! Z0H      agg. roughness length for heat transfer considering snow
! ALGR     albedo of bare ground (soil)
! EMGR     emissivity of bare ground (soil)
! PSNVH    fraction of HIGH vegetation covered by snow
! PSNVHA   fraction of HIGH vegetation covered by snow as seen from
!          the ATMOSPHERE
! ALVA     average vegetation albedo
! LAIVA    average vegetation LAI
! CVPA     average thermal coefficient of vegetation considering effect
!          of LAI
! EVA      average vegetation emissivity
! Z0HA      AVERAGED Local roughness associated with exposed (no snow)
!           vegetation only (also no orography), for heat transfer
! RGLA     average vegetation parameter stomatal resistance
! STOMRA   average minimum stomatal resistance
! GAMVA    average stomatal resistance param.
!
include "isbapar.cdk"


      integer I

      real LAMI, CI, DAY

      real ADRYSAND, AWETSAND, ADRYCLAY, AWETCLAY
      real EDRYSAND, EWETSAND, EDRYCLAY, EWETCLAY

      real, dimension(n) :: a, b, cnoleaf, cva, laivp, lams, lamsv, &
           zcs, zcsv, z0_snow_low


!***********************************************************************





!                                    Define some constants for
!                                    the ice

!                                    NOTE:  these definitions should
!                                           be put in a COMMON

      LAMI   = 2.22
      CI     = 2.106E3
      DAY    = 86400.
!                       Albedo values from literature
      ADRYSAND = 0.35
      AWETSAND = 0.24
      ADRYCLAY = 0.15
      AWETCLAY = 0.08
!                       Emissivity values from van Wijk and Scholte Ubing (1963)
      EDRYSAND = 0.95
      EWETSAND = 0.98
      EDRYCLAY = 0.95
      EWETCLAY = 0.97






!        1.     THE HEAT CAPACITY OF BARE-GROUND
!               --------------------------------

!                          Actually, the 'C' coefficients in
!                          ISBA do not represent heat capacities,
!                          but rather their inverse.  So in the
!                          following formulation, CG is large
!                          when W2 is small, thus leading to small
!                          values for the heat capacity.  In other
!                          words, a small amount of energy will
!                          result in great temperature variations
!                          (for drier soils).

      do I=1,N
        CG(I) = CGSAT(I) * ( WSAT(I,1)/ max(WD(I,1)+WF(I,1),0.001))** &
                      ( 0.5*BCOEF(I,1)/log(10.) )
        CG(I) = min( CG(I), 2.0E-5 )

      end do








!*       2.     THE HEAT CAPACITY OF THE SNOW
!               -----------------------------

!                          First calculate the conductivity of snow

      do I=1,N
        LAMS(I) = LAMI * RHOS(I)**1.88
        ZCS(I) = 2.0 * sqrt( PI/( LAMS(I) * 1000* RHOS(I) *CI*DAY) )

        LAMSV(I) = LAMI * RHOSV(I)**1.88
        ZCSV(I) = 2.0 * sqrt( PI/(LAMSV(I)* 1000*RHOSV(I) *CI*DAY) )

      end do




!*       3.     FRACTIONS OF SNOW COVERAGE
!               --------------------------

      do I=1,N


!                        average snow cover fraction of bare ground and low veg
         if(SNM(I).ge.CRITSNOWMASS ) then

            ! use z0=0.03m for bare ground, 0.1m for low veg
             z0_snow_low(i) = exp (  (  (1-VEGH(I) -VEGL(I)) * log( 0.03) &
                                  +  VEGL(I) * log(0.1) ) / ( 1 - VEGH(I) ) )


             PSNGRVL(I) = min( SNM(I) / (SNM(I) + RHOS(I)* 5000.* z0_snow_low(i) ) , 1.0)

         else

            PSNGRVL(I) = 0.0

         endif




         if(SVM(I).ge.CRITSNOWMASS ) then

!                       SNOW FRACTION AS SEEN FROM THE GROUND


            PSNVH(I)  =  min( SVM(I) / (SVM(I)+ RHOSV(I)*5000.*0.1 ), 1.0)
!                       SNOW FRACTION AS SEEN FROM THE SPACE
!                       NEED TO ACCOUNT FOR SHIELDING OF LEAVES/TREES

            PSNVHA(I) = (EVER(I) * 0.2 + DECI(I) * max(LAI0 - LAID(I), 0.2)) * PSNVH(I)

         else
            PSNVH(I)  = 0.0
            PSNVHA(I) = 0.0
         endif




      end do



!*      4.      FRACTIONS OF SVS TYPES AS SEEN FROM SPACE
!               --------------------------

      call weights_svs(VEGH,VEGL,PSNGRVL,PSNVHA,N,WTA)



!*      5.      AGGREGATED  VEGETATION FIELDS
!              ----------------------------------

!                        Here the snow cover of low vegetation
!                        is taken into account
      do I=1,N

         if((VEGH(I)+VEGL(I)*(1-PSNGRVL(I))).ge.EPSILON_SVS)then

!                        LEAF AREA INDEX

            LAIVA(I) = max(   ( VEGH(I) * LAIVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * LAIVL(I))  &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , EPSILON_SVS)


!                        LEAF AREA INDEX CONSIDERING WOODY PART

            LAIVP(I) =max(  ( VEGH(I)*max(LAIVH(I),0.1)+ VEGL(I)*(1.-PSNGRVL(I))*LAIVL(I)) &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , EPSILON_SVS )

!                        THERMAL COEFFICIENT
!                        -- Impose min. of 1.0E-5 consistent with look-up table in inicove_svs

            CVA(I) = max(  ( VEGH(I) * CVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * CVL(I)) &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) ,  1.0E-5 )

!                        ALBEDO
!                        -- Impose min. of 0.12 consistent with look-up table in inicove_svs

            ALVA(I) = max( ( VEGH(I) * ALVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * ALVL(I)) &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , 0.12 )


!                        PARAMETER STOMATAL RESISTANCE
!                        -- Impose min. of 30. consistent with look-up table in inicove_svs

            RGLA(I) = max( ( VEGH(I) * RGLVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * RGLVL(I)) &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , 30. )

!                        LOCAL ROUGHNESS

            Z0HA(I) = ( VEGH(I)                  * log(Z0MVH(I) * Z0M_TO_Z0H) &
                     + VEGL(I) * (1.-PSNGRVL(I)) * log(Z0MVL(I) * Z0M_TO_Z0H)) &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))

            Z0HA(I) = exp(Z0HA(I))

!                        MINIMUM STOMATAL RESISTANCE
!                        -- Impose min. of 100. consistent with look-up table in inicove_svs

            STOMRA(I) = max( ( VEGH(I) * STOMRVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * STOMRVL(I))  &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , 100.0 )

!                        STOMATAL RESISTANCE PARAMETER

            GAMVA(I) = ( VEGH(I) * GAMVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * GAMVL(I)) &
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))

         else

!                       Set the coefficients parameters to lowest physical values found in look-up table
!                       to have physical values and not bogus values. Use "EPSILON" instead of 0.0 below.


            LAIVA(I) = EPSILON_SVS
            LAIVP(I) = EPSILON_SVS
            CVA(I)   = 1.0E-5
            ALVA(I)  = 0.12
            RGLA(I)  = 30.
            Z0HA(I) = Z0M_TO_Z0H * Z0(I)
            STOMRA(I)= 100.
            GAMVA(I) = EPSILON_SVS
         endif

      end do


!*      6.      THERMAL COEFFICIENT for VEGETATION
!               ----------------------------------
      do I=1,N

        if((VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))).ge.EPSILON_SVS) then


!                       CNOLEAF = thermal coefficient for the leafless
!                       portion of vegetation areas
              CNOLEAF(I) = (VEGL(I)*(1.-PSNGRVL(I)) +  &
                          VEGH(I)*(1.-PSNVH(I))) / CG(I) + PSNGRVL(I)/ZCS(I) &
                         + PSNVH(I)/ZCSV(I)

!                       CVPA = thermal coefficient of vegetation that
!                       considers the effect of LAI


             if(CNOLEAF(I).gt.0.0) then
                CNOLEAF(I) = (VEGL(I)+VEGH(I))/CNOLEAF(I)
             else
                CNOLEAF(I) = 0.0
             endif

             CVPA(I) = min( LAIVP(I),LAI0 ) *  CVA(I)/LAI0  &
                     + max( LAI0 - LAIVP(I), 0.0 ) * CNOLEAF(I)/LAI0

          else


!                       Set the coefficients to lowest physical values found in look-up table
!                       to avoid division by zero

             CVPA(I)=1.0E-5
             CNOLEAF(I)=1.0E-5

          endif

      enddo






!       7.     THERMAL ROUGHNESS LENGTH THAT CONSIDERS THE SNOW
!               ----------------------------------------

!                      The roughness length for heat transfers
!                      is calculated here, in ISB.  For Z0H,
!                      the roughness length used for heat transfers,
!                      only the local roughness is used (associated
!                      with vegetation), i.e., no orography effect

!                      Note that the effect of snow is only on the
!                      thermal roughness length

      do I=1,N
         Z0H(I)=(1-PSNVHA(I)) *      VEGH(I)              * log(Z0MVH(I)*Z0M_TO_Z0H)+   &
                (1-PSNGRVL(I))*      VEGL(I)              * log(Z0MVL(I)*Z0M_TO_Z0H)+   &
                (1-PSNGRVL(I))*(1-VEGH(I)-VEGL(I))        * log(Z0MBG   *Z0M_TO_Z0H)+   &
                (VEGH(I)*PSNVHA(I)+(1-VEGH(I))*PSNGRVL(I))* log(Z0HSNOW)


        Z0H(I)=exp(Z0H(I))

      end do


!        8.     BARE GROUND ALBEDO   EMISSIVITY
!               ---------------------------------------

!                      The bare ground albedo   emissivity for the energy budget
!                      of the bare ground only is calculated here.
!                      It is a rough approximation based on a bi-linear
!                      interpolation of four "extreme" soil values
!                      mainly "dry-sand", "wet-sand", "dry-clay"
!                      and "wet-clay". The albedo   emissivity of these four
!                      categories is estimated from existing
!                      literature. Use superficial volumetric
!                      water content to assess soil wetness as
!                      want *surface* albedo & emissivity.

!                      compute relative texture and soil wetness coefficients A & B

!                       N.B. Even when mapping soil parameters from texture
!                       database unto model layer, use 1st layer texture from
!                       database as is here...

      do I=1,N
!                      A few constraints
         if((CLAY(I)+SAND(I)).gt.0.0) then
            A(I)= SAND(I) / ( CLAY(I) + SAND(I) )
         else
            A(I) = 0.0
         endif
!                      If  superficial soil layer dryer than wilting point
!                      set it to wilting points ...and so get 0.0 for B(I)

         if((WSAT(I,1)-WWILT(I,1)).gt.0.0.and.WD(I,1).ge.WWILT(I,1)) then
            B(I) = ( WD(I,1) - WWILT(I,1) ) / ( WSAT(I,1) - WWILT(I,1))
         else
            B(I) = 0.0
         endif

!                      Added security for coefficients
         if(A(I).gt.1.0) A(I)=1.0
         if(B(I).gt.1.0) B(I)=1.0

      end do

!                      Compute soil albedo
      do I=1,N

         ALGR(I)      = ADRYCLAY * ( 1. - A(I) ) * ( 1. - B(I) ) &
                      + ADRYSAND *        A(I)   * ( 1. - B(I) ) &
                      + AWETCLAY * ( 1. - A(I) ) *        B(I)  &
                      + AWETSAND *        A(I)   *        B(I)

      enddo



!         9.     NO-SNOW EMISSIVITIES
!               ---------------------------------------                     Compute soil albedo & emissivity

!               If use "read-in"  emissivity, then stamp all emissivities with "read at entry value",
!               otherwise calculate individual emissivities


      if ( read_emis ) then

         do I = 1, N
            EMISVH(I)        = ETG(I)
            EMISVL(I)        = ETG(I)
            EMGR  (I)        = ETG(I)
            EVA   (I)        = ETG(I)
         enddo

      else

         do I = 1, N

!                      AGGREGATED VEG. EMISSIVITY  (VEGETATION NOT COVERED BY SNOW)

            if((VEGH(I)+VEGL(I)*(1-PSNGRVL(I))).ge.EPSILON_SVS)  then
               EVA(I) = ( VEGH(I) * EMISVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * EMISVL(I)) &
                    / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))

            else
               EVA(I) = 1.0
            endif

!                      BARE GROUND EMISSIVITY

            EMGR(I)       = EDRYCLAY * ( 1. - A(I) ) * ( 1. - B(I) )  &
                          + EDRYSAND *        A(I)   * ( 1. - B(I) ) &
                          + EWETCLAY * ( 1. - A(I) ) *        B(I) &
                          + EWETSAND *        A(I)   *       B(I)


         enddo

      endif


      return
      end
