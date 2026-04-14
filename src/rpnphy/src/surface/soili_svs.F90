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
module soili_svs_mod
  implicit none
  public
contains
      SUBROUTINE SOILI_SVS (WD, &
           WF, SNM, SVM, RHOS, RHOSV, & 
           VEGH, VEGL, &  
           CGSAT, WSAT, WWILT, BCOEF, & 
           CVH, CVL, ALVH, ALVL,  & 
           EMISVH, EMISVL, ETG, RGLVH , RGLVL, STOMRVH, STOMRVL,  & 
           GAMVH,GAMVL, &  
           LAIVH, LAIVL, &   
           Z0MVH, Z0MVL, Z0, CLAY, SAND, DECI, EVER,LAID,  &  
           CONDDRY, CONDMINFAC, CONDSLD, &  
           WTA, CG, PSNGRVL,  & 
           Z0H, ALGR, EMGR, PSNVH, PSNVHA,   &
           ALVA, LAIVA, CVPA, EVA, Z0HA, Z0MVG, RGLA, STOMRA ,&  
           GAMVA, SOILHCAPZ, SOILCONDZ, SIG_TOPO, N )
         !
        use tdpack_const, only: PI
        use svs_configs
        use sfc_options, only: read_emis, svs_urban_params, soil_cond, svs_snowfrac_ground, svs_cg_depth, svs_cg_ice
     implicit none
!!!#include <arch_specific.hf>

!
      INTEGER N 
      REAL WD(N,NL_SVS),WF(N,NL_SVS)

      REAL SNM(N), RHOS(N)
      REAL RHOSV(N), Z0MVH(N), VEGH(N), VEGL(N), SVM(N)
      REAL CGSAT(N,NL_SVS), WSAT(N,NL_SVS), WWILT(N,NL_SVS), BCOEF(N,NL_SVS)
      REAL Z0(N)
      REAL CG(N), WTA(N,svs_tilesp1)
      REAL PSNGRVL(N)
      REAL Z0H(N), ALGR(N), CLAY(N), SAND(N), SIG_TOPO(N)
      REAL DECI(N), EVER(N), LAID(N)
      REAL EMGR(N), PSNVH(N), PSNVHA(N),  LAIVH(N)
      REAL ALVA(N), LAIVA(N), CVPA(N), EVA(N)
      REAL LAIVL(N), CVH(N), CVL(N), ALVL(N), ALVH(N)
      REAL EMISVH(N), EMISVL(N), ETG(N), Z0MVL(N), RGLVH(N), RGLVL(N)
      REAL Z0HA(N), Z0MVG(N), RGLA(N), STOMRA(N), STOMRVH(N), STOMRVL(N)
      REAL GAMVL(N), GAMVH(N), GAMVA(N)
      REAL SOILHCAPZ(N,NL_SVS), SOILCONDZ(N,NL_SVS)
      REAL CONDDRY(N,NL_SVS), CONDMINFAC(N,NL_SVS), CONDSLD(N,NL_SVS)
      
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
! SIG_TOPO Standard deviation of subgrid topo [m]
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
! CLAY     percentage of clay of surface soil layer 
! SAND     percentage of sand of surface soil layer   
! DECI     fraction of high vegetation that is deciduous
! EVER     fraction of high vegetation that is evergreen
! LAID     LAI of deciduous trees
! CONDSLD  Soil thermal conductivity
! CONDDRY  Dry thermal conductivity
! CONDMINFAC  Factor applied to the thermal conductivity of minerals (only when optsoilcond = 1)
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
! Z0HA     AVERAGED Local roughness associated with exposed (no snow)
      !           vegetation only (also no orography), for heat transfer
! Z0MVG    AVERAGED momentum roughness for exposed vegetation (no snow)    
! RGLA     average vegetation parameter stomatal resistance
! STOMRA   average minimum stomatal resistance
! GAMVA    average stomatal resistance param. 
! SOILHCAPZ soil heat capacity
! SOILCONDZ soil thermal conductivity

include "isbapar.cdk"

!
      INTEGER I, K
!    
      REAL LOG_CONDI, LOG_CONDW, XF, XU, WORK1, WORK2, WORK3, CONDSAT
      REAL LAM_ZERO, k_min, k_air, k_ice 
      REAL SATDEG, KERSTEN
      REAL C_ICE
!
      real, dimension(n) :: a, b, cnoleaf, cva, laivp, lams, lamsv, &
           zcs, zcsv, z0_snow_low, sd_opn, sd_for

      real, dimension(n,nl_svs) :: cgk

      REAL :: CVAMIN = 1.0E-5

      IF (SVS_URBAN_PARAMS) THEN
         CVAMIN = 0.3E-5   ! matches value of CVDAT(21) reset in inicover_svs.F90
      ENDIF
!
!***********************************************************************
!
!
!
!        1.     THE HEAT CAPACITY OF BARE-GROUND
!               --------------------------------
!
!                          Actually, the 'C' coefficients in
!                          ISBA do not represent heat capacities,
!                          but rather their inverse.  So in the
!                          following formulation, CG is large
!                          when W2 is small, thus leading to small
!                          values for the heat capacity.  In other
!                          words, a small amount of energy will
!                          result in great temperature variations
!                          (for drier soils).
!
!
      IF (svs_cg_ice .EQ. 'BELAIR2003') THEN
        DO I=1,N
          DO K=1,NL_SVS
            CGK(I,K) = CGSAT(I,K) * ( WSAT(I,K)/ MAX(WD(I,K)+WF(I,K),0.001))** &
                    ( 0.5*BCOEF(I,K)/LOG(10.) )
            CGK(I,K) = MIN( CGK(I,K), 2.0E-5 )
          END DO
        END DO
      ELSE IF (svs_cg_ice .EQ. 'BOONE2000') THEN
        ! Boone et al. (2000), eq. B3
        C_ICE = 2 * (PI / (LAMI*CICE*RHOI*DAY))**0.5
        DO I=1,N
          DO K=1,NL_SVS
            ! Boone et al. (2000), eq. B2
            CGK(I,K) = CGSAT(I,K) * ((WSAT(I,K)-WF(I,K)) / MAX(WD(I,K),0.001))** &
                    ( 0.5*BCOEF(I,K)/LOG(10.) )
            CGK(I,K) = MIN( CGK(I,K), 2.0E-5 )
            CGK(I,K) = (1-WF(I,K)) * CGK(I,K) + WF(I,K) * C_ICE
          END DO
        END DO
      ELSE
         ! Unknown option for svs_cg_ice
         STOP
      END IF
!
!     Average value of Cg computed using weights for each layer
!     that decrease exponentially with a constant characteristic length
!     defined by option svs_cg_depth. Weights are constant
!     and computed in svs_configs.
      DO I=1,n
        CG(I) = CGK(I,1)
        IF (svs_cg_depth.GT.0.) THEN
          DO K=2,NL_SVS
            CG(I) = CG(I)*(1.-WEIGHT_CG(K)) + CGK(I,K)*WEIGHT_CG(K)
          END DO
        END IF
      END DO
!
!
!
!
!
!
!
!
!*       2.     THE HEAT CAPACITY OF THE SNOW
!               -----------------------------
!
!                          First calculate the conductivity of snow
!
      DO I=1,N
        LAMS(I) = LAMI * RHOS(I)**1.88
        ZCS(I) = 2.0 * SQRT( PI/( LAMS(I) * 1000* RHOS(I) *CICE*DAY) )
!
        LAMSV(I) = LAMI * RHOSV(I)**1.88
        ZCSV(I) = 2.0 * SQRT( PI/(LAMSV(I)* 1000*RHOSV(I) *CICE*DAY) )
!
      END DO    
!
!
!
!
!*       3.     FRACTIONS OF SNOW COVERAGE
!               --------------------------
!
      DO I=1,N


!                        average snow cover fraction of bare ground and low veg
         IF(SNM(I).GE.CRITSNOWMASS ) THEN


             IF(SVS_SNOWFRAC_GROUND=='NIL') THEN   
               ! use z0=0.03m for bare ground, 0.1m for low veg
               z0_snow_low(i) = exp (  (  (1-VEGH(I) -VEGL(I)) * log( 0.03) &
                                  +  VEGL(I) * log(0.1) ) / ( 1 - VEGH(I) ) )
             

               PSNGRVL(I) = MIN( SNM(I) / (SNM(I) + RHOS(I)* 5000.* z0_snow_low(i) ) , 1.0)

             ELSE IF(SVS_SNOWFRAC_GROUND=='LA23') THEN 

               ! Snow depth in open terrain [m]      
               SD_OPN(I) = SNM(I) / (RHOS(I)*1000.)

               ! Roughness parameter accounting for the fraction of low veg. 
               Z0_SNOW_LOW(I) = EXP (  (  (1-VEGH(I) -VEGL(I)) * LOG(Z0BG_LA23) &
                                  +  VEGL(I) * LOG(Z0LV_LA23)) / ( 1 - VEGH(I) ) )

               ! Snow cover fraction accounting for sugbrid topography as in Lalande et al. (2023)
               PSNGRVL(I) = MIN(1.0, TANH(SD_OPN(I)/                                          &
                          (2.5 * Z0_SNOW_LOW(I)  * (RHOS(I)*1000./RHON_LA23)**MFAC_LA23)      & 
                            +  BETA_LA23 * SIG_TOPO(I) * (RHOS(I)*1000./RHON_LA23)**NFAC_LA23)) 

             ELSE IF(SVS_SNOWFRAC_GROUND=='AR25') THEN

               ! Snow depth in open terrain [m]  
               SD_OPN(I) = SNM(I) / (RHOS(I)*1000.)

               ! Snow cover fraction optimized for a resolution of 2.5 km taken from Abolafia-Rosenzweig et al. (2025)
               PSNGRVL(I) = MIN(1.0, TANH(SD_OPN(I)/( SCF_AR25  * (RHOS(I)*1000./RHON_AR25)**MFAC_AR25)))

             ENDIF               

         ELSE

             PSNGRVL(I) = 0.0
           
         ENDIF
            
       

         IF(SVM(I).GE.CRITSNOWMASS ) THEN
!
!                       SNOW FRACTION AS SEEN FROM THE GROUND
!
            IF(SVS_SNOWFRAC_GROUND=='NIL') THEN   

               ! Snow cover fraction seen from the ground for snow below high vegetation [-]
               PSNVH(I)  =  MIN( SVM(I) / (SVM(I)+ RHOSV(I)*5000.*0.1 ), 1.0)

            ELSE IF(SVS_SNOWFRAC_GROUND=='LA23') THEN 

               ! Snow depth below high vegetation [m]                   
               SD_FOR(I) = SVM(I) / (RHOSV(I)*1000.) 

               ! Snow cover fraction seen from the ground for snow below high vegetation [-]
               ! Effect of subgrid topo are not taken into account in the forest
               PSNVH(I) = MIN(1.0, TANH(SD_FOR(I)/(2.5 * Z0HV_LA23  * (RHOSV(I)*1000./RHON_LA23)**MFAC_LA23)))

             ELSE IF(SVS_SNOWFRAC_GROUND=='AR25') THEN 

               ! Snow depth below high vegetation [m]   
               SD_FOR(I) = SVM(I) / (RHOSV(I)*1000.) 

               ! Snow cover fraction seen from the ground for snow below high vegetation [-]
               PSNVH(I) = MIN(1.0, TANH(SD_FOR(I)/( SCF_AR25  * (RHOSV(I)*1000./RHON_AR25)**MFAC_AR25)))
               
            ENDIF

!                       SNOW FRACTION AS SEEN FROM THE SPACE
!                       NEED TO ACCOUNT FOR SHIELDING OF LEAVES/TREES
!
            IF(SVS_SNOWFRAC_GROUND=='NIL') THEN
               ! Original formulation used in SVS     
               PSNVHA(I) = (EVER(I) * 0.2 + DECI(I) * MAX(LAI0 - LAID(I), 0.2)) * PSNVH(I)

            ELSE IF(SVS_SNOWFRAC_GROUND=='LA23' .OR. SVS_SNOWFRAC_GROUND=='AR25') THEN 
               ! Revised formulation to be compatible with LA23 and AR25     
               PSNVHA(I) = (EVER(I) * 0.15 + DECI(I) * MAX(0.8 - LAID(I), 0.15)) * PSNVH(I)
            ENDIF

         ELSE
            PSNVH(I)  = 0.0
            PSNVHA(I) = 0.0
         ENDIF



!        
      END DO
!

!
!*      4.      FRACTIONS OF SVS TYPES AS SEEN FROM SPACE
!               --------------------------

      call weights_svs(VEGH,VEGL,PSNGRVL,PSNVHA,N,WTA)


!
!*      5.      AGGREGATED  VEGETATION FIELDS
!              ----------------------------------
!
!                        Here the snow cover of low vegetation
!                        is taken into account
      DO I=1,N

         IF((VEGH(I)+VEGL(I)*(1-PSNGRVL(I))).GE.EPSILON_SVS)THEN
!
!                        LEAF AREA INDEX
!
            LAIVA(I) = MAX(   ( VEGH(I) * LAIVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * LAIVL(I))  & 
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , EPSILON_SVS)

!
!                        LEAF AREA INDEX CONSIDERING WOODY PART
!
            LAIVP(I) =MAX(  ( VEGH(I)*MAX(LAIVH(I),0.1)+ VEGL(I)*(1.-PSNGRVL(I))*LAIVL(I)) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , EPSILON_SVS )
!
!                        THERMAL COEFFICIENT 
!                        -- Impose min. of 1.0E-5 consistent with look-up table in inicove_svs
!
            CVA(I) = MAX(  ( VEGH(I) * CVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * CVL(I)) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , CVAMIN )
!
!                        ALBEDO
!                        -- Impose min. of 0.12 consistent with look-up table in inicove_svs
!
            ALVA(I) = MAX( ( VEGH(I) * ALVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * ALVL(I)) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , 0.12 ) 
!
! 
!                        PARAMETER STOMATAL RESISTANCE
!                        -- Impose min. of 30. consistent with look-up table in inicove_svs
!
            RGLA(I) = MAX( ( VEGH(I) * RGLVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * RGLVL(I)) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , 30. ) 
!
!                        LOCAL VEG. MOMENTUM ROUGNESS
            Z0MVG(I)= ( VEGH(I)                   * LOG(Z0MVH(I)) &
                      + VEGL(I) * (1.-PSNGRVL(I)) * LOG(Z0MVL(I))) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))

            Z0MVG(I)= EXP(Z0MVG(I))

!                        LOCAL VEG. THERMAL ROUGNESS            
!
            Z0HA(I) = ( VEGH(I)                  * LOG(Z0MVH(I) * Z0M_TO_Z0H) &
                     + VEGL(I) * (1.-PSNGRVL(I)) * LOG(Z0MVL(I) * Z0M_TO_Z0H)) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))

            Z0HA(I) = EXP(Z0HA(I))
!
!                        MINIMUM STOMATAL RESISTANCE
!                        -- Impose min. of 100. consistent with look-up table in inicove_svs
!
            STOMRA(I) = MAX( ( VEGH(I) * STOMRVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * STOMRVL(I))  &    
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))) , 100.0 )
!
!                        STOMATAL RESISTANCE PARAMETER
!
            GAMVA(I) = ( VEGH(I) * GAMVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * GAMVL(I)) &     
                              / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))
!
         ELSE

!                       Set the coefficients parameters to lowest physical values found in look-up table
!                       to have physical values and not bogus values. Use "EPSILON" instead of 0.0 below.
!
!               
            LAIVA(I) = EPSILON_SVS
            LAIVP(I) = EPSILON_SVS
            CVA(I)   = CVAMIN
            ALVA(I)  = 0.12
            RGLA(I)  = 30.
            Z0HA(I) = Z0M_TO_Z0H * Z0(I)  
            STOMRA(I)= 100.
            GAMVA(I) = EPSILON_SVS
         ENDIF
!                        
      END DO
!
!
!*      6.      THERMAL COEFFICIENT for VEGETATION
!               ----------------------------------
      DO I=1,N
!
        IF((VEGH(I)+VEGL(I)*(1.-PSNGRVL(I))).GE.EPSILON_SVS) THEN
!
!
!                       CNOLEAF = thermal coefficient for the leafless   
!                       portion of vegetation areas
              CNOLEAF(I) = (VEGL(I)*(1.-PSNGRVL(I)) +  &
                          VEGH(I)*(1.-PSNVH(I))) / CG(I) + PSNGRVL(I)/ZCS(I) & 
                         + PSNVH(I)/ZCSV(I)
!
!                       CVPA = thermal coefficient of vegetation that
!                       considers the effect of LAI
!
!
             IF(CNOLEAF(I).gt.0.0) THEN
                CNOLEAF(I) = (VEGL(I)+VEGH(I))/CNOLEAF(I)
             ELSE
                CNOLEAF(I) = 0.0
             ENDIF

             CVPA(I) = MIN( LAIVP(I),LAI0 ) *  CVA(I)/LAI0  & 
                     + MAX( LAI0 - LAIVP(I), 0.0 ) * CNOLEAF(I)/LAI0 
             
          ELSE

!         
!                       Set the coefficients to lowest physical values found in look-up table 
!                       to avoid division by zero
!
             CVPA(I) = CVAMIN
             CNOLEAF(I) = CVAMIN
!     
          ENDIF
!      
      ENDDO
!
!
!
!
!
!
!       7.     THERMAL ROUGHNESS LENGTH THAT CONSIDERS THE SNOW
!               ----------------------------------------
!
!                      The roughness length for heat transfers
!                      is calculated here, in ISB.  For Z0H, 
!                      the roughness length used for heat transfers, 
!                      only the local roughness is used (associated
!                      with vegetation), i.e., no orography effect
!
!                      Note that the effect of snow is only on the 
!                      thermal roughness length
!
      DO I=1,N         
         Z0H(I)=(1-PSNVHA(I)) *      VEGH(I)              * LOG(Z0MVH(I)*Z0M_TO_Z0H)+   &
                (1-PSNGRVL(I))*      VEGL(I)              * LOG(Z0MVL(I)*Z0M_TO_Z0H)+   &
                (1-PSNGRVL(I))*(1-VEGH(I)-VEGL(I))        * LOG(Z0MBG   *Z0M_TO_Z0H)+   &
                (VEGH(I)*PSNVHA(I)+(1-VEGH(I))*PSNGRVL(I))* LOG(Z0HSNOW)
 

        Z0H(I)=EXP(Z0H(I))

      END DO
!
!
!        8.     BARE GROUND ALBEDO   EMISSIVITY
!               ---------------------------------------
!
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
!
!                      compute relative texture and soil wetness coefficients A & B
!
!                       N.B. Even when mapping soil parameters from texture 
!                       database unto model layer, use 1st layer texture from 
!                       database as is here... 
!
      DO I=1,N
!        A few constraints
!	     Take the texture of the surface soil layer 
         IF((CLAY(I)+SAND(I)).gt.0.0) THEN         
            A(I)= SAND(I) / ( CLAY(I) + SAND(I) ) 
         ELSE         
            A(I) = 0.0           
         ENDIF
!                      If  superficial soil layer dryer than wilting point
!                      set it to wilting points ...and so get 0.0 for B(I)
!
         IF((WSAT(I,1)-WWILT(I,1)).gt.0.0.and.WD(I,1).ge.WWILT(I,1)) THEN         
            B(I) = ( WD(I,1) - WWILT(I,1) ) / ( WSAT(I,1) - WWILT(I,1))
         ELSE
            B(I) = 0.0 
         ENDIF
! 
!                      Added security for coefficients
         IF(A(I).gt.1.0) A(I)=1.0
         IF(B(I).gt.1.0) B(I)=1.0
!
      END DO

!                      Compute soil albedo
      DO I=1,N
!
         ALGR(I)      = ADRYCLAY * ( 1. - A(I) ) * ( 1. - B(I) ) & 
                      + ADRYSAND *        A(I)   * ( 1. - B(I) ) &   
                      + AWETCLAY * ( 1. - A(I) ) *        B(I)  &   
                      + AWETSAND *        A(I)   *        B(I)   

      ENDDO

!
!      
!         9.     NO-SNOW EMISSIVITIES
!               ---------------------------------------                     Compute soil albedo & emissivity 
!
!               If use "read-in"  emissivity, then stamp all emissivities with "read at entry value", 
!               otherwise calculate individual emissivities


      if ( read_emis ) then

         DO I = 1, N
            EMISVH(I)        = ETG(I)
            EMISVL(I)        = ETG(I)
            EMGR  (I)        = ETG(I)
            EVA   (I)        = ETG(I)
         ENDDO
  
      else

         DO I = 1, N
!
!                      AGGREGATED VEG. EMISSIVITY  (VEGETATION NOT COVERED BY SNOW)
!
            IF((VEGH(I)+VEGL(I)*(1-PSNGRVL(I))).GE.EPSILON_SVS)  THEN
               EVA(I) = ( VEGH(I) * EMISVH(I) + VEGL(I) * (1.-PSNGRVL(I)) * EMISVL(I)) &     
                    / (VEGH(I)+VEGL(I)*(1.-PSNGRVL(I)))
!
            ELSE
               EVA(I) = 1.0
            ENDIF
!
!                      BARE GROUND EMISSIVITY
!
            EMGR(I)       = EDRYCLAY * ( 1. - A(I) ) * ( 1. - B(I) )  &  
                          + EDRYSAND *        A(I)   * ( 1. - B(I) ) &   
                          + EWETCLAY * ( 1. - A(I) ) *        B(I) &    
                          + EWETSAND *        A(I)   *       B(I) 
            !

         ENDDO

      endif

!
!
!       10.      SOIL THERMAL PROPERTIES PROFILE (either PL98 or Tian2016)
!               ------------------------------------------
!
!
       LOG_CONDI = LOG(LAMI)
       LOG_CONDW = LOG(LAMW)
       DO I=1,N

          DO K=1,NL_SVS

              IF (soil_cond == 'PL1998') THEN     !use the model from Peters-Lidard et al. (1998) for frozen soil that involve the Kersten number (Johanssen, 1975)
    
!                Kersten parameter for thermal conductivity
                 XF = WF(I,K) / (WF(I,K) + MAX(WD(I,K),0.001))
                 XU = (1.0-XF) * WSAT(I,K)
    
                 WORK1   = LOG(CONDSLD(I,K))*(1.0-WSAT(I,K))
                 WORK2   = LOG_CONDI*(WSAT(I,K)-XU)
                 WORK3   = LOG_CONDW*XU
                 CONDSAT = EXP(WORK1+WORK2+WORK3)
                 SATDEG  = MAX(0.1, (WF(I,K) + WD(I,K))/WSAT(I,K))  ! degree of saturation
                 SATDEG  = MIN(1.0,SATDEG)
                 KERSTEN  = LOG10(SATDEG) + 1.0               ! Kersten number
    
!                Put in a smooth transition from thawed to frozen soils:
!                simply linearly weight Kersten number by frozen fraction
!                in soil:
                 KERSTEN  = (1.0-XF)*KERSTEN + XF *SATDEG
    
!                Thermal conductivity
                 SOILCONDZ(I,K) = KERSTEN*(CONDSAT-CONDDRY(I,K)) + CONDDRY(I,K)

              ELSE IF (soil_cond == 'TIAN2016') THEN     !use the physical model from Tian et al. (2016) [https://doi.org/10.1111/ejss.12366]
                  IF (WD(I,K) > 0.2*WSAT(I,K)) THEN
                      LAM_ZERO = LAMW
                  ELSE
                      IF (WF(I,K) > 0.5*WSAT(I,K)) THEN
                          LAM_ZERO = LAMI
                      ELSE
                          LAM_ZERO = 0.025 !thermal conductivity of air
                      ENDIF
                  ENDIF

                  k_min = (2/3)*(1+(CONDSLD(I,K)/LAM_ZERO-1)*CONDMINFAC(I,K))**(-1) + &
                          (1/3)*(1+(CONDSLD(I,K)/LAM_ZERO-1)*(1-2*CONDMINFAC(I,K)))**(-1)

                  k_ice = (2/3)*(1+(LAMI/LAM_ZERO-1)*0.333*(1-WF(I,K)/WSAT(I,K)))**(-1) + (1/3)*(1+(LAMI/LAM_ZERO-1)*(1-2*0.333*(1-WF(I,K)/WSAT(I,K))))**(-1)
                  
                  k_air = (2/3)*(1+(0.025/LAM_ZERO-1)*0.333*(1-(WSAT(I,K)-WD(I,K)-WF(I,K))/WSAT(I,K)))**(-1) + &
                          (1/3)*(1+(0.025/LAM_ZERO-1)*(1-2*0.333*(1-(WSAT(I,K)-WD(I,K)-WF(I,K))/WSAT(I,K))))**(-1)

                  SOILCONDZ(I,K) = (WD(I,K)*LAMW + k_ice*WF(I,K)*LAMI + k_air*(WSAT(I,K) - WD(I,K) - WF(I,K))*0.025 + k_min*(1 - WSAT(I,K))*CONDSLD(I,K)) / &
                                  (WD(I,K) +k_ice*WF(I,K)+k_air*(WSAT(I,K)-WD(I,K)-WF(I,K))+k_min*(1-WSAT(I,K)))
              ENDIF

!             Heat capacity (J m-3 K-1)
              SOILHCAPZ(I,K) = (1. - WSAT(I,K)) * 2700. * 733. + WD(I,K) * CWAT * RHOW &
                             + WF(I,K) * CICE * RHOI

          END DO

       END DO

      RETURN
    END SUBROUTINE SOILI_SVS
  end module soili_svs_mod
