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
      SUBROUTINE VEGI_SVS ( RG, T, TVEG, HU, PS, &
           WD , RGL, LAI, LAIH, RSMIN, GAMMA, WWILT, WFC, &   
           SUNCOS, DRZ, D50, D95, PSNGRVL, VEGH, VEGL, RS, SKYVIEW, VTR, &    
           FCD, ACROOT, WRMAX, N  )
!
        use tdpack
        use svs_configs
      implicit none
!!!#include <arch_specific.hf>
!
      INTEGER N 
      REAL WD(N,NL_SVS), FCD(N, NL_SVS)
      REAL RG(N), T(N), HU(N), PS(N), TVEG(N)
      REAL SUNCOS(N), LAIH(N), PSNGRVL(N), VEGH(N), VEGL(N)
      REAL RGL(N), LAI(N), RSMIN(N), GAMMA(N), WWILT(N,NL_SVS)
      REAL WFC(N,NL_SVS), RS(N), SKYVIEW(N), VTR(N), DRZ(N)
      REAL D50(N), D95(N), ACROOT(N,NL_SVS) , WRMAX(N)  
!
!Author
!          S. Belair, M.Abrahamowicz,S.Z.Husain (June 2015)
!Revisions
! 001      V. Fortin (2016) - Active roots !
!
!Object
!
!          Calculates the surface stomatal resistance Rs,
!          and vegetation sky view factor and transmissivity
!
!
!Arguments
!
!
!          - Input -
! RG       solar radiation
! T        low-level temperature of air
! TVEG     skin (surface) temperature of vegetation
! HU       low-level specific humidity of air
! PS       surface pressure
! RGL      AVERAGE constant for the calculation of the stomatal resistance
! LAI      AVERAGE Leaf area index
! LAIH     Leaf area index for high vegetation only
! RSMIN    AVERAGE minimum stomatal resistance
! GAMMA    other constant for RS (AVERAGED)
! WWILT    volumetric water content at the wilting point
! WFC      volumetric water content at the field capacity
! SUNCOS   cosinus of solar angle at LOCAL hour
! DRZ      root-zone depth (rooting depth) [m]
! D50      depth above which 50% of roots are located [m]
! D95      depth above whihch 95% of the roots are located [m]
! WD(NL_SVS)   soil volumetric water content in soil layer (NL_SVS soil layers)
! PSNGRVL  fraction of the bare soil or low veg. covered by snow
! VEGH     fraction of HIGH vegetation
! VEGL     fraction of LOW vegetation
!      
!          - Output -
! RS       Surface or stomatal resistance
! SKYVIEW  Sky view factor for tall/high vegetation
! VTR      (HIGH) Vegetation transmissivity 
! FCD(NL_SVS)  Root fraction within soil layer (NL_SVS soil layers)
! ACROOT(NL_SVS) Active fraction of roots (0-1) in the soil layer
! WRMAX    Max volumetric water content retained on vegetation
!
!
      INTEGER I, K
      REAL CSHAPE, f2_k(nl_svs)
     
      real, dimension(n) :: extinct, f, f1, f2, f3, f4, qsat

!
!***********************************************************************
!
!
!
!
!
       DO I=1,N      
          IF ( (VEGH(I)+VEGL(I)*(1-PSNGRVL(I))).GE.EPSILON_SVS ) THEN
             ! VEGETATION PRESENT


!
!*       1.     THE 'ZF1' FACTOR
!               ---------------
!                      This factor measures the influence
!                      of the photosynthetically active radiation
!
             
             F(I)  = 0.55*2.*RG(I) / (RGL(I)+1.E-6) /  &
                  ( LAI(I)+1.E-6 )
             F1(I) = ( F(I) + RSMIN(I)/5000. ) / ( 1. + F(I) )
!
!
!
!
!*       2.     THE 'ZF2' FACTOR  & ACTIVE ROOT FRACTION
!               ----------------
!                      This factor takes into account the effect
!                      of the water stress on the surface
!                      resistance
!
!                      NOTE that W2 (liquid portion of soil water) is
!                      used here instead of W2+WF.  Thus, when soil water
!                      freezes (W2 --> 0), ZF2 becomes small and the
!                      surface increases increases (no transpiration when
!                      soils are frozen).
!
!           Calculation of root-zone soil moisture
!            
!
            ! Shape parameter CSHAPE
             CSHAPE=LOG10(0.05/0.95)/LOG10(D95(I)/D50(I))       
!        
!                 Root fractions at different depths
!                 The equation has depths in cm, but ratio, so leave in meters                 
             DO K=1,NL_SVS-1
                if(DL_SVS(K).lt.DRZ(I)) then
                   FCD(I,K) = 1./ (1.+ (DL_SVS(K)/D50(I))**CSHAPE)
                else
                   FCD(I,K) = 1.
                endif
             ENDDO
             FCD(I,NL_SVS)=1.
      
             !                 Calculate f2 -- weighted mean using root fraction above
             !                 1.E-5 (dry soil--> large resistance) < f2 < 1.0 (humid soil --> no impact on resistance)
             
             ! f2 for each layer           
             DO K=1,NL_SVS
                ! calculate moisture stress on roots
                f2_k(k) =  min( 1.0,  max( 1.E-5  ,  &
                     max( wd(i,k) - wwilt(i,k) , 0.0) / (wfc(i,k) - wwilt(i,k)) ) )
             ENDDO

             ! root fraction weighted mean
             ! k=1
             f2(i) = f2_k(1) * FCD(I,1)
             DO K=2,NL_SVS
                f2(i) = f2(i) +   f2_k(k) *  ( FCD(I,K) - FCD(I,K-1) ) 
             ENDDO

             ! active fraction of roots
             ! k=2
             acroot(i,1)=f2_k(1)*FCD(I,1)/f2(i)
             DO K=2,NL_SVS
                acroot(i,k)= f2_k(k) *  ( FCD(I,K) - FCD(I,K-1) ) /f2(i)
             ENDDO
    
!
!
!
!*       3.     THE 'ZF3' FACTOR
!               ----------------
!                           This factor represents the effect of
!                           vapor pressure deficit of the atmosphere.
!                           For very humid air, the stomatal resistance
!                           is small, whereas it increases as the
!                           air becomes drier.
!
!
!
             QSAT(I) = FOQST( TVEG(I), PS(I) )
!
             F3(I) = MAX( 1. - GAMMA(I)*( QSAT(I) - HU(I) )*1000. , 1.E-3 )
!
!
!
!*       4.     THE 'ZF4' FACTOR
!               ----------------
!                  This factor introduces an air temperature
!                  dependance on the surface stomatal resistance
!
             F4(I) = MAX( 1.0 - 0.0016*(298.15-T(I))**2, 1.E-3 )             
!
!
!*       5.     THE SURFACE STOMATAL RESISTANCE
!               -------------------------------
!
             RS(I) = RSMIN(I) / ( LAI(I)+1.E-6 ) & 
                  / F1(I) / F2(I) / F3(I) / F4(I)
!
             RS(I) = MIN( RS(I),5000.  )
             RS(I) = MAX( RS(I), 1.E-4 )
       !
!
!*       6.     TRANSMISSIVITY OF CANOPY (HIGH VEG ONLY)
!               ----------------------------------------
!
!
!                 Calculate the extinction coefficient... 
!                 the constant 0.5 is a first approximation
!                 based on Sicart et al. (2004), the calculation
!                 of the coefficient of extinction could/should
!                 be refined according to vegetation type
!
             EXTINCT(I)   =  0.5  / SUNCOS(I)
!
!                 Calculate the transmissivity
!
             VTR(I) = EXP( -1.0 * EXTINCT(I) * LAIH(I) )
!
!
!
!
!*       7.     SKYVIEW FACTOR (HIGH VEG ONLY)
!               ------------------------------
!
!
!                 According to Verseghy et al. (1993), the skyview
!                 factor for NEEDLELEAF is exp(-0.5*LAI) and
!                 exp(-1.5*LAI) for BROADLEAF (and exp(-0.8*LAI) for 
!                 crops). Here as a first approximation, we take the
!                 skyview factor for tall/high vegetation to be exp(-1*LAI).
!
             SKYVIEW(I) = EXP( -1.0 * LAIH(I) )
!
!
!
!*       8.     MAXIMUM VOLUMETRIC WATER CONTENT RETAINED ON VEGETATION (m3/m3)
!               ------------------------------
!
     
          
             WRMAX(I) = 0.2 * LAI(I)

          ELSE
          !  NO VEGETATION
             f1(i)=0.0
             f2(i)=0.0
             f3(i)=0.0
             f4(i)=0.0
             qsat(i) = 0.0
             do K=1,NL_SVS      
                ! set root and active root fraction to zero 
                FCD(I,K) = 0.0 
                acroot(i,k) = 0.0
             enddo
             RS(I) = 5000.
             EXTINCT(I) = 0.0
             VTR(I) = 1.0
             SKYVIEW(i) = 1.0
             WRMAX(I)=EPSILON_SVS  ! To avoid division by zero
          ENDIF
       ENDDO
!
!
      RETURN
      END
