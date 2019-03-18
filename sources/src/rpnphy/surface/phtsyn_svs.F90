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
      SUBROUTINE PHTSYN_SVS ( LAI_NCLASS, VEGFRAC, &
               TCAN,  PRESSG,  RESAVG,  QA, QSWV1,  WD, &
               FCD, COSZS, WFC, WWILT, MASKLAT, &
!--------------------------INPUTS ABOVE AND OUTPUTS BELOW --------------
               FCANC, AILCG,  RC,CO2I1,  AVG_GWSOL, NCLASS, N)
!
!
        use svs_configs
      implicit none
#include <arch_specific.hf>

      INTEGER N, NCLASS
      REAL LAI_NCLASS(N,NCLASS),VEGFRAC(N,NCLASS)
      REAL WD(N,NL_SVS),FCD(N,NL_SVS),WFC(N,NL_SVS), WWILT(N,NL_SVS)
   
      INTEGER KK,ICC,IG, IC, L2MAX, NN, MASKLAT(N), GEXP
      PARAMETER (KK=12)  ! PRODUCT OF CLASS PFTs AND L2MAX (4 x 3 = 12)
      PARAMETER (ICC=9)  ! 9 types of vegetation CCMA land surface model
      PARAMETER (IG=3)  ! 
      PARAMETER (IC=4)  ! Number of plant functional types PFTS:
      ! NEEDLE LEAF , BROAD LEAF , CROPS, GRASSES
      PARAMETER (L2MAX=3)  !
      PARAMETER (GEXP=2) ! !     EXPONENT FOR SOIL MOISTURE STRESS. FOR SN EQUAL TO 1, PHOTOSYNTHESIS
!     DECREASES LINEARLY WITH SOIL MOISTURE, AND OF COURSE NON-LINEARLY
!     FOR VALUES HIGHER THAN 1. WHEN SN IS ABOUT 10, PHOTOSYNTHESIS DOES
!     NOT START DECREASING UNTIL ABOUT SOIL MOISTURE IS HALF WAY BETWEEN
!     WILTING POINT AND FIELD CAPACITY.
!!!!!!
      INTEGER KK_TO_ICC(ICC) 
      REAL FC_TEST(N), SIGMA(N),  TGAMMA(N), KC(N), KO(N), IPAR(N), GB(N)
      REAL RH(N), VPD(N), O2_CONC(N), CO2A(N), USEBB(ICC)
!
      REAL BETA_WSOL(N,NL_SVS), G_WSOL(N,NL_SVS)
      REAL AVG_GWSOL(N), VMAXC(N,ICC) 
      REAL VMUNS1(N,ICC), VMUNS2(N,ICC), VMUNS3(N,ICC), VMUNS(N,ICC), VM(N,ICC)
      REAL CO2I(N,ICC), PREV_CO2I(N,ICC), FPAR(N,ICC), JC(N,ICC),  JC1(N,ICC)
      REAL JC2(N,ICC), JC3(N,ICC), JE(N,ICC), JE1(N,ICC), JE2(N,ICC), JS(N,ICC)
      REAL A_VEG(N,ICC), RC_VEG(N,ICC), GCTU(N,ICC), GCMIN(N,ICC), GCMAX(N,ICC)
      REAL VPD_TERM(N,ICC), CO2LS(N,ICC), GC(N,ICC)
!!!!!!----------------------------------------------------------------------

      INTEGER NOL2PFTS(IC)
      INTEGER I, J, ICOUNT, K1, K2, M, REQITER, K
      INTEGER PS_COUP, IT_COUNT
!
      REAL  TCAN(N), PRESSG(N), CFLUX(N), RESAVG(N), QA(N)
      REAL  QSWV(N),  QSWV1(N)
      REAL  COSZS(N),  RC(N),  FC(N)
!
!
      REAL  RML_VEG(N,ICC), AN_VEG(N,ICC)
      REAL  CO2I1(N,ICC), AILCG(N,ICC), FCANC(N,ICC)
!!
      REAL   Q10, Q10_FUNCN, Q10_FUNC, Q10_FUNCD
      REAL   TEMP_B, TEMP_C, TEMP_R, TEMP_Q1, TEMP_Q2,  TEMP_JP
      REAL   BETA1,  BETA2,   GAMMA_W, GAMMA_M, CO2IMAX
!
      REAL VPD0(KK), OMEGA(KK),VMAX(KK)
      REAL TUP(KK), TLOW(KK), KN(KK)
      REAL INICO2I(KK), ALPHA(KK), RMLCOEFF(KK), BB(KK), MM(KK)
      REAL CO2CC, DELTA_CO2, N_EFFECT
      INTEGER ISC4(KK)

!
!Author
!          S. Zhang (March 2013)
!Revisions
! 001      M. Abrahamowicz (May 2015)
!            Re-ogarnize the subroutine, pass NCLASS as argument, rename variables
!
!Object
!          The initial version was from CTEM of CCCma
!    
!
!Method
!
!
! ???
!
!Arguments
!
!          - Input/Output -
!
!
!          - Input -
!
! LAI_NCLASS  Leaf Area Index for the NCLASS classes
! VEGFRAC     Fraction of vegetation for NCLASS classes 
! TCAN        skin (surface) temperature of vegetation
! PRESSG      surface pressure
! RESAVG      aerodynamical surface resistance for vegetation 
! QA          low-level specific humidity of air
! QSWV1       net radiation (aggregated over snow, vegetation, soil)
! WD(NL)      soil volumetric water content in soil layer (NL soil layers)
! FCD(NL)     Root fraction within soil layer (NL soil layers)
! COSZS       cosinus of solar angle at LOCAL HOUR
! WFC(NL)     volumetric water content at the field capacity (NL soil layers)
! WWILT(NL)   volumetric water content at the wilting point (NL soil layers)
! MASKLAT     mask = 1.0 if ABS(latitude) <= 50.0deg, 0.0 otherwise
! NCLASS      number of classes for natural covers 
!          - Output -
!
! FCANC      vegetation fraction for CTEM veg types 1-9
! AILCG      Mean LAI for CTEM veg types 1-9
! RC         stomatal resistance 
! CO2I1      intercellular CO2 concentration 
! AVG_GWSOL  average soil moisture stress term G (common to all CTEM veg types)
      REAL  TFREZ, STD_PRESS, ZERO
      REAL  CA, CB, EA, EASAT
      LOGICAL SMOOTH, MIN2, MIN3

   

!***********************************************************************
!        Vegetation Classes correspondance
!     -----------------------------------------------------------------
!
!     CONSTANTS AND PARAMETERS USED IN THE PHOTOSYNTHESIS MODEL. ALSO
!     NOTE THE STRUCTURE OF VECTORS WHICH CLEARLY SHOWS THE CLASS
!     PFTs (ALONG ROWS) AND CTEM SUB-PFTs (ALONG COLUMNS)
!
!     NEEDLE LEAF |  EVG       DCD       ---
!     BROAD LEAF  |  EVG   DCD-CLD   DCD-DRY
!     CROPS       |   C3        C4       ---
!     GRASSES     |   C3        C4       ---
!
!
!      1    NEEDLE LEAF EVERYGREEN            
!      2    NEEDLE LEAF DECIDUOUS
!      3    BROAD  LEAF EVERYGREEN
!      4    BROAD  LEAF COLD DECIDUOUS
!      5    BROAD  LEAF DRY  DECIDUOUS
!      6    C3  CROPS
!      7    C4  CROPS
!      8    C3  GREEN GRASSES
!      9    C4  GREEN GRASSES
!
!
!                           the geophysical fields determined from
!                           vegetation are done so using the following
!                           classification:
!
!    SVS_PFTs     Vegetation type                        CTEM_PFTs
!     =====       ===============
!       1         (salt) water
!       2         ice
!       3         inland lake
!       4         evergreen needleleaf trees              1    C3 
!       5         evergreen broadleaf trees               3    C3
!       6         deciduous needleleaf trees              2    C3
!       7         deciduous broadleaf trees               4    C3
!       8         tropical broadleaf trees                3    C3
!       9         drought deciduous trees                 5    C3
!       10        evergreen broadleaf shrub               3    C3
!       11        deciduous shrub                         2    C3
!       12        thorn shrubs                            2    C3
!       13        short grass and forbs                   8    C3
!       14        long grass                              9    C4
!       15        crops                                   6    C3
!       16        rice                                    6    C3
!       17        sugar                                   7    C4
!       18        maize                                   7    C4
!       19        cotton                                  6    C3
!       20        irrigated crops                         6    C3
!       21        urban
!       22        tundra                                  6    C3
!       23        swamp                                   6    C3
!       24        desert
!       25        mixed wood forests                      2    C3
!       26        mixed shrubs                            2    C3
!
!
!!!! NOL2PFTS (4)
     DATA  NOL2PFTS/2, 3, 2, 2/
!
!
!     CANOPY LIGHT/NITROGEN EXTINCTION COEFFICIENT - THIS BASICALLY
!     ASSUMES THAT MEAN PROFILE OF NITROGEN IS SAME AS THAT FOR
!     TIME MEAN PROFILE OF RADIATION - THE ASSUMPTION MADE BY SINGLE
!     BIG-LEAF MODELS
      DATA   KN/0.50, 0.50, 0.00,   &
                0.50, 0.50, 0.50,   &
                0.40, 0.48, 0.00,   &
                0.46, 0.44, 0.00/
!      
!     LOWER AND UPPER TEMPERATURE LIMITS FOR PHOTOSYNTHESIS, KELVIN
      DATA TLOW/268.1, 268.1, 0.000, &
                273.1, 273.1, 273.1, &
                270.1, 278.1, 0.000, &
                272.1, 283.1, 0.000/
!
!     UPPER LIMIT IN KELVIN
      DATA  TUP/320.0, 320.0, 0.000,  &
                320.0, 320.0, 325.0,  &
                325.0, 325.0, 0.000,  &
                325.0, 330.0, 0.000/
!     ARRAY TELLING WHICH VEGETATION TYPE IS C4
      DATA  ISC4/0, 0, 0,     &
                 0, 0, 0,     &
                 0, 1, 0,     &
                 0, 1, 0/
!
!     QUANTUM EFFICIENCIES, VALUES OF 0.08 & 0.04 ARE USED FOR C3 AND
!     C4 PLANTS, RESPECTIVELY
      DATA  ALPHA/0.08, 0.08, 0.00,  &
                  0.08, 0.08, 0.08,  &
                  0.08, 0.04, 0.00,  &
                  0.08, 0.04, 0.00/
!!
!     LEAF SCATERRING COEFFICIENTS, VALUES OF 0.15 & 0.17 ARE USED
!     FOR C3 AND C4 PLANTS, RESPECTIVELY
      DATA  OMEGA/0.15, 0.15, 0.00,   &
                  0.15, 0.15, 0.15,   &
                  0.15, 0.17, 0.00,   &
                  0.15, 0.17, 0.00/
!                  
!     PARAMETER M USED IN BWB PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
!     COUPLING. 
!
      DATA  MM/9.0, 9.0, 0.0,   &
              12.0,12.0,12.0,   &
              12.0, 6.0, 0.0,   &
              12.0, 6.0, 0.0/
!
!     PARAMETER B USED IN BWB PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
!     COUPLING.
      DATA  BB/0.01, 0.01, 0.00,   &
               0.01, 0.01, 0.01,   &
               0.01, 0.04, 0.00,   &
               0.01, 0.04, 0.00/
!
!     PARAMETER VPD0 USED IN LEUNING TYPE PHOTOSYNTHESIS - STOMATAL
!     CONDUCTANCE COUPLING, IN PASCALS
      DATA VPD0/2000., 2000., 0.000,   &
                2000., 2000., 2000.,   &
                1500., 1500., 0.000,   &
                1500., 1500., 0.000/
!
! !     EXPONENT FOR SOIL MOISTURE STRESS. FOR SN EQUAL TO 1, PHOTOSYNTHESIS
! !     DECREASES LINEARLY WITH SOIL MOISTURE, AND OF COURSE NON-LINEARLY
! !     FOR VALUES HIGHER THAN 1. WHEN SN IS ABOUT 10, PHOTOSYNTHESIS DOES
! !     NOT START DECREASING UNTIL ABOUT SOIL MOISTURE IS HALF WAY BETWEEN
! !     WILTING POINT AND FIELD CAPACITY.
! !
! !     PFT 3 HAS A HIGHER VALUE TO EMULATE DEEP ROOTS GIVING ACCESS TO
! !     GROUND WATER. JM 26.06.2012
!       DATA SN/2.0, 2.0, 0.0,   &
!               20.0, 2.0, 2.0,  &
!               2.0, 2.0, 0.0,   &
!               2.0, 2.0, 0.0/
! !     ADDITIONAL CONSTRAIN OF SOIL MOISTURE STRESS ON PHOTOSYNTHESIS.
! !     THIS CAN BE USED TO SIMULATE THE EFFECT OF IRRIGATION FOR CROPS.
! !
!       DATA SMSCALE/0.0, 0.0, 0.0,   &
!                    0.0, 0.0, 0.0,   &
!                    0.1, 0.1, 0.0,   &
!                    0.0, 0.0, 0.0/
! !
!     MAX. PHOTOSYNTHETIC RATE, MOL CO2 M^-2 S^-1
!     VALUES ARE MAINLY DERIVED FROM KATTHE ET AL. 2009 WHICH 
!     DOESN'T INCLUDE C4
      DATA VMAX/35.0E-06,  40.0E-06, 0.00E-06,   &
                51.0E-06,  67.0E-06, 40.0E-06,   &
                100.0E-06, 40.0E-06, 0.00E-06,   &
                75.0E-06, 15.0E-06, 0.00E-06/
!     NO. OF ITERATIONS FOR CALCULATING INTERCELLULAR CO2 CONCENTRATION
      DATA  REQITER/4/
!
!     MAX. INTERCELLULAR CO2 CONCENTRATION, PASCALS
      DATA CO2IMAX/2000.00/
!
!     PHOTOSYNTHESIS COUPLING OR CURVATURE COEFFICIENTS
      DATA BETA1/0.950/
      DATA BETA2/0.990/
!
!     PARAMETER TO INITIALIZE INTERCELLULAR CO2 CONC.
      DATA  INICO2I/0.65, 0.65, 0.00,     &
                    0.65, 0.65, 0.65,     &
                    0.65, 0.37, 0.00,     &
                    0.65, 0.37, 0.00/
!
!     LEAF MAINTENANCE RESPIRATION COEFFICIENTS
      DATA  RMLCOEFF/0.015, 0.017, 0.000,  &
                     0.020, 0.015, 0.015,  &
                     0.015, 0.025, 0.000,  &
                     0.013, 0.025, 0.000/
!
!     FREEZING TEMPERATURE
      DATA TFREZ/273.16/
!
!     STANDARD ATMOS. PRESSURE
      DATA STD_PRESS/101325.0/
!
!     ZERO
      DATA ZERO/1E-20/
!
!     PHOTOSYNTHESIS DOWN REGULATION PARAMETERS
!     EQUIVALENT CO2 FERTILIZATION EFFECT THAT WE WANT MODEL TO YIELD
      DATA GAMMA_W/0.45/
!
!     EQUIVALENT CO2 FERTILIZATION EFFECT THAT MODEL ACTUALLY GIVES
!     WITHOUT ANY PHOTOSYNTHESIS DOWN-REGULATION
      DATA GAMMA_M/0.95/
! 
!     --------------------------------------------------------------
!     DECIDE IF WE WANT TO USE BWB (1) OR LEUNING TYPE (2) PHOTOSYNTHESIS
!     STOMATAL CONDUCTANCE COUPLING
!
      PS_COUP=2
!
!     DECIDE IF WE WANT TO ESTIMATE PHOTOSYNTHETIC RATE AS A SMOOTHED
!     AVERAGE OF THE 3 LIMITING RATES, OR AS A MINIMUM OF THE 3 LIMITING
!     RATES, OR AS A MINIMUM OF THE RUBSICO AND LIGHT LIMITING RATES
!
      SMOOTH=.TRUE.
      MIN2=.FALSE.
      MIN3=.FALSE.
!!!!
!       NOTE: SINCE WE DON'T HAVE INPUT, CO2CONC IS SET TO CONSTANT C02CC 375.0
!       EXPERIMENTS INDICATES THE IMPACT OF THE VALUES OF CO2CONC IS SMALL
      CO2CC = 375.0
!!!!
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      DO I=1, n
         ! calculate aerodynamic conductance 
         CFLUX(I)=1./RESAVG(I)
      ENDDO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         INITIALIZATION
!          INITIALIZE REQUIRED ARRAYS TO ZERO FOR SINGLE LEAF MODEL

      CO2A           = 0.0
      IPAR           = 0.0
      SIGMA          = 0.0
      TGAMMA         = 0.0
      KC             = 0.0
      KO             = 0.0
      GB             = 0.0
      RC             = 0.0
      FC_TEST        = 0.0     
      FPAR           = 0.0
      VMAXC          = 0.0
      VMUNS1         = 0.0
      VMUNS2         = 0.0
      VMUNS3         = 0.0
      VMUNS          = 0.0
      VM             = 0.0
      JC1            = 0.0
      JC2            = 0.0
      JC3            = 0.0
      JC             = 0.0
      JE             = 0.0
      JE1            = 0.0
      JE2            = 0.0
      JS             = 0.0
      A_VEG          = 0.0
      RML_VEG        = 0.0
      AN_VEG         = 0.0
      CO2LS          = 0.0
      GC             = 0.0
      GCTU           = 0.0
      RC_VEG         = 5000.0
      USEBB          = 0.0
      !SM_FUNC        = 0.0   
!
!     GENERATE THE KK_TO_ICC INDEX FOR CORRESPONDENCE BETWEEN 9 PFTs AND THE
!     12 VALUES IN THE PARAMETER VECTORS
!
      ICOUNT=0
      DO J = 1, IC
         DO M = 1, NOL2PFTS(J)
            NN = (J-1)*L2MAX + M
            ICOUNT = ICOUNT + 1
            KK_TO_ICC(ICOUNT)=NN
         ENDDO
      ENDDO
!
!     INITIALIZATION ENDS
!
!     IF LAI IS LESS THAN SLAI THAN WE USE STORAGE LAI TO PHOTOSYNTHESIZE.
!     HOWEVER, WE DO NOT USE THE STOMATAL RESISTANCE ESTIMATED IN THIS CASE,
!     BECAUSE STORAGE LAI IS AN IMAGINARY LAI, AND WE SET STOMATAL RESISTANCE
!     TO ITS MAX. NOTE THAT THE CONCEPT OF
!     STORAGE/IMAGINARY LAI IS USED FOR PHENOLOGY PURPOSES AND THIS
!     IMAGINARY LAI ACTS AS MODEL SEEDS.
!
      DO I=1,N
         ! FCANC(I,1) = VEGFRAC(I,4)
         ! FCANC(I,2) = VEGFRAC(I,6)+VEGFRAC(I,11)+VEGFRAC(I,12)+VEGFRAC(I,25)  &
         !      +VEGFRAC(I,26)
         ! FCANC(I,3) = VEGFRAC(I,5)+VEGFRAC(I,8)+VEGFRAC(I,10)
         ! FCANC(I,4) = VEGFRAC(I,7)
         ! FCANC(I,5) = VEGFRAC(I,9)
         ! FCANC(I,6) = VEGFRAC(I,15)+VEGFRAC(I,16)+VEGFRAC(I,19)+VEGFRAC(I,20)  &
         !      +VEGFRAC(I,22)+VEGFRAC(I,23)
         ! FCANC(I,7) = VEGFRAC(I,17)+VEGFRAC(I,18)
         ! FCANC(I,8) = VEGFRAC(I,13)
         ! FCANC(I,9) = VEGFRAC(I,14)
!
         ! AILCG(I,1) = LAI_NCLASS(I,4) 
         ! AILCG(I,2) = LAI_NCLASS(I,6)+LAI_NCLASS(I,11)+LAI_NCLASS(I,12) &
         !      +LAI_NCLASS(I,25)+LAI_NCLASS(I,26)
         ! AILCG(I,3) = LAI_NCLASS(I,5)+LAI_NCLASS(I,8)+LAI_NCLASS(I,10)
         ! AILCG(I,4) = LAI_NCLASS(I,7)
         ! AILCG(I,5) = LAI_NCLASS(I,9)
         ! AILCG(I,6) = LAI_NCLASS(I,15)+LAI_NCLASS(I,16)+LAI_NCLASS(I,19) &
         !      +LAI_NCLASS(I,20)+LAI_NCLASS(I,22)+LAI_NCLASS(I,23) 
         ! AILCG(I,7) = LAI_NCLASS(I,17)+LAI_NCLASS(I,18)
         ! AILCG(I,8) = LAI_NCLASS(I,13)
         ! AILCG(I,9) = LAI_NCLASS(I,14)

         ! class 10 evergeen broadleaf shrub... mostly arid regions i think, so in class 5 even if not evergreen, same for class 12:thorn shrubs. 

         ! long grass in c3 grass
         ! urban in 0.6 grass 0.2 broadleef cold, 0.2 needle leaf evergeen 
         ! 25 split between class 4 broadleaf (60%)  1 needleleaf evergreen (40%)
         ! 26 split between class 4 broadleaf (80%) and needleleaf evergreen (20%) below 50 lat
         ! 26 split between needleleaf evergreen (50%) grass (20%) broadleaf cold (30%)

         FCANC(I,1) =                 VEGFRAC(I, 4) + &
                                0.2 * VEGFRAC(I,21) + &
                                0.4 * VEGFRAC(I,25) + &
             REAL(MASKLAT(I)) * 0.2 * VEGFRAC(I,26) + &
        REAL( 1 - MASKLAT(I)) * 0.5 * VEGFRAC(I,26) 

         FCANC(I,2) =                 VEGFRAC(I, 6)
        

         FCANC(I,3) =                 VEGFRAC(I, 5) + &
                                      VEGFRAC(I, 8)


         FCANC(I,4) =                 VEGFRAC(I, 7) + &
                                      VEGFRAC(I,11) + &
                                0.2 * VEGFRAC(I,21) + &
                                0.6 * VEGFRAC(I,25) + &
             REAL(MASKLAT(I)) * 0.8 * VEGFRAC(I,26) + &
        REAL( 1 - MASKLAT(I)) * 0.3 * VEGFRAC(I,26) 


         FCANC(I,5) =                 VEGFRAC(I, 9) + &
                                      VEGFRAC(I,10) + &
                                      VEGFRAC(I,12)


         FCANC(I,6) =                 VEGFRAC(I,15) + &
                                      VEGFRAC(I,16) + &
                                      VEGFRAC(I,19) + &
                                      VEGFRAC(I,20) + &
                                      VEGFRAC(I,22) + &
                                      VEGFRAC(I,23)

         FCANC(I,7) =                 VEGFRAC(I,17) + &
                                      VEGFRAC(I,18)

         FCANC(I,8) =                 VEGFRAC(I,13) + &
                                      VEGFRAC(I,14) + &
                                0.6 * VEGFRAC(I,21) + &
        REAL( 1 - MASKLAT(I)) * 0.2 * VEGFRAC(I,26) 

         FCANC(I,9) = 0.0

!    FOR LAI, LAI_NCLASS has already seen vegfrac...
! NEED TO DIVIDE BY VEGETATION FRACTION



         AILCG(I,1) =                 LAI_NCLASS(I, 4) + &
                                0.2 * LAI_NCLASS(I,21) + &
                                0.4 * LAI_NCLASS(I,25) + &
             REAL(MASKLAT(I)) * 0.2 * LAI_NCLASS(I,26) + &
        REAL( 1 - MASKLAT(I)) * 0.5 * LAI_NCLASS(I,26) 

         AILCG(I,2) =                 LAI_NCLASS(I, 6)
        

         AILCG(I,3) =                 LAI_NCLASS(I, 5) + &
                                      LAI_NCLASS(I, 8)


         AILCG(I,4) =                 LAI_NCLASS(I, 7) + &
                                      LAI_NCLASS(I,11) + &
                                0.2 * LAI_NCLASS(I,21) + &
                                0.6 * LAI_NCLASS(I,25) + &
             REAL(MASKLAT(I)) * 0.8 * LAI_NCLASS(I,26) + &
        REAL( 1 - MASKLAT(I)) * 0.3 * LAI_NCLASS(I,26) 


         AILCG(I,5) =                 LAI_NCLASS(I, 9) + &
                                      LAI_NCLASS(I,10) + &
                                      LAI_NCLASS(I,12)


         AILCG(I,6) =                 LAI_NCLASS(I,15) + &
                                      LAI_NCLASS(I,16) + &
                                      LAI_NCLASS(I,19) + &
                                      LAI_NCLASS(I,20) + &
                                      LAI_NCLASS(I,22) + &
                                      LAI_NCLASS(I,23)

         AILCG(I,7) =                 LAI_NCLASS(I,17) + &
                                      LAI_NCLASS(I,18)

         AILCG(I,8) =                 LAI_NCLASS(I,13) + &
                                      LAI_NCLASS(I,14) + &
                                0.6 * LAI_NCLASS(I,21) + &
        REAL( 1 - MASKLAT(I)) * 0.2 * LAI_NCLASS(I,26) 

         AILCG(I,9) = 0.0


         DO M=1,9

            if (FCANC(i,m).gt.0) then
               !print * , 'for i=',i,'m=',m
               !print *, 'lai=',ailcg(i,m),'fcanc=',fcanc(i,m)
               ailcg(i,m) = max( ailcg(i,m) / fcanc(i,m), 1.E-6)
            else
               ailcg(i,m) = 1.E-6
            endif
         ENDDO

         
         QSWV(I)=max(0.0, 0.5*QSWV1(I))

      ENDDO
      DO I=1,N
         FC(I)=FCANC(I,1)+FCANC(I,2)+FCANC(I,3)+FCANC(I,4)+FCANC(I,5)    &
              +FCANC(I,6)+FCANC(I,7)+FCANC(I,8)+FCANC(I,9)
      ENDDO
!     SET MIN. AND MAX. VALUES FOR STOMATAL CONDUCTANCE. WE MAKE SURE
!     THAT MAX. STOMATAL RESISTANCE IS AROUND 5000 S/M AND MIN. STOMATAL
!     RESISTANCE IS 51 S/M. 
!
       DO J = 1, ICC 
         DO I=1, N
            GCMIN(I,J)=0.0002 * (TFREZ/TCAN(I)) * (1./0.0224) *     &
                 (PRESSG(I)/STD_PRESS)
!
            GCMAX(I,J)=0.1    * (TFREZ/TCAN(I)) * (1./0.0224) *   &
                 (PRESSG(I)/STD_PRESS)        
!
         ENDDO
       ENDDO
!     IF WE ARE USING LEUNING TYPE PHOTOSYNTHESIS-STOMATAL CONDUCTANCE
!     COUPLING WE NEED VAPOR PRESSURE DEFICIT AS WELL. CALCULATE THIS
!     FROM THE RH AND AIR TEMPERATURE WE HAVE. WE FIND E_SAT, E, AND VPD
!     IN PASCALS.
!
       IF(PS_COUP.EQ.1)THEN

          DO I = 1, N
             IF(TCAN(I).GE.TFREZ) THEN
                CA=17.269
                CB=35.86
             ELSE
                CA=21.874
                CB=7.66
             ENDIF
             EA     = QA(I)*PRESSG(I)/(0.622+0.378*QA(I))
             EASAT  = 611.0*EXP(CA*(TCAN(I)-TFREZ)/(TCAN(I)-CB))
             RH(I)  = EA/EASAT
             
          ENDDO

       ELSE IF(PS_COUP.EQ.2)THEN  
          
          VPD_TERM  =0.0
          
          DO I = 1, N
             VPD(I)=0.0
             IF(TCAN(I).GE.TFREZ) THEN
                CA=17.269
                CB=35.86
             ELSE
                CA=21.874
                CB=7.66
             ENDIF
             EA     = QA(I)*PRESSG(I)/(0.622+0.378*QA(I))
             EASAT  = 611.0*EXP(CA*(TCAN(I)-TFREZ)/(TCAN(I)-CB))
             RH(I)  = EA/EASAT
             VPD(I) = EASAT-EA
             VPD(I)=MAX(0.0,VPD(I))
         ENDDO
!
         K1 = 0
         DO J = 1, IC
            IF(J.EQ.1) THEN
               K1 = K1 + 1
            ELSE
               K1 = K1 + NOL2PFTS(J-1)
            ENDIF
            K2 = K1 + NOL2PFTS(J) - 1
            DO M = K1, K2
               DO I = 1, N
                  VPD_TERM(I,M)=1.0/( 1.0 +( VPD(I)/VPD0(KK_TO_ICC(M)) ) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!     ESTIMATE PARTIAL PRESSURE OF CO2 AND IPAR
!
      DO I = 1, N
         IF(FC(I).GT.0.)THEN
!       ----CONVERT CO2CC FROM PPM TO PASCALS
            CO2A(I)=CO2CC * PRESSG(I)/1E+06
!
!       ---CHANGE PAR FROM W/M^2 TO MOL/M^2.S
!
            IPAR(I) = QSWV(I)*4.6*1E-06
!
!       SUNLIT PART GETS BOTH DIRECT AND DIFFUSED, WHILE
!       THE SHADED PART GETS ONLY DIFFUSED
!
         ENDIF
      ENDDO
!!!!
      K1 = 0
      DO J = 1, IC
         IF(J.EQ.1) THEN
            K1 = K1 + 1
         ELSE
            K1 = K1 + NOL2PFTS(J-1)
         ENDIF
         K2 = K1 + NOL2PFTS(J) - 1
         DO M = K1, K2
            DO I = 1, N
               IF(FCANC(I,M).GT.ZERO)THEN
!
!         FIND FPAR - FACTOR FOR SCALING PHOTOSYNTHESIS TO CANOPY
!         BASED ON ASSUMPTION THAT NITROGEN IS OPTIMALLY DISTRIBUTED.
!         IN THE END ADD CONDUCTANCE AND NET PHOTOSYNTHESIS
!         TO GET THE TOTAL.
!
                  FPAR(I,M)=(1.0/KN(KK_TO_ICC(M)))*(1.0-EXP(-KN(KK_TO_ICC(M))            &
                       *AILCG(I,M)))
!
!         FIND Vmax,canopy, THAT IS Vmax SCALED BY LAI FOR THE SINGLE
!         LEAF MODEL
!
                  VMAXC(I,M)=VMAX(KK_TO_ICC(M)) * FPAR(I,M)
!
!         FIND Vm,unstressed (DUE TO WATER) BUT STRESSED DUE TO TEMPERATURE
!
                  Q10 = 2.00
                  Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
!
                  IF(COSZS(I).GT.0.0)THEN           
                     VMUNS1(I,M) = VMAXC(I,M) * Q10_FUNC
                  ENDIF
!
!         ASSUMING THAT SUNLIT AND SHADED TEMEPERATURES ARE SAME
!
                  VMUNS2(I,M) = (1. + EXP(0.3*(TCAN(I)-TUP(KK_TO_ICC(M)))))
                  VMUNS3(I,M) = (1. + EXP(0.3*(TLOW(KK_TO_ICC(M))-TCAN(I))))
!
                  VMUNS(I,M)=VMUNS1(I,M) / (VMUNS2(I,M)*VMUNS3(I,M))         
!
               ENDIF
            ENDDO
         ENDDO
      ENDDO


      




!!!!
!!!
!
!     CALCULATE SOIL MOIS STRESS TO ACCOUNT FOR REDUCTION IN PHOTOSYN
!     DUE TO LOW SOIL MOISTURE
!
!     FOR EACH SOIL LAYER LOOK AT SOIL MOISTURE  VS. FIELD CAPACITY AND WILTING POINT
!     NOTE THAT WE USE AVERAGE ROOT PROFILE CALCULATE IN VEGI SO THAT THE RESISTANCE IS
!     CONSISTANT WITH THE ROOT PROFILE USED IN HYDRO
!     IN ORIGINAL CODE, ONE ROOT PROFILE FOR EACH PFT CLASS. HERE ALL CLASSES SHARE THE SAME
!     ROOT PROFILE

      DO I=1,N

         DO K=1,NL_SVS
            ! beta between 0. an 1.  , 
            BETA_WSOL(I,K) =  min( max( wd(i,k) - wwilt(i,k) , 0.0) / (wfc(i,k) - wwilt(i,k)) , 1.0)
            ! soil moisture stress term per layer
            G_WSOL(i,k) = 1.0 - ( 1.0 - BETA_WSOL(I,K) ) ** GEXP
         ENDDO
         ! average soil moisture term ... weighted by root fractions 
         ! Total roots = 1.0 if vegetation present ... set the term=0.0 if no vegetation
         if (FCD(I,NL_SVS).gt.0.999) then
             avg_gwsol(i) = g_wsol(i,1) * fcd(i,1)
             do k=2,NL_SVS
                avg_gwsol(i) = avg_gwsol(i)  + g_wsol(i,k) * ( fcd(i,k) - fcd(i,k-1) ) 
             enddo
         else
            ! no vegetation
            avg_gwsol(i) = 0.0
         endif
         
      ENDDO






!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!
!     USE SOIL MOISTURE FUNCTION TO MAKE Vm,unstressed -> Vm STRESSED
!
      DO J = 1, ICC
         DO I = 1, N
            IF((FCANC(I,J).GT.ZERO).AND.(COSZS(I).GT.0.0)) THEN
               VM(I,J) = VMUNS(I,J) * AVG_GWSOL(I)
            ENDIF
         ENDDO
      ENDDO
!
!     FIND TEMPERATURE DEPENDENT PARAMETER VALUES
!
      DO I = 1, N
!
!       FIND RUBISCO SPECIFICITY FOR CO2 RELATIVE TO O2 - SIGMA
!
         Q10 = 0.57
         Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
         SIGMA(I) = 2600.0 * Q10_FUNC
!
!       FIND CO2 COMPENSATION POINT USING RUBISCO SPECIFICITY - TGAMMA.
!       KEEP IN MIND THAT CO2 COMPENSATION POINT FOR C4 PLANTS IS ZERO,
!       SO THE FOLLOWING VALUE IS RELEVANT FOR C3 PLANTS ONLY
!
         O2_CONC(I) = .2095 * PRESSG(I)
         TGAMMA(I) = O2_CONC(I) / (2.0*SIGMA(I))
!
!       ESTIMATE MICHELIS-MENTON CONSTANTS FOR CO2 (Kc) and O2 (Ko) TO
!       BE USED LATER FOR ESTIMATING RUBISCO LIMITED PHOTOSYNTHETIC RATE
!
         Q10 = 2.10
         Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
         KC(I) = 30.0 * Q10_FUNC
!
         Q10 = 1.20
         Q10_FUNC = Q10**(0.1*(TCAN(I)-298.16))
         KO(I) = 30000.0 * Q10_FUNC
!
      ENDDO
!
!     CHOOSE A VALUE OF INTERCELLULAR CO2 CONCENTRATION (CO2i) 
!     ITERATIVE PROCESS, NUMBER OF ITERATION SET ABOVE 
!
      DO IT_COUNT = 1, REQITER
!
         DO J = 1, ICC
            DO I = 1,N               
               CO2I(I,J) = CO2I1(I,J)
               IF(CO2I(I,J).LE.ZERO)  CO2I(I,J) = INICO2I(KK_TO_ICC(J)) * CO2A(I)
            ENDDO
         ENDDO
!
         !     ESTIMATE RUBISCO LIMITED PHOTOSYNTHETIC RATE
!
         DO J = 1, ICC
            DO I = 1,N
               IF((FCANC(I,J).GT.ZERO).AND.(COSZS(I).GT.0.0)) THEN
                     
                  JC1(I,J) = CO2I(I,J) - TGAMMA(I)
                  JC2(I,J) = KC(I) * (1.0 + (O2_CONC(I)/KO(I)) )
                  JC3(I,J) = CO2I(I,J) + JC2(I,J)
                        !
                  IF (ISC4(KK_TO_ICC(J)).EQ.1) THEN
                     JC(I,J)  = VM(I,J)
                  ELSE IF (ISC4(KK_TO_ICC(J)).EQ.0) THEN
                     JC(I,J)  = VM(I,J) * (JC1(I,J)/JC3(I,J))
                  ENDIF
                     
               ENDIF
!               
            ENDDO
         ENDDO
!!!!!!!!!!!!!!!!
!
!     ESTIMATE PHOTOSYNTHETIC RATE LIMITED BY AVAILABLE LIGHT
!
         DO  J = 1, ICC
            DO  I = 1, N
               IF((FCANC(I,J).GT.ZERO).AND.(COSZS(I).GT.0.0)) THEN
                  
                  JE1(I,J)= FPAR(I,J)*ALPHA(KK_TO_ICC(J))*(1.0-OMEGA(KK_TO_ICC(J)))
                  JE2(I,J)=( CO2I(I,J)-TGAMMA(I) ) /              &
                       ( CO2I(I,J)+ (2.0*TGAMMA(I)) )
!
                  IF (ISC4(KK_TO_ICC(J)).EQ.1) THEN
                     JE(I,J) = IPAR(I) * JE1(I,J)
                  ELSE IF (ISC4(KK_TO_ICC(J)).EQ.0) THEN
                     JE(I,J) = IPAR(I) * JE1(I,J) * JE2(I,J)
                  ENDIF
                  
               ENDIF
!
            ENDDO
         ENDDO
!!!!
!
!     ESTIMATE PHOTOSYNTHETIC RATE LIMITED BY TRANSPORT CAPACITY
!
         DO J = 1, ICC
            DO I = 1,N
               IF((FCANC(I,J).GT.ZERO).AND.(COSZS(I).GT.0.0)) THEN
                   
                  IF (ISC4(KK_TO_ICC(J)).EQ.1) THEN
                     JS(I,J) = 20000.0 * VM(I,J)*(CO2I(I,J)/PRESSG(I))
                  ELSE IF (ISC4(KK_TO_ICC(J)).EQ.0) THEN
                     JS(I,J) = 0.5 * VM(I,J)
                  ENDIF
                   
               ENDIF
            ENDDO
         ENDDO
!
!!!!!!!!!!!!!!!!
!
!     INCLUDE NUTRIENT LIMITATION EFFECT BY DOWN-REGULATING PHOTOSYNTHESIS
!     N_EFFECT DECREASES FROM 1.0 AS CO2 INCREASES ABOVE 288 PPM.
         !
         DELTA_CO2  = CO2CC - 288.0
         TEMP_R = LOG(1.0 + (DELTA_CO2 /288.0))
         TEMP_B = 1.0 + TEMP_R * (GAMMA_W)
         TEMP_C = 1.0 + TEMP_R * (GAMMA_M)
         !       LIMIT N_EFFECT TO MAX OF 1.0 SO THAT NO UP-REGULATION OCCURS
         N_EFFECT  = MIN(  TEMP_B / TEMP_C  , 1.0 )
!
!     FIND THE SMOOTHED AVERAGE OF THREE PHOTOSYNTHETIC RATES JC, JE,
!     AND JS USING COLLATZ'S TWO QUADRATIC EQUATIONS, OR FIND THE MIN.
!     OF THIS TWO RATES OR FIND MIN. OF JC AND JE.
!
      DO J = 1, ICC
         DO I = 1, N
            IF((FCANC(I,J).GT.ZERO).AND.(COSZS(I).GT.0.0)) THEN
                  
               IF(SMOOTH)THEN
                  !
                  TEMP_B  = JC(I,J) + JE(I,J)
                  TEMP_C  = JC(I,J) * JE(I,J)
                  TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA1 * TEMP_C), 0.0)
                  TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA1)
                  TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA1)
                  TEMP_JP = MIN(TEMP_Q1, TEMP_Q2)
                  !
                  TEMP_B  = TEMP_JP + JS(I,J)
                  TEMP_C  = TEMP_JP * JS(I,J)
                  TEMP_R  = MAX( (TEMP_B**2 - 4. * BETA2 * TEMP_C), 0.0)
                  TEMP_Q1 = (TEMP_B + SQRT(TEMP_R)) / (2. * BETA2)
                  TEMP_Q2 = (TEMP_B - SQRT(TEMP_R)) / (2. * BETA2)
                  A_VEG(I,J) = MIN(TEMP_Q1, TEMP_Q2)
               ELSEIF(MIN2) THEN
                  A_VEG(I,J)=MIN(JC(I,J), JE(I,J))
               ELSEIF(MIN3) THEN
                  A_VEG(I,J)=MIN(JC(I,J), JE(I,J), JS(I,J))
               ENDIF
!           DOWN-REGULATE PHOTOSYNTHESIS FOR C3 PLANTS
               IF(ISC4(KK_TO_ICC(J)).EQ.0) A_VEG(I,J) = A_VEG(I,J) * N_EFFECT
         
               A_VEG(I,J) = MAX(0.0, A_VEG(I,J))
            ENDIF
               !
         ENDDO
      ENDDO
!
!     ESTIMATE LEAF MAINTENANCE RESPIRATION RATES AND NET PHOTOSYNTHETIC
!     RATE. THIS NET PHOSYNTHETIC RATE IS /M^2 OF VEGETATED LAND.
!
      DO J = 1, ICC
         DO I = 1, N
            IF(FCANC(I,J).GT.ZERO)THEN
               !
!         RECENT STUDIES SHOW RmL IS LESS TEMPERATURE SENSITIVE THAN
!         PHOTOSYNTHESIS DURING DAY, THAT'S WHY A SMALL Q10 VALUE IS
!         USED DURING DAY.
!
               Q10_FUNCN = 2.00**(0.1*(TCAN(I)-298.16))
               Q10_FUNCD = 1.30**(0.1*(TCAN(I)-298.16))
               !         
               IF(COSZS(I).GT.0.0)THEN
                  RML_VEG(I,J) = RMLCOEFF(KK_TO_ICC(J))*VMAXC(I,J)*Q10_FUNCD
               ELSE
                  RML_VEG(I,J) = RMLCOEFF(KK_TO_ICC(J))*VMAXC(I,J)*Q10_FUNCN 
               ENDIF
               AN_VEG(I,J) = A_VEG(I,J) - RML_VEG(I,J)
               !
            ENDIF
         ENDDO
      ENDDO
!
!     FIND CO2 CONCENTRATION AT LEAF SURFACE FOR ALL VEGETATION TYPES.
!     ALTHOUGH WE ARE FINDING CO2 CONC AT THE LEAF SURFACE SEPARATELY
!     FOR ALL VEGETATION TYPES, THE BIG ASSUMPTION HERE IS THAT THE
!     AERODYNAMIC CONDUCTANCE IS SAME OVER ALL VEGETATION TYPES. CLASS
!     FINDS AERODYNAMIC RESISTANCE OVER ALL THE 4 SUB-AREAS, BUT NOT
!     FOR DIFFERENT VEGETATION TYPES WITHIN A SUB-AREA.
!     ALSO CHANGE AERODYNAMIC CONDUCTANCE, CFLUX, FROM M/S TO MOL/M^2/S
!
      DO J = 1, ICC
         DO I = 1, N
            IF(FCANC(I,J).GT.ZERO)THEN
               !
               GB(I)=CFLUX(I) * (TFREZ/TCAN(I))*(PRESSG(I)/STD_PRESS) &
                    *(1./0.0224)
               GB(I)= MIN(10.0, MAX(0.1, GB(I)) )  
               ! 
               CO2LS(I,J) = 0.5*(CO2LS(I,J)+                              &
                    (CO2A(I)-( (AN_VEG(I,J)*1.37*PRESSG(I)) / GB(I)))  )
               CO2LS(I,J) = MAX (1.05*TGAMMA(I) , CO2LS(I,J))             
!
            ENDIF
         ENDDO
      ENDDO
      !
!     FIND STOMATAL CONDUCTANCE AS PER BALL-WOODROW-BERRY FORMULATION
!     USED BY COLLATZ ET AL. OR USE THE LEUNING TYPE FORMULATION WHICH
!     USES VPD INSTEAD OF RH
!
      DO J = 1, ICC
         DO I = 1, N
            IF(FCANC(I,J).GT.ZERO)THEN
               !
!         IF LIGHT IS TOO LESS MAKE PARAMETER BB VERY SMALL
               IF(QSWV(I).LT.2.0) THEN
                  USEBB(J)=0.001
               ELSE
                  USEBB(J)=BB(KK_TO_ICC(J))
               ENDIF
               RH(I)=MAX(0.3, RH(I))          !FOLLOWING IBIS
               !
               
               IF(PS_COUP.EQ.1)THEN
                  GC(I,J)=( (MM(KK_TO_ICC(J))*RH(I)*PRESSG(I)*AN_VEG(I,J) )/        &
                       CO2LS(I,J) ) + USEBB(J)*AILCG(I,J)*AVG_GWSOL(I)

               ELSE IF(PS_COUP.EQ.2) THEN
                  GC(I,J)=( (MM(KK_TO_ICC(J))*VPD_TERM(I,J)*PRESSG(I)*AN_VEG(I,J) )/ &
                       (CO2LS(I,J)-TGAMMA(I)) ) +                             &
                       USEBB(J)*AILCG(I,J)*AVG_GWSOL(I)

               ENDIF
!
               GC(I,J)=MAX(GCMIN(I,J),                                    &
                    USEBB(J)*AILCG(I,J)*AVG_GWSOL(I), GC(I,J))

               GC(I,J)=MIN(GCMAX(I,J), GC(I,J))   
!
            ENDIF
         ENDDO
      ENDDO
      !     FIND THE INTERCELLULAR CO2 CONCENTRATION BASED ON ESTIMATED
!     VALUE OF GC
!
      DO J = 1, ICC
         DO I = 1, N
            IF(FCANC(I,J).GT.ZERO)THEN
!     
               CO2I(I,J)= 0.5*(CO2I(I,J) + (CO2LS(I,J) -               & 
                    ( (AN_VEG(I,J)*1.65*PRESSG(I))/GC(I,J) ) ) )
               CO2I(I,J)=MAX(1.05*TGAMMA(I), MIN(CO2IMAX, CO2I(I,J)))
               PREV_CO2I(I,J)=CO2I(I,J)          
!
            ENDIF
         ENDDO
      ENDDO
!
      DO J = 1, ICC
         DO I = 1, N
            IF(FCANC(I,J).GT.ZERO)THEN
               CO2I1(I,J) = PREV_CO2I(I,J)
            ENDIF
         ENDDO
      ENDDO
!
   ENDDO ! ITERATIVE LOOP
!
!     WHEN REQUIRED NO. OF ITERATIONS HAVE BEEN PERFORMED THEN FIND
!     STOMATAL CONDUCTANCES FOR ALL VEGETATION TYPES IN M/S AND THEN
!     USE CONDUCTANCES TO FIND RESISTANCES. GCTU IMPLIES GC IN TRADITIONAL
!     UNITS OF M/S
!
   DO J = 1, ICC
      DO I = 1, N
         IF(FCANC(I,J).GT.ZERO) THEN
!         
            GCTU(I,J)= GC(I,J) * (TCAN(I)/TFREZ) * &
                      (STD_PRESS/PRESSG(I))*0.0224
            RC_VEG(I,J) = 1./GCTU(I,J)
           
         ENDIF
      ENDDO
   ENDDO
!
!
!     AND FINALLY TAKE WEIGHTED AVERAGE OF RC_VEG BASED ON FRACTIONAL
!     COVERAGE OF OUR 4 VEGETATION TYPES
!
   DO I = 1, N
      DO J = 1, ICC
         RC(I)=RC(I)+FCANC(I,J)*RC_VEG(I,J)
      ENDDO
!
      FC_TEST(I)=FCANC(I,1)+FCANC(I,2)+FCANC(I,3)+FCANC(I,4)+  &
                 FCANC(I,5)+FCANC(I,6)+FCANC(I,7)+FCANC(I,8)+  &
                 FCANC(I,9)
!
      IF(FC_TEST(I).GT.ZERO)THEN
         RC(I)=RC(I)/FC_TEST(I)
      ELSE
         RC(I) = 5000.0     
      ENDIF

   ENDDO
!
!
!     CONVERT AN_VEG AND RML_VEG TO u-MOL CO2/M2.SEC
!
      DO J = 1, ICC
        DO I = 1, N
          AN_VEG(I,J)=AN_VEG(I,J)*1.0E+06
          RML_VEG(I,J)=RML_VEG(I,J)*1.0E+06
        ENDDO
      ENDDO
!

      RETURN
    END SUBROUTINE PHTSYN_SVS
