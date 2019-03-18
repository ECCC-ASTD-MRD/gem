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
!**S/P MCICA_CLD_GENERATOR - PART OF THE ISCCP CLOUD SIMULATOR PACKAGE
!
      SUBROUTINE MCICA_CLD_GENERATOR(ZF, AVG_CF, AVG_QLW,AVG_QIW, &
                                     RLC_CF, RLC_CW,SIGMA_QCW, &
                                     ILG, IL1, IL2, LAY, &
                                     ICEWCIN_S,LIQWCIN_S, NCLDY)

      implicit none
#include <arch_specific.hf>
#include "mcica.cdk"
#include "tdpack_const.hf"
!#include "xcwdata.cdk"

! INPUT DATA
!

      REAL, INTENT(IN) :: &
       ZF(ILG,LAY),             &! Full-level (layer midpoint) altitude (km)
       AVG_CF(ILG,LAY),         &! Cloud fraction for each layer
       AVG_QLW(ILG,LAY),        &! Cloud mean liquid water mixing ratio (g/m^3)
       AVG_QIW(ILG,LAY),        &! Cloud mean ice mixing ratio          (g/m^3)
       SIGMA_QCW(ILG,LAY),      &! Normalized cloud condensate std. dev.
       RLC_CF(ILG,LAY),         &! Cloud fraction decorrelation length
       RLC_CW(ILG,LAY)         ! Cloud condensate decorrelation length

      INTEGER, INTENT(IN) :: &
       ILG, &
       IL1, &
       IL2, &
       LAY

!
! OUTPUT DATA
!

      REAL,INTENT(OUT) :: &
       ICEWCIN_S(ILG,LAY,NX_LOC),      &! Column ice water mixing ratio profile    (g/m^3)
       LIQWCIN_S(ILG,LAY,NX_LOC)      ! Column liquid water mixing ratio profile (g/m^3)

      INTEGER, INTENT(OUT) :: &
       NCLDY(ILG)                     ! Number of cloudy subcolumns

!Author
!        Jason Cole, MSC/Cloud Physics (Summer 2005)

!Revisions
! 001    ...

!Object
!        Generate a NX column by LAY layer cloud field. Input profiles
!        must be from top of the model to the bottom as is the output.
!        A Method and model evaluation is described in:

!        Raisanen et al, 2004, Stochastic generation of subgrid-scale
!        cloudy columns for large-scale models, Quart. J. Roy. Meteorol.
!        Soc., 2004 , 130 , 2047-2068


!Arguments
!                           - Input -
! ZF       (ILG,LAY)        Full-level (layer midpoint) altitude (km)
! AVG_CF   (ILG,LAY)        Cloud fraction for each layer
! AVG_QLW  (ILG,LAY)        Cloud mean liquid water mixing ratio (g/m^3)
! AVG_QIW  (ILG,LAY)        Cloud mean ice mixing ratio          (g/m^3)
! SIGMA_QCW(ILG,LAY)        Normalized cloud condensate std. dev.
! RLC_CF   (ILG,LAY)        Cloud fraction decorrelation length
! RLC_CW   (ILG,LAY)        Cloud condensate decorrelation length
! ILG, IL1, IL2             Horizontal dimension
! LAY                       Number of layers in GCM
! NX_LOC                    Number of columns to generate

!                           - Ouput -
! ICEWCIN_S(ILG,LAY,NX_LOC) Column ice water mixing ratio profile    (g/m^3)
! LIQWCIN_S(ILG,LAY,NX_LOC) Column liquid water mixing ratio profile (g/m^3)
! NCLDY    (ILG)            Number of cloudy subcolumns

!Implicites

      REAL, PARAMETER :: &
       CUT = 0.001

!
! LOCAL DATA
!

      INTEGER :: &
       IL,          &! Counter over ILG GCM columns
       I,           &! Counter
       K,           &! Counter over LAY vertical layers
       K_TOP,       &! Index of top most cloud layer
       K_BASE,      &! Index of lowest cloud layer
       IND1,        &! Index in variability calculation
       IND2        ! Index in variability calculation

      REAL :: &
       ALPHA(ILG,LAY),  &! Fraction of maximum/random cloud overlap
       RCORR(ILG,LAY),  &! Fraction of maximum/random cloud condensate overlap
       RIND1,           &! Real index in variability calculation
       RIND2,           &! Real index in variability calculation
       ZCW             ! Ratio of cloud condensate miximg ratio for this cell to its layer cloud-mean value

      REAL ::  &! Random number vectors
       X(ILG), &
       Y(ILG), &
       X1(ILG), &
       Y1(ILG), &
       X2(ILG), &
       Y2(ILG)

      INTEGER :: &
       I_LOC(ILG) ! Place the new subcolumns into the arrays starting from the front

      LOGICAL :: &
       L_TOP, &
       L_BOT, &
       L_CLD(ILG), &
       L_CLDSUB(ILG)

! Initialize the arrays

      DO IL = IL1, IL2
         I_LOC(IL) = 1
         L_CLD(IL) = .FALSE.
      END DO

! Find uppermost cloudy layer

      L_TOP = .FALSE.
      K_TOP = 0
      DO K=1,LAY
         DO IL = IL1, IL2
            K_TOP = K
            IF (AVG_CF(IL,K) .GT. CUT) L_TOP = .TRUE.
         END DO ! IL
         IF (L_TOP) EXIT
      END DO ! K

! If no cloudy layers in any GCM column in this group, exit

      IF (K_TOP .EQ. 0) THEN
         RETURN
      END IF

! Find lowermost cloudy layer

      L_BOT = .FALSE.

      K_BASE = 0
      DO K=LAY,1,-1
         DO IL = IL1, IL2
            K_BASE = K
            IF (AVG_CF(IL,K) .GT. CUT) L_BOT = .TRUE.
         END DO ! IL
         IF (L_BOT) EXIT
      END DO ! K

! Calculate overlap factors ALPHA for cloud fraction and RCORR for cloud
! condensate based on layer midpoint distances and decorrelation depths

      DO K=1, LAY !K_TOP,K_BASE-1
         DO IL = IL1, IL2
!            IF (RLC_CF(IL,K) .GT. 0.0) THEN
!              ALPHA(IL,K) = EXP(-(ZF(IL,K) - ZF(IL,K+1)) / RLC_CF(IL,K))
               ALPHA(IL,K) = -(ZF(IL,K) - ZF(IL,K+1)) / RLC_CF(IL,K)
!            ELSE
!               ALPHA(IL,K) = 0.0
!            END IF
!            IF (RLC_CW(IL,K) .GT. 0.0) THEN
!              RCORR(IL,K) = EXP(-(ZF(IL,K) - ZF(IL,K+1)) / RLC_CW(IL,K))
              RCORR(IL,K) = -(ZF(IL,K) - ZF(IL,K+1)) / RLC_CW(IL,K)
!            ELSE
!               RCORR(IL,K) = 0.0
!            END IF
         END DO ! IL
      END DO ! K

      CALL VSEXP(ALPHA,ALPHA,LAY*(IL2-IL1+1))
      CALL VSEXP(RCORR,RCORR,LAY*(IL2-IL1+1))

      DO I=1,NX_LOC

! Generate all subcolumns for latitude chain

         DO K = K_TOP, K_BASE

            IF (K .EQ. K_TOP) THEN
               CALL RANDOM_NUMBER(X)
               CALL RANDOM_NUMBER(Y)
            END IF

            CALL RANDOM_NUMBER(X1)
            CALL RANDOM_NUMBER(Y1)

            CALL RANDOM_NUMBER(X2)
            CALL RANDOM_NUMBER(Y2)

            DO IL = IL1, IL2

! Maximum-random overlap
               IF (LMAXRAN) THEN
                  IF (X(IL) .LT. 1.0-AVG_CF(IL,K-1)) THEN !It is clear above
                     X(IL) = X1(IL) * (1.0 - AVG_CF(IL,K-1))
                     Y(IL) = Y1(IL)
                  END IF
! Generalized overlap
               ELSE
                  IF (X1(IL) .GT. ALPHA(IL,K-1)) X(IL) = X2(IL)
                  IF (Y1(IL) .GT. RCORR(IL,K-1)) Y(IL) = Y2(IL)
               END IF

! Treatment of cloudy cells
               IF (X(IL) .GT. 1.0-AVG_CF(IL,K)) THEN ! Generate cloud in this layer

                  IF (LPPH) THEN ! Homogeneous clouds
                     ZCW = 1.0
                  ELSE
! Horizontally variable clouds:
! Determine ZCW = ratio of cloud condensate miximg ratio QC for this cell to
! its mean value for all cloudy cells in this layer.
! Use bilinear interpolation of ZCW tabulated in array XCW as a function of
!    * cumulative probability Y
!    * relative standard deviation SIGMA
! Take care that the definition of RIND2 is consistent with subroutine
! TABULATE_XCW

                        RIND1 = Y(IL) * (N1 - 1) + 1.0
                        IND1  = MAX(1, MIN(INT(RIND1), N1-1))
                        RIND1 = RIND1 - IND1
                        RIND2 = 40.0 * SIGMA_QCW(IL,K) - 3.0
                        IND2  = MAX(1, MIN(INT(RIND2), N2-1))
                        RIND2 = RIND2 - IND2

                        ZCW = (1.0-RIND1) * (1.0-RIND2) * XCW(IND1,IND2) &
                          + (1.0-RIND1) * RIND2       * XCW(IND1,IND2+1) &
                          + RIND1 * (1.0-RIND2)       * XCW(IND1+1,IND2) &
                        + RIND1 * RIND2             * XCW(IND1+1,IND2+1)
                  END IF

! A horizontally constant IWC/LWC ratio is assumed for each layer so far
                  L_CLD(IL)             = .TRUE.
                  LIQWCIN_S(IL,K,I_LOC(IL)) = ZCW * AVG_QLW(IL,K)
                  ICEWCIN_S(IL,K,I_LOC(IL)) = ZCW * AVG_QIW(IL,K)
               END IF

            END DO              ! IL
         END DO                 ! K
! Need to check if a cloudy subcolumn was generated
         DO IL = IL1, IL2
            IF (L_CLD(IL)) THEN
               I_LOC(IL) = I_LOC(IL) + 1
               L_CLD(IL) = .FALSE.
            END IF
         END DO
      END DO                    ! I

! Record the number of cloudy subcolumns generated
      DO IL = IL1, IL2
         NCLDY(IL) = I_LOC(IL) - 1
      END DO ! IL
      RETURN
      END SUBROUTINE MCICA_CLD_GENERATOR

      SUBROUTINE PREP_MCICA(RLC_CF, RLC_CW, SIGMA_QCW, NUAGE, ILG, IL1, &
                            IL2, LAY)

! --------------------------------------------------------------------
! Subroutine to define the horizontal variability and overlap for the
! stochastic cloud generator.
! --------------------------------------------------------------------

      implicit none
#include <arch_specific.hf>

!
! INPUT DATA
!

      INTEGER, INTENT(IN) :: &
       ILG, &
       IL1, &
       IL2, &
       LAY

      REAL, INTENT(IN) :: &
       NUAGE(ILG,LAY)     ! Cloud amount

!
! OUTPUT DATA
!

      REAL, INTENT(OUT) :: &
       RLC_CF(ILG,LAY),    &! Cloud fraction decorrelation length               (km)
       RLC_CW(ILG,LAY),    &! Cloud condensate decorrelation length             (km)
       SIGMA_QCW(ILG,LAY) ! Normalized standard deviation of cloud condensate (unitless)

!
! LOCAL DATA
!

      INTEGER :: &
       IL, &
       KK

      REAL :: &
       ANU

! ZERO OUT THE OUTPUT ARRAYS
      DO KK = 1, LAY
         DO IL = IL1, IL2
            SIGMA_QCW(IL,KK) = 0.0
            RLC_CF(IL,KK)    = 0.0
            RLC_CW(IL,KK)    = 0.0
         END DO ! IL
      END DO ! KK

! FOR NOW SET THE DECORRELATION LENGTHS TO BE REASONABLE ESTIMATES
! UPDATE WHEN HAVE SOME CLUE ABOUT HOW TO DO THIS BASED ON LARGE-SCALE
! VARIABLES

      DO KK = 1, LAY
         DO IL = IL1, IL2
            RLC_CF(IL,KK) = 2.0
            RLC_CW(IL,KK) = 1.0
         END DO
      END DO

! DEFINE THE HORIZONTAL VARIABILITY USING METHOD CURRENTLY USED IN GEM
! WILL CHANGE SOON

      DO KK = 1, LAY
         DO IL = IL1, IL2
            IF (NUAGE(IL,KK) .LE. 0.9)                              THEN
               ANU  =  1.0
            ELSEIF (NUAGE(IL,KK) .GT. 0.9 .AND. &
                    NUAGE(IL,KK) .LT. 1.0)                          THEN
               ANU  =  2.0
            ELSEIF (NUAGE(IL,KK) .GT. 0.0)                          THEN
               ANU  =  4.0
            ENDIF

! THE NORMALIZED VARIABILITY (STD DEV/MEAN) CAN BE APPROXIMATED AS
! 1.0/SQRT(ANU)
            IF (NUAGE(IL,KK) .GT. 0.0) THEN
               SIGMA_QCW(IL,KK) = 1.0/SQRT(ANU)
            ELSE
               SIGMA_QCW(IL,KK) = 0.0
            END IF
         END DO ! IL
      END DO ! KK

      RETURN

      END SUBROUTINE PREP_MCICA

      SUBROUTINE McICA_CLD_GEN(NUAGE, LIQWCIN, ICEWCIN, RLC_CF, RLC_CW, &
                               SIGMA_QCW, T, S, PS, ILG, IL1, IL2, LAY, &
                               NCLDY, LIQWCIN_S, ICEWCIN_S, CLDTOT)

! --------------------------------------------------------------------
! Driver to call stochastic cloud generator and produce diagnostic cloud
! properties that are consistent with McICA radiative transfer routine.
! --------------------------------------------------------------------

      implicit none
#include <arch_specific.hf>

!       EXTERNAL MCICA_CLD_GENERATOR

#include "mcica.cdk"
#include "tdpack_const.hf"

! Note: LAY    => Number of layers
! Note: NX_LOC => Number of subcolumns to generate

!
! PARAMETER
!

      REAL, PARAMETER :: &
       M2KM = 1.0/1000.0,  &! Convert meters to kilometers
       CUT  = 0.001

!
! INPUT DATA
!

      REAL, INTENT(IN) :: &
       NUAGE(ILG,LAY),        &! Column cloud fraction
       ICEWCIN(ILG,LAY),      &! Column in-cloud ice water mixing ratio profile    (kg/kg)
       LIQWCIN(ILG,LAY)      ! Column in-cloud liquid water mixing ratio profile (kg/kg)

      INTEGER, INTENT(IN) ::  &! Counters and array sizes
       ILG, &
       IL1, &
       IL2, &
       LAY

      REAL, INTENT(IN) :: &
       RLC_CF(ILG,LAY),     &! Cloud fraction decorrelation length               (km)
       RLC_CW(ILG,LAY),     &! Cloud condensate decorrelation length             (km)
       SIGMA_QCW(ILG,LAY),  &! Normalized standard deviation of cloud condensate (unitless)
       T(ILG,LAY),          &! Column temperature at layer midpoint              (K)
       PS(ILG),             &! Surface pressure                                  (Pa)
       S(ILG,LAY)          ! Column hybrid coord at layer midpoint             (unitless)

!
! OUTPUT DATA
!

      REAL, INTENT(OUT) :: &
       CLDTOT(ILG),                &! McICA vertical projected total cloud fraction  (unitless)
       ICEWCIN_S(ILG,LAY,NX_LOC),  &! Column ice water mixing ratio profile          (g/m^3)
       LIQWCIN_S(ILG,LAY,NX_LOC)  ! Column liquid water mixing ratio profile       (g/m^3)

      INTEGER,INTENT(OUT) :: &
       NCLDY(ILG)                 ! Number of cloudy subcolumns

!
! LOCAL DATA
!

      INTEGER :: &
       IL,        &! Counter over GCM columns
       II,        &! Counter over NX subcolumns
       KK        ! Counter over lay vertical layers

      REAL :: &
       RHO                     ! Density of air                                 (g/m^3)

      REAL :: &
       P,                       &! GCM column pressure at layer midpoint        (Pa)
       ROG                     ! Total gas constant/gravity

      REAL :: &
       DMULT, &
       Q_TOT

      REAL :: &
       X1(ILG)

      REAL :: &
       QI_PROF(ILG,LAY), &
       QC_PROF(ILG,LAY), &
       ZM(ILG,LAY), &
       DZ(ILG,LAY)

      INTEGER :: &
       K, &
       IND1

! Zero out fields
      DO IL = IL1,IL2
         CLDTOT(IL) = 0.0
         NCLDY(IL)  = 0
      END DO ! IL

      DO KK = 1, LAY
         DO IL = IL1,IL2
            QC_PROF(IL,KK) = 0.0
            QI_PROF(IL,KK) = 0.0
         END DO ! IL
      END DO ! KK

      DO II = 1, NX_LOC
         DO KK = 1 , LAY
            DO IL = IL1,IL2
                  ICEWCIN_S(IL,KK,II) = 0.0
                  LIQWCIN_S(IL,KK,II) = 0.0
            END DO
         END DO
      END DO

! Compute the heights of mid-layers

      ROG=RGASD/GRAV

      DO KK = 1, LAY
         DO IL = IL1, IL2
            P = S(IL,KK)*PS(IL)
            ZM(IL,KK) = ROG*T(IL,KK)*LOG(PS(IL)/P)*M2KM
         END DO
      END DO

! Convert the cloud condensate from kg/kg to g/m^3 (is the cloud condensate cloud or domain mean)?
      DO KK = 1, LAY
         DO IL = IL1, IL2
            P = S(IL,KK)*PS(IL)
! Compute layer height
            IF (NUAGE(IL,KK) .GT. CUT) THEN
! RHO in kg/m^3
               RHO            = 1000.0*P/(RGASD*T(IL,KK))
               DMULT          = RHO
               QI_PROF(IL,KK) = ICEWCIN(IL,KK)*DMULT
               QC_PROF(IL,KK) = LIQWCIN(IL,KK)*DMULT
            END IF
         END DO ! IL
      END DO ! KK

! Call cloud generator
      CALL MCICA_CLD_GENERATOR(ZM, NUAGE, QC_PROF, QI_PROF, &
                               RLC_CF, RLC_CW, SIGMA_QCW, &
                               ILG, IL1, IL2, LAY, &
                               ICEWCIN_S,LIQWCIN_S, NCLDY)

      DO IL = 1, ILG
         CLDTOT(IL) = REAL(NCLDY(IL))/REAL(NX_LOC)
      END DO !IL

      RETURN

      END SUBROUTINE McICA_CLD_GEN

