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

subroutine MCICA_CLD_GENERATOR(ZF, AVG_CF, AVG_QLW,AVG_QIW, &
     RLC_CF, RLC_CW,SIGMA_QCW, &
     ILG, IL1, IL2, LAY, &
     ICEWCIN_S,LIQWCIN_S, NCLDY)
   use tdpack_const, only: GRAV, RGASD
   implicit none
!!!#include <arch_specific.hf>
#include "mcica.cdk"
   !#include "xcwdata.cdk"

   ! INPUT DATA
   !

   integer, intent(IN) :: &
        ILG, &
        IL1, &
        IL2, &
        LAY

   real, intent(IN) :: &
        ZF(ILG,LAY),             &! Full-level (layer midpoint) altitude (km)
        AVG_CF(ILG,LAY),         &! Cloud fraction for each layer
        AVG_QLW(ILG,LAY),        &! Cloud mean liquid water mixing ratio (g/m^3)
        AVG_QIW(ILG,LAY),        &! Cloud mean ice mixing ratio          (g/m^3)
        SIGMA_QCW(ILG,LAY),      &! Normalized cloud condensate std. dev.
        RLC_CF(ILG,LAY),         &! Cloud fraction decorrelation length
        RLC_CW(ILG,LAY)         ! Cloud condensate decorrelation length

   !
   ! OUTPUT DATA
   !

   real,intent(OUT) :: &
        ICEWCIN_S(ILG,LAY,NX_LOC),      &! Column ice water mixing ratio profile    (g/m^3)
        LIQWCIN_S(ILG,LAY,NX_LOC)      ! Column liquid water mixing ratio profile (g/m^3)

   integer, intent(OUT) :: &
        NCLDY(ILG)                     ! Number of cloudy subcolumns

   !@Author Jason Cole, MSC/Cloud Physics (Summer 2005)
   !@Object  PART OF THE ISCCP CLOUD SIMULATOR PACKAGE
   !        Generate a NX column by LAY layer cloud field. Input profiles
   !        must be from top of the model to the bottom as is the output.
   !        A Method and model evaluation is described in:

   !        Raisanen et al, 2004, Stochastic generation of subgrid-scale
   !        cloudy columns for large-scale models, Quart. J. Roy. Meteorol.
   !        Soc., 2004 , 130 , 2047-2068
   !@Arguments
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

   real, parameter :: &
        CUT = 0.001

   integer :: &
        IL,          &! Counter over ILG GCM columns
        I,           &! Counter
        K,           &! Counter over LAY vertical layers
        K_TOP,       &! Index of top most cloud layer
        K_BASE,      &! Index of lowest cloud layer
        IND1,        &! Index in variability calculation
        IND2        ! Index in variability calculation

   real :: &
        ALPHA(ILG,LAY),  &! Fraction of maximum/random cloud overlap
        RCORR(ILG,LAY),  &! Fraction of maximum/random cloud condensate overlap
        RIND1,           &! Real index in variability calculation
        RIND2,           &! Real index in variability calculation
        ZCW             ! Ratio of cloud condensate miximg ratio for this cell to its layer cloud-mean value

   real ::  &! Random number vectors
        X(ILG), &
        Y(ILG), &
        X1(ILG), &
        Y1(ILG), &
        X2(ILG), &
        Y2(ILG)

   integer :: &
        I_LOC(ILG) ! Place the new subcolumns into the arrays starting from the front

   logical :: &
        L_TOP, &
        L_BOT, &
        L_CLD(ILG)

! Initialize the arrays

      do IL = IL1, IL2
         I_LOC(IL) = 1
         L_CLD(IL) = .false.
      end do

! Find uppermost cloudy layer

      L_TOP = .false.
      K_TOP = 0
      do K=1,LAY
         do IL = IL1, IL2
            K_TOP = K
            if (AVG_CF(IL,K) .gt. CUT) L_TOP = .true.
         end do ! IL
         if (L_TOP) exit
      end do ! K

! If no cloudy layers in any GCM column in this group, exit

      if (K_TOP .eq. 0) then
         return
      end if

! Find lowermost cloudy layer

      L_BOT = .false.

      K_BASE = 0
      do K=LAY,1,-1
         do IL = IL1, IL2
            K_BASE = K
            if (AVG_CF(IL,K) .gt. CUT) L_BOT = .true.
         end do ! IL
         if (L_BOT) exit
      end do ! K

! Calculate overlap factors ALPHA for cloud fraction and RCORR for cloud
! condensate based on layer midpoint distances and decorrelation depths

      do K=1, LAY !K_TOP,K_BASE-1
         do IL = IL1, IL2
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
         end do ! IL
      end do ! K

      call VSEXP(ALPHA,ALPHA,LAY*(IL2-IL1+1))
      call VSEXP(RCORR,RCORR,LAY*(IL2-IL1+1))

      do I=1,NX_LOC

! Generate all subcolumns for latitude chain

         do K = K_TOP, K_BASE

            if (K .eq. K_TOP) then
               call random_number(X)
               call random_number(Y)
            end if

            call random_number(X1)
            call random_number(Y1)

            call random_number(X2)
            call random_number(Y2)

            do IL = IL1, IL2

! Maximum-random overlap
               if (LMAXRAN) then
                  if (X(IL) .lt. 1.0-AVG_CF(IL,K-1)) then !It is clear above
                     X(IL) = X1(IL) * (1.0 - AVG_CF(IL,K-1))
                     Y(IL) = Y1(IL)
                  end if
! Generalized overlap
               else
                  if (X1(IL) .gt. ALPHA(IL,K-1)) X(IL) = X2(IL)
                  if (Y1(IL) .gt. RCORR(IL,K-1)) Y(IL) = Y2(IL)
               end if

! Treatment of cloudy cells
               if (X(IL) .gt. 1.0-AVG_CF(IL,K)) then ! Generate cloud in this layer

                  if (LPPH) then ! Homogeneous clouds
                     ZCW = 1.0
                  else
! Horizontally variable clouds:
! Determine ZCW = ratio of cloud condensate miximg ratio QC for this cell to
! its mean value for all cloudy cells in this layer.
! Use bilinear interpolation of ZCW tabulated in array XCW as a function of
!    * cumulative probability Y
!    * relative standard deviation SIGMA
! Take care that the definition of RIND2 is consistent with subroutine
! TABULATE_XCW

                        RIND1 = Y(IL) * (N1 - 1) + 1.0
                        IND1  = max(1, min(int(RIND1), N1-1))
                        RIND1 = RIND1 - IND1
                        RIND2 = 40.0 * SIGMA_QCW(IL,K) - 3.0
                        IND2  = max(1, min(int(RIND2), N2-1))
                        RIND2 = RIND2 - IND2

                        ZCW = (1.0-RIND1) * (1.0-RIND2) * XCW(IND1,IND2) &
                          + (1.0-RIND1) * RIND2       * XCW(IND1,IND2+1) &
                          + RIND1 * (1.0-RIND2)       * XCW(IND1+1,IND2) &
                        + RIND1 * RIND2             * XCW(IND1+1,IND2+1)
                  end if

! A horizontally constant IWC/LWC ratio is assumed for each layer so far
                  L_CLD(IL)             = .true.
                  LIQWCIN_S(IL,K,I_LOC(IL)) = ZCW * AVG_QLW(IL,K)
                  ICEWCIN_S(IL,K,I_LOC(IL)) = ZCW * AVG_QIW(IL,K)
               end if

            end do              ! IL
         end do                 ! K
! Need to check if a cloudy subcolumn was generated
         do IL = IL1, IL2
            if (L_CLD(IL)) then
               I_LOC(IL) = I_LOC(IL) + 1
               L_CLD(IL) = .false.
            end if
         end do
      end do                    ! I

! Record the number of cloudy subcolumns generated
      do IL = IL1, IL2
         NCLDY(IL) = I_LOC(IL) - 1
      end do ! IL
      return
      end subroutine MCICA_CLD_GENERATOR

      subroutine PREP_MCICA(RLC_CF, RLC_CW, SIGMA_QCW, NUAGE, ILG, IL1, &
                            IL2, LAY)

! --------------------------------------------------------------------
! Subroutine to define the horizontal variability and overlap for the
! stochastic cloud generator.
! --------------------------------------------------------------------

      implicit none
!!!#include <arch_specific.hf>

!
! INPUT DATA
!

      integer, intent(IN) :: &
       ILG, &
       IL1, &
       IL2, &
       LAY

      real, intent(IN) :: &
       NUAGE(ILG,LAY)     ! Cloud amount

!
! OUTPUT DATA
!

      real, intent(OUT) :: &
       RLC_CF(ILG,LAY),    &! Cloud fraction decorrelation length               (km)
       RLC_CW(ILG,LAY),    &! Cloud condensate decorrelation length             (km)
       SIGMA_QCW(ILG,LAY) ! Normalized standard deviation of cloud condensate (unitless)

!
! LOCAL DATA
!

      integer :: &
       IL, &
       KK

      real :: &
       ANU

! ZERO OUT THE OUTPUT ARRAYS
      do KK = 1, LAY
         do IL = IL1, IL2
            SIGMA_QCW(IL,KK) = 0.0
            RLC_CF(IL,KK)    = 0.0
            RLC_CW(IL,KK)    = 0.0
         end do ! IL
      end do ! KK

! FOR NOW SET THE DECORRELATION LENGTHS TO BE REASONABLE ESTIMATES
! UPDATE WHEN HAVE SOME CLUE ABOUT HOW TO DO THIS BASED ON LARGE-SCALE
! VARIABLES

      do KK = 1, LAY
         do IL = IL1, IL2
            RLC_CF(IL,KK) = 2.0
            RLC_CW(IL,KK) = 1.0
         end do
      end do

! DEFINE THE HORIZONTAL VARIABILITY USING METHOD CURRENTLY USED IN GEM
! WILL CHANGE SOON

      do KK = 1, LAY
         do IL = IL1, IL2
            if (NUAGE(IL,KK) .le. 0.9)                              then
               ANU  =  1.0
            elseif (NUAGE(IL,KK) .gt. 0.9 .and. &
                    NUAGE(IL,KK) .lt. 1.0)                          then
               ANU  =  2.0
            elseif (NUAGE(IL,KK) .gt. 0.0)                          then
               ANU  =  4.0
            endif

! THE NORMALIZED VARIABILITY (STD DEV/MEAN) CAN BE APPROXIMATED AS
! 1.0/SQRT(ANU)
            if (NUAGE(IL,KK) .gt. 0.0) then
               SIGMA_QCW(IL,KK) = 1.0/sqrt(ANU)
            else
               SIGMA_QCW(IL,KK) = 0.0
            end if
         end do ! IL
      end do ! KK

  return
end subroutine PREP_MCICA



subroutine McICA_CLD_GEN(NUAGE, LIQWCIN, ICEWCIN, RLC_CF, RLC_CW, &
     SIGMA_QCW, T, S, PS, ILG, IL1, IL2, LAY, &
                               NCLDY, LIQWCIN_S, ICEWCIN_S, CLDTOT)
   use tdpack_const, only: GRAV, RGASD
   implicit none
!!!#include <arch_specific.hf>

! --------------------------------------------------------------------
! Driver to call stochastic cloud generator and produce diagnostic cloud
! properties that are consistent with McICA radiative transfer routine.
! --------------------------------------------------------------------

!       EXTERNAL MCICA_CLD_GENERATOR

#include "mcica.cdk"

! Note: LAY    => Number of layers
! Note: NX_LOC => Number of subcolumns to generate

!
! PARAMETER
!

      real, parameter :: &
       M2KM = 1.0/1000.0,  &! Convert meters to kilometers
       CUT  = 0.001

!
! INPUT DATA
!

      integer, intent(IN) ::  &! Counters and array sizes
       ILG, &
       IL1, &
       IL2, &
       LAY

      real, intent(IN) :: &
       NUAGE(ILG,LAY),        &! Column cloud fraction
       ICEWCIN(ILG,LAY),      &! Column in-cloud ice water mixing ratio profile    (kg/kg)
       LIQWCIN(ILG,LAY)      ! Column in-cloud liquid water mixing ratio profile (kg/kg)


      real, intent(IN) :: &
       RLC_CF(ILG,LAY),     &! Cloud fraction decorrelation length               (km)
       RLC_CW(ILG,LAY),     &! Cloud condensate decorrelation length             (km)
       SIGMA_QCW(ILG,LAY),  &! Normalized standard deviation of cloud condensate (unitless)
       T(ILG,LAY),          &! Column temperature at layer midpoint              (K)
       PS(ILG),             &! Surface pressure                                  (Pa)
       S(ILG,LAY)          ! Column hybrid coord at layer midpoint             (unitless)

!
! OUTPUT DATA
!

      real, intent(OUT) :: &
       CLDTOT(ILG),                &! McICA vertical projected total cloud fraction  (unitless)
       ICEWCIN_S(ILG,LAY,NX_LOC),  &! Column ice water mixing ratio profile          (g/m^3)
       LIQWCIN_S(ILG,LAY,NX_LOC)  ! Column liquid water mixing ratio profile       (g/m^3)

      integer,intent(OUT) :: &
       NCLDY(ILG)                 ! Number of cloudy subcolumns

!
! LOCAL DATA
!

      integer :: &
       IL,        &! Counter over GCM columns
       II,        &! Counter over NX subcolumns
       KK        ! Counter over lay vertical layers

      real :: &
       RHO                     ! Density of air                                 (g/m^3)

      real :: &
       P,                       &! GCM column pressure at layer midpoint        (Pa)
       ROG                     ! Total gas constant/gravity

      real :: &
       DMULT

      real :: &
       QI_PROF(ILG,LAY), &
       QC_PROF(ILG,LAY), &
       ZM(ILG,LAY)

! Zero out fields
      do IL = IL1,IL2
         CLDTOT(IL) = 0.0
         NCLDY(IL)  = 0
      end do ! IL

      do KK = 1, LAY
         do IL = IL1,IL2
            QC_PROF(IL,KK) = 0.0
            QI_PROF(IL,KK) = 0.0
         end do ! IL
      end do ! KK

      do II = 1, NX_LOC
         do KK = 1 , LAY
            do IL = IL1,IL2
                  ICEWCIN_S(IL,KK,II) = 0.0
                  LIQWCIN_S(IL,KK,II) = 0.0
            end do
         end do
      end do

! Compute the heights of mid-layers

      ROG=RGASD/GRAV

      do KK = 1, LAY
         do IL = IL1, IL2
            P = S(IL,KK)*PS(IL)
            ZM(IL,KK) = ROG*T(IL,KK)*log(PS(IL)/P)*M2KM
         end do
      end do

! Convert the cloud condensate from kg/kg to g/m^3 (is the cloud condensate cloud or domain mean)?
      do KK = 1, LAY
         do IL = IL1, IL2
            P = S(IL,KK)*PS(IL)
! Compute layer height
            if (NUAGE(IL,KK) .gt. CUT) then
! RHO in kg/m^3
               RHO            = 1000.0*P/(RGASD*T(IL,KK))
               DMULT          = RHO
               QI_PROF(IL,KK) = ICEWCIN(IL,KK)*DMULT
               QC_PROF(IL,KK) = LIQWCIN(IL,KK)*DMULT
            end if
         end do ! IL
      end do ! KK

! Call cloud generator
      call MCICA_CLD_GENERATOR(ZM, NUAGE, QC_PROF, QI_PROF, &
                               RLC_CF, RLC_CW, SIGMA_QCW, &
                               ILG, IL1, IL2, LAY, &
                               ICEWCIN_S,LIQWCIN_S, NCLDY)

      do IL = 1, ILG
         CLDTOT(IL) = real(NCLDY(IL))/real(NX_LOC)
      end do !IL

   return
end subroutine McICA_CLD_GEN

