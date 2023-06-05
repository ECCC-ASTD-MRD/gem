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
!**S/R LIQWC
!
      SUBROUTINE LIQWC(LWC, SIG, T, PSOL, LMX, LEV, MM, SATUCO)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
!
      INTEGER LMX,LEV,MM,J,K
      REAL CLWC,GAMSAT,GAMW,GB
      REAL RHO,WC,WCCON,WCFAC,WCOBS,WCSTA
          REAL LWC(LMX,LEV), SIG(lmx,LEV), T(MM ,LEV), REC_T(MM ,LEV)
      REAL P(LMX,LEV),GA(LMX,LEV),EXP1_CAPPA,EXP2_CAPPA
      REAL PSOL(LMX)   ,HT(LMX,LEV)
      LOGICAL LHIG(LMX,LEV), SATUCO
      REAL P0, REC_CAPPA, REC_RGASD

!Author
!          L.Garand (1989)
!
!Revision
! 001      G.Pellerin(Mar90)Standard documentation
! 002      N. Brunet  (May91)
!                New version of thermodynamic functions
!                and file of constants
! 003      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 004      R. Benoit (Aug 93) Local Sigma
! 005      L. Garand (Apr 1995) Output in kg/kg to input cldoptx
! 006      M. Lepine (March 2003) -  CVMG... Replacements
! 007      D. Talbot (June 2003) - IBM conversion
!                - calls to exponen4 (to calculate power function '**')
!                - calls to optimized routine mfoqst3
!                - divisions replaced by reciprocals
! 008      B. Bilodeau (Aug 2003) - exponen4 replaces vspown1
!
!Object
!          to calculate the liquid water content (kg/kg) as a
!          function of temperature and the input parameters for
!          calulating solar radiation (Refer to Betts et
!          Harshvardan JGR 1987)
!
!Arguments
!
!          - Output -
! LWC      liquid water content (kg water per kg air)
!
!          - Input -
! SIG      sigma levels
! T        temperature in Kelvins
! PSOL     surface pressure (N/m**2)
! LMX      number of points to process
! LEV      number of layers
! MM       1st dimension of the temperature field
! SATUCO   .TRUE. if water/ice phase for saturation
!          .FALSE. if water phase only for saturation
!
!*

      P0=101325.
      REC_CAPPA=1./CAPPA
      REC_RGASD=1./RGASD
!   MULTIPLICATEUR POUR NUAGES CONVECTIFS
      WCFAC=1.
      DO J=1,LEV
      DO K=1,LMX
! . . . . MID LAYER PRESSURE IN [PASCAL].
              P(K,J) = PSOL(K) * SIG(k,J)
! . . . . LOCATE UPPER TROPOSPHERE ABOVE CONVECTIVE REGIONS.
              LHIG(k,j) = (SIG(k,J).GT.0.6)
! . . . . EVALUATE L.W.C. FROM PAPER BY BETTS AND HARSVARDAN (1987).
!         ASSUMING THAT ADIABATIC LIFTING OCCURE OVER HALF OF THE LAYER.
!         GAM SAT IS FROM CONVEC ROUTINE = T * (1 - GAM S/GAM D)
              HT(K,J) = 3.1364E+3 - 2.3682 * T(K,J)
!         On calcule la reciproque de T car T est souvent au denominateur plus bas
              REC_T(K,J)=1./T(K,J)
      ENDDO
      ENDDO

      IF(SATUCO)THEN
         CALL mfoqst3 (GA,T,P,LMX,LEV,LMX)
      ELSE
         CALL mfoqsa3 (GA,T,P,LMX,LEV,LMX)
      END IF

      DO J=1,LEV
      DO K=1,LMX
         exp1_cappa = (p0/p(K,J))**CAPPA
         exp2_cappa = (p(K,J)/p0)**CAPPA
              GA(K,J) = HT(K,J)*GA(K,J)/(CAPPA*T(K,J))
              GB = GA(K,J) * EPS1 * HT(K,J) * REC_T(K,J)
              GAMSAT = T(K,J) * (GB - GA(K,J)) / (1. + GB)
! . . . . GAMMA W = (D THETA/DP) ON THETA ES CONST (BETTS EQ(4))
              GAMW = EXP1_CAPPA*CAPPA/P(K,J)*GAMSAT
! . . . . AIR DENSITY [KG/M3].
              RHO = P(K,J)*REC_RGASD*REC_T(K,J)
! . . . . LWC = (CP*1000[G/KG]/L)*(T/THETA)*GAMSAT*RHO*DP*MIXING.
      CLWC=1000./HT(K,J)
      WCOBS = CLWC*EXP2_CAPPA*GAMW*RHO*1693.
!             UNITS [G/M3]. (DP = 30 MB AND MIXING = 0.56433)
! . . . . FIT OF MEAN WC FROM OBSERVATIONS (FEIGELSON, 1977)
!             WC OBS = (4.35E-5 * T(K,J) - 0.0174) * T(K,J) + 1.734
!             WC OBS = AMAX1 (WC OBS, 0.01)
! . . . . FIT OF MINIMUM LIQUID WATER CONTENT (HEYMSFIELD, 1977)
!             WC LOW = 8.64E-18 * EXP (0.1462 * T(K,J))
! . . . . SET LIQUID WATER CONTENT IN STRATIFORM CLOUDS.
!             WC STA = AMIN1 (WC OBS, WC LOW)
              WCSTA = WCOBS
! . . . . SET LIQUID WATER CONTENT IN CONVECTIVE REGIONS.
!             WC CON = AMIN1 (WC OBS, WC LOW * WC FAC)
              WCCON = WCOBS * WCFAC
! . . . . SELECT STRATIFORM OR CONVECTIVE VALUES OF LWC
              if ( LHIG(k,j) ) then
                 WC = WCCON
              else
                 WC = WCSTA
              endif
!  MULT PAR 0.5 ASSUMANT DEMI-EFFICACITE
      LWC(K,J)=WC  * 0.001 *.5 / RHO
      ENDDO
      ENDDO
      RETURN
      END
