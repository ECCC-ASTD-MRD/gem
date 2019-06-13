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
      SUBROUTINE SUN_RADFIX1(SIGMA,COSAS1,COSAS2,PSOL,HEAT,TOIT,NI,NK)
!
      implicit none
!!!#include <arch_specific.hf>
!
      INTEGER NI,NK
      REAL    SIGMA(NI,NK)
      REAL    COSAS1(NI),COSAS2(NI),PSOL(NI),HEAT(NI,NK)
      LOGICAL TOIT(NI)

!Author    FRANCOIS LEMAY (CMDN), SEPTEMBRE 2003
!
!Revisions
! 001      A-M.LEDUC (Sept. 2003)   - IBM conversion
! 002      B. Bilodeau (March 2004) - Add double loop INSIDE the subroutine
!
!Object    THE SOLAR HEATING IS FIXED FOR THE STRATOSPHERIC LEVELS
!          (ABOVE 85MB). A NEW ROUTINE IS CREATED FOR THIS PURPOSE
!          RATHER THAN HAVING THE SAME RADFIX EXPRESSION IN
!          BOTH NEWRAD (non-radiative timesteps) AND SUN (radiative
!          timesteps) ROUTINES.
!Arguments
! SIGMA    sigma coordinate
! COSAS1,2 cosine of the solar zenith angle
! PSOL     surface pressure [N/m2]
! HEAT     heating rate [Kelvin/second]
! TOIT     .true. if the sigma value corresponds to the roof of the model

!Note:     the reason there is two cosines of the sun angle comes from
!          the newrad3.ftn code where the average of the cosine between
!          two radiative timesteps (cosas1) and the cosine itself (cosas2)
!          are involved in the radfix code. To get bit validation with the old code,
!          we introduced two variables but this is not very relevant.
!
!*
!
      INTEGER I,K
      REAL    SIGMA2,VV,VV2,A,ZZ,Y
      REAL    REC_86400

!     LES TAUX SONT FIXES ET NON CALCULES AU DESSUS DE 85 MB

      DO K=1,NK
      DO I=1,NI

      VV=AMAX1(SIGMA(I,K),0.002)
!     vertical coordinate independent of the surface pressure:
      VV2=VV*PSOL(I)*1.E-5
!     VV2=VV

!     EXPOSANT VARIANT AVEC HAUTEUR AU-DESSUS DE 85MB

      A=0.4

!     vertical coordinate independent of the surface pressure:
      SIGMA2=SIGMA(I,K)*PSOL(I)*1.E-5

      IF (SIGMA2.LT.0.020) A=.4+.274*(0.02-SIGMA2)*50
      IF (SIGMA2.LT.0.085.AND.SIGMA2.GT.0.0019) THEN
!        ZZ=COSAS1(I)**A
         ZZ=exp(A*log(COSAS1(I)))
      ELSE
         ZZ=COSAS1(I)
      ENDIF

      REC_86400=1./86400.

!     LE FIT EST DU SECOND DEGRE EN 1/SIGMA POUR SOLEIL AU NADIR
!     Y= (0.327+0.0642/VV -7.082E-5/VV**2)/86400. *ZZ

!     THE NEW FIT IS A LINEAR FUNCTION WITH A FLOOR VALUE EXPRESSED
!     IN TERMS OF PRESSURE RATHER THAN SIGMA (TOPOGR. INDEPENDANT)
!     Y= MAX((-75.*VV  + 5.25),2.0)/86400. *ZZ
!     Y= MAX((-75.*VV2 + 5.25),2.0)/86400. *ZZ
      Y= MAX((-75.*VV2 + 5.25),2.0)*REC_86400 *ZZ

!     SOLAR HEATING IS FIXED AT THE TOP
      IF (TOIT(K)) Y= (4.5)*REC_86400 *ZZ

      IF (COSAS2(I).LE.0.00001) Y=0.
      IF (SIGMA2.LT.0.085) HEAT(I,K)=Y

      HEAT(I,K)=AMAX1(HEAT(I,K),0.)
!
      END DO
      END DO

      RETURN
      END
