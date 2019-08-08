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
!**S/P LSCTROL
!
      SUBROUTINE LSCTROL ( ilab, OMEGAP, SIGMA, NI, NK )
      implicit none
!!!#include <arch_specific.hf>
!
!
      INTEGER NI , NK
      INTEGER ilab(NI,NK)
      REAL OMEGAP(NI,NK), SIGMA(NI,NK)
!
!Author
!          Bernard Bilodeau
!
!Revision
! 001      B. Bilodeau (Jan 2001) - Automatic arrays
!
!Object
!
!Arguments
!
!          - Output -
! ilab     label array from Kuo schemes
!
!          - Input -
! OMEGAP   vertical velocity in pressure coordinates
! SIGMA    sigma levels
! NI       1st horizontal dimension
! NK       vertical dimension
!
!Notes
!*
      LOGICAL LO,LO1
      INTEGER JK,JL
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      INTEGER, dimension(NI) :: KM1
      INTEGER, dimension(NI) :: KM2
      INTEGER, dimension(NI) :: KP1
      INTEGER, dimension(NI) :: KP2
!
      REAL, dimension(NI) :: SIGMAK1
      REAL, dimension(NI) :: SIGMAK2
      REAL, dimension(NI) :: KWW1
      REAL, dimension(NI) :: KWW2
!
!***********************************************************************
!
!
!     ------------------------------------------------------------------
!
!

!     The moisture accession is set to zero for all levels when the
!     vertical velocity in pressure coordinates OMEGAP is positive
!     (downward) at both sigma levels SIGMAK1 and SIGMAK2.
!
!     Since SIGMAK1 and 2 do not necessarly coincide with sigma levels
!     of the model, OMEGAP is linearly interpolated at levels SIGMAK1 and 2
!     using weighting factors KWW1 and KWW2. KP1 or KP2 and KM1 or KM2 are
!     the indices of the model's sigma levels from which the interpolation
!     takes place.

!
!     SIGMAK1 and SIGMAK2 are the sigma levels at which OMEGAP is tested
!
      DO JK=1,NK
         DO JL=1,NI
            ilab(jl,jk) = 1
         END DO
      END DO
!
      DO JL = 1,NI
         SIGMAK1(JL) = 0.9
         SIGMAK2(JL) = 0.7
         KM1    (JL) = NK
         KM2    (JL) = NK
      END DO
!
!     general case
      DO JK = 1,NK
!
         DO JL = 1,NI
!
            IF (SIGMA(JL,JK) .LE. SIGMAK1(JL)) KM1(JL) = JK
            IF (SIGMA(JL,JK) .LE. SIGMAK2(JL)) KM2(JL) = JK
!
            KP1(JL) = KM1(JL) + 1
            KP2(JL) = KM2(JL) + 1
!
            KWW1(JL)=(SIGMAK1(JL)-SIGMA(JL,KM1(JL)))/ &
                     (SIGMA(JL,KP1(JL))-SIGMA(JL,KM1(JL)))
!
            KWW2(JL)=(SIGMAK2(JL)-SIGMA(JL,KM2(JL)))/ &
                     (SIGMA(JL,KP2(JL))-SIGMA(JL,KM2(JL)))
!
         END DO
!
      END DO
!
!     special cases
      DO JL = 1,NI
!
         IF (SIGMA(JL,1).GT.SIGMAK1(JL)) THEN
            SIGMAK1(JL) = SIGMA(JL,1)
            KM1(JL)     = 1
            KP1(JL)     = 1
            KWW1(JL)    = 1.0
         ENDIF

         IF (SIGMA(JL,1).GT.SIGMAK2(JL)) THEN
            SIGMAK2(JL) = SIGMA(JL,1)
            KM2(JL)     = 1
            KP2(JL)     = 1
            KWW2(JL)    = 1.0
         ENDIF

         IF ( KM1(JL)  .EQ. NK )         THEN
            KP1(JL)     = KM1(JL)
            KWW1(JL)    = 1.0
         ENDIF

         IF ( KM2(JL)  .EQ. NK )         THEN
            KP2(JL)     = KM2(JL)
            KWW2(JL)    = 1.0
         ENDIF
!
      END DO
!
!
      DO JK=1,NK
         DO JL=1,NI
            LO=(OMEGAP(JL,KP1(JL))*KWW1(JL) + &
                OMEGAP(JL,KM1(JL))*(1-KWW1(JL))).GT.0.
            LO1=(OMEGAP(JL,KP2(JL))*KWW2(JL)+ &
                 OMEGAP(JL,KM2(JL))*(1-KWW2(JL))).GT.0.
!
            if( LO.and.LO1 ) ilab(jl,jk) = 0
!
         END DO
      END DO
!
      RETURN
      END
