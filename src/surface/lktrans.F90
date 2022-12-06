      SUBROUTINE LKTRANS (CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,   &
                 BLAK,IL1,IL2,ILG,CQ1BI,CQ2BI,CQ3BI)
!======================================================================
!     * DEC  7/07 - M.MACKAY.  	COMPUTES LIGHT EXTINCTION COEFFICIENTS
!     *                         THIS ROUTINE KEPT FOR FUTURE DYNAMIC
!				CHANGES TO EXTINCTION
!
      IMPLICIT NONE
!
! ----* INPUT FIELDS *------------------------------------------------
!
      REAL,DIMENSION(ILG) :: BLAK, CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,CQ1BI,CQ2BI,CQ3BI
      INTEGER IL1,IL2,ILG
!
! ----* LOCAL VARIABLES *------------------------------------------------
!
      INTEGER I

!======================================================================
!                CQ1A,    CQ1B,    CQ2A,    CQ2B,    CQ3A,    CQ3B 
! Rayner,1980    0.54,    0.561,   0.30,    6.89,    0.16,    69.0 
!======================================================================
      DO 100 I=IL1,IL2
! FIXED WATER VALUES
      CQ1A(I)=0.5817
      CQ2A(I)=0.4183
      CQ2B(I)=6.89
      CQ3A(I)=0.0
      CQ3B(I)=69.0
! FIXED ICE VALUES (from Patterson and Hamblin, 1988, L&O)
      CQ1BI(I)=1.5
!mdm  CQ1BI(I)=3.75    !test - high extinction for snow ice
      CQ2BI(I)=20.0
      CQ3BI(I)=69.0

!======================================================================
! CQ1B NOW READ IN .INI FILE
!----------------------------------------------------------------------
      CQ1B(I)=BLAK(I)

100   CONTINUE

      RETURN
      END
