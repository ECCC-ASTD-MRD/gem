      SUBROUTINE NDRCOEFL (CDMN,CDHN,VA,ZREFM,ZREFH,FICE,FLS,ILG,IL1,IL2)
!=======================================================================
!     * AUG 14/18 - M.MACKAY.  	NEUTRAL TURBULENT TRANSFER COEFFICIENTS
!     *                         BASED ON DRCOEFL.F90
!=======================================================================
!
      IMPLICIT NONE
!
! ----* INPUT FIELDS *------------------------------------------------
!
      INTEGER ILG,IL1,IL2
      REAL,DIMENSION(ILG) :: VA,ZREFM,ZREFH,FICE,FLS
!
! ----* OUTPUT FIELDS *------------------------------------------------
!
      REAL,DIMENSION(ILG) :: CDHN, CDMN
!
! ----* CLASS COMMON BLOCKS *------------------------------------------
!
      REAL RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
!
! ----* LOCAL VARIABLES *---------------------------------------------
!
      INTEGER I,J,K
      REAL G,TOL,CHARN,RESID,VSQ,Cold,PI, CDMNW,CDHNW,CDMNI,CDHNI,XICE
!
! ----* LOCAL PARAMETERS *--------------------------------------------
!
      G=GRAV
      TOL=1.0E-8
      CHARN=0.0175
      PI=3.14159
!
!======================================================================
      DO 100 I=IL1,IL2
!-----------------------------------------------------------------------
! NEUTRAL DRAG COEFFICIENTS  
!  Iterative routine to solve for CDN based on Charnock relation for 
!  roughness
!
        CDHNW=1.35E-3
        CDMNW=1.0E-3        !initial trial value
        VSQ=VA(I)*VA(I)
        RESID=999.0
        DO 
          IF (ABS(RESID) .LE. TOL ) EXIT
          Cold=CDMNW
          CDMNW=VKC*VKC/(LOG(ZREFM(I)*G/(CHARN*Cold*VSQ))*LOG(ZREFM(I)*G/  &
           (CHARN*Cold*VSQ)))
          RESID=CDMNW-Cold
        END DO   
!
        CDMNI=(VKC/(LOG(ZREFM(I)/0.002)))**2
        CDHNI=(VKC/(LOG(ZREFH(I)/0.00067)))**2
        IF(FICE(I).GT.(FLS(I)+0.001))                THEN
           XICE=MAX(FICE(I)-FLS(I),0.0)
           CDMN(I)=(XICE*CDMNI+(1.0-FICE(I))*CDMNW)/(1.0-FLS(I))
           CDHN(I)=(XICE*CDHNI+(1.0-FICE(I))*CDHNW)/(1.0-FLS(I))
        ELSE
           CDMN(I)=CDMNW
           CDHN(I)=CDHNW
        ENDIF
!-----------------------------------------------------------------------
100   CONTINUE

      RETURN 
      END        
