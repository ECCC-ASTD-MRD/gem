      SUBROUTINE SNOADD(ALBSNO,TSNOW,RHOSNO,ZSNOW,HCPSNO,HTCS,    &
                        FI,S,TS,RHOSNI,WSNOW,ILG,IL1,IL2,JL)
!
!     * NOV 17/11 - M.LAZARE.   CHANGE SNOW ALBEDO REFRESHMENT 
!     *                         THRESHOLD (SNOWFALL IN CURRENT
!     *                         TIMESTEP) FROM 0.005 TO 1.E-4 M. 
!     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
!     * SEP 10/04 - R.HARVEY/D.VERSEGHY. INCREASE SNOW ALBEDO 
!     *                         REFRESHMENT THRESHOLD; ADD
!     *                         "IMPLICIT NONE" COMMAND.
!     * JUL 26/02 - D.VERSEGHY. CHANGE RHOSNI FROM CONSTANT TO
!     *                         VARIABLE.
!     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
!     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
!     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
!     *                         COMPLETION OF ENERGY BALANCE
!     *                         DIAGNOSTICS.
!     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
!     *                                  NEW DIAGNOSTIC FIELDS.
!     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
!     *                                  REVISED AND VECTORIZED CODE
!     *                                  FOR MODEL VERSION GCM7.
!     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
!     *                         CLASS VERSION 2.0 (WITH CANOPY).
!     * APR 11/89 - D.VERSEGHY. ACCUMULATION OF SNOW ON GROUND.
!
      IMPLICIT NONE
!
!     * INTEGER CONSTANTS.
!
      INTEGER ILG,IL1,IL2,JL,I
!
!     * INPUT/OUTPUT ARRAYS.
!
      REAL ALBSNO(ILG),   TSNOW (ILG),   RHOSNO(ILG),   ZSNOW (ILG),  &
           HCPSNO(ILG),   HTCS  (ILG)
!
!     * INPUT ARRAYS.
!
      REAL FI    (ILG),   S     (ILG),   TS    (ILG),   RHOSNI(ILG),  &
           WSNOW (ILG)

!     * TEMPORARY VARIABLES.
!
      REAL SNOFAL,HCPSNP
!                                                                                 
!     * COMMON BLOCK PARAMETERS.
!
      REAL DELT,TFREZ,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,         &
           SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,       &
           CLHVAP
!
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,         &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,          &
                      TCGLAC,CLHMLT,CLHVAP
!-----------------------------------------------------------------------
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0. .AND. S(I).GT.0.)   THEN
              HTCS  (I)=HTCS(I)-FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT
              SNOFAL=S(I)*DELT                      
              IF(SNOFAL.GE.1.E-4) THEN 
!   ipy test
!              IF(SNOFAL.GE.0.005)                               THEN 
                  ALBSNO(I)=0.84 
              ELSE IF(.NOT.(ZSNOW(I).GT.0.)) THEN
                  ALBSNO(I)=0.50  
              ENDIF    
              HCPSNP=HCPICE*RHOSNI(I)/RHOICE
              TSNOW (I)=((TSNOW(I)+TFREZ)*ZSNOW(I)*HCPSNO(I) +     &
                         (TS   (I)+TFREZ)*SNOFAL  *HCPSNP)/        &
                        (ZSNOW(I)*HCPSNO(I) + SNOFAL*HCPSNP) -     &
                         TFREZ
              RHOSNO(I)=(ZSNOW(I)*RHOSNO(I) + SNOFAL*RHOSNI(I))/   &
                        (ZSNOW(I)+SNOFAL)                          
              ZSNOW (I)=ZSNOW(I)+SNOFAL        
              HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I))
              HTCS  (I)=HTCS(I)+FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/DELT
          ENDIF                                                 
  100 CONTINUE
!                                                                                  
      RETURN 
      END    
