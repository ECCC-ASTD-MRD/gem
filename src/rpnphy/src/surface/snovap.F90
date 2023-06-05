      SUBROUTINE SNOVAP(RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAP,QFN,QFG,HTCS, &
                        WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,           &
                        FI,R,S,RHOSNI,WSNOW,ILG,IL1,IL2,JL)
!
!     * AUG 25/11 - D.VERSEGHY. CORRECT CALCULATION OF TRUNOF
!     *                         AND TOVRFL.
!     * FEB 22/07 - D.VERSEGHY. NEW ACCURACY LIMITS FOR R AND S.
!     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
!     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
!     * JUL 26/02 - D.VERSEGHY. CHANGE RHOSNI FROM CONSTANT TO
!     *                         VARIABLE.
!     * APR 11/01 - M.LAZARE.   CHECK FOR EXISTENCE OF SNOW BEFORE
!     *                         PERFORMING CALCULATIONS.
!     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
!     *                         PASS IN NEW "CLASS4" COMMON BLOCK.
!     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.
!     *                         COMPLETION OF ENERGY BALANCE
!     *                         DIAGNOSTICS.
!     * AUG 16/95 - D.VERSEGHY. CLASS - VERSION 2.4.
!     *                         INCORPORATE DIAGNOSTIC ARRAY "WLOST". 
!     * DEC 22/94 - D.VERSEGHY. CLASS - VERSION 2.3.
!     *                         ADDITIONAL DIAGNOSTIC CALCULATION -
!     *                         UPDATE HTCS.
!     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
!     *                                  NEW DIAGNOSTIC FIELDS.
!     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
!     *                                  REVISED AND VECTORIZED CODE
!     *                                  FOR MODEL VERSION GCM7.
!     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
!     *                         CLASS VERSION 2.0 (WITH CANOPY).
!     * APR 11/89 - D.VERSEGHY. SUBLIMATION FROM SNOWPACK.
!                                          
      IMPLICIT NONE
!
!     * INTEGER CONSTANTS.
!
      INTEGER ILG,IL1,IL2,JL,I
!
!     * INPUT/OUTPUT ARRAYS.
!
      REAL RHOSNO(ILG),   ZSNOW (ILG),   HCPSNO(ILG),   TSNOW (ILG),    &
           EVAP  (ILG),   QFN   (ILG),   QFG   (ILG),   HTCS  (ILG),    &
           WLOST (ILG),   TRUNOF(ILG),   RUNOFF(ILG),   TOVRFL(ILG),    &
           OVRFLW(ILG)
!
!     * INPUT ARRAYS.
!
      REAL FI    (ILG),   R     (ILG),   S     (ILG),   RHOSNI(ILG),    &
           WSNOW (ILG)   
!
!     * TEMPORARY VARIABLES.
!
      REAL ZADD,ZLOST,ZREM
!
!     * COMMON BLOCK PARAMETERS.
!
      REAL DELT,TFREZ,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,      &
           SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
!                                       
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,           &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,            &
                      TCGLAC,CLHMLT,CLHVAP
!-----------------------------------------------------------------------
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0. .AND. (S(I).LT.1.0E-11 .OR. R(I).LT.1.0E-11)  &
                      .AND. ZSNOW(I).GT.0.)                       THEN
              HTCS(I)=HTCS(I)-FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*        &
                      ZSNOW(I)/DELT
              IF(EVAP(I).LT.0.)                             THEN 
                  ZADD=-EVAP(I)*DELT*RHOW/RHOSNI(I)
                  RHOSNO(I)=(ZSNOW(I)*RHOSNO(I)+ZADD*RHOSNI(I))/       &
                            (ZSNOW(I)+ZADD)                          
                  ZSNOW (I)=ZSNOW(I)+ZADD                
                  HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/     &
                      (RHOW*ZSNOW(I))
                  EVAP  (I)=0.0                       
              ELSE                                       
                  ZLOST=EVAP(I)*DELT*RHOW/RHOSNO(I)
                  IF(ZLOST.LE.ZSNOW(I))                     THEN 
                      ZSNOW(I)=ZSNOW(I)-ZLOST                  
                      EVAP (I)=0.0                            
                      HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/ &
                          (RHOW*ZSNOW(I))
                  ELSE                              
                      ZREM=(ZLOST-ZSNOW(I))*RHOSNO(I)/RHOW
                      ZSNOW(I)=0.0                                
                      HCPSNO(I)=0.0
                      EVAP(I)=ZREM*(CLHMLT+CLHVAP)/(CLHVAP*DELT)
                      WLOST(I)=WLOST(I)-ZREM*RHOW*CLHMLT/CLHVAP
                      IF(RUNOFF(I).GT.0. .OR. WSNOW(I).GT.0.)          &
                       TRUNOF(I)=(TRUNOF(I)*RUNOFF(I)+(TSNOW(I)+TFREZ)*&
                            WSNOW(I)/RHOW)/(RUNOFF(I)+WSNOW(I)/RHOW)
                      RUNOFF(I)=RUNOFF(I)+WSNOW(I)/RHOW
                      IF(OVRFLW(I).GT.0. .OR. WSNOW(I).GT.0.)          &
                       TOVRFL(I)=(TOVRFL(I)*OVRFLW(I)+(TSNOW(I)+TFREZ)*&
                            FI(I)*WSNOW(I)/RHOW)/(OVRFLW(I)+FI(I)*     &
                            WSNOW(I)/RHOW)
                      OVRFLW(I)=OVRFLW(I)+FI(I)*WSNOW(I)/RHOW
                      TSNOW(I)=0.0 
                      WSNOW(I)=0.0
                      QFN(I)=QFN(I)-FI(I)*ZREM*RHOW/DELT
                      QFG(I)=QFG(I)+FI(I)*EVAP(I)*RHOW
                  ENDIF    
              ENDIF 
              HTCS(I)=HTCS(I)+FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*  &
                      ZSNOW(I)/DELT
          ENDIF                                           
  100 CONTINUE
!                                                                                  
      RETURN                                      
      END        
