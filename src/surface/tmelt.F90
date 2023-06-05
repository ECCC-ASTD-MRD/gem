      SUBROUTINE TMELT(ZSNOW,TSNOW,QMELT,R,TR,GZERO,RALB,            &
                       HMFN,HTCS,HTC,FI,HCPSNO,RHOSNO,WSNOW,         &
                       ISAND,IG,ILG,IL1,IL2,JL) 
!                                                                                  
!     * JAN 06/09 - D.VERSEGHY/M.LAZARE. SPLIT 100 LOOP INTO TWO.
!     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.
!     * SEP 24/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
!     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.
!     * JUN 20/97 - D.VERSEGHY. CLASS - VERSION 2.7.
!     *                         MODIFICATIONS TO ALLOW FOR VARIABLE
!     *                         SOIL PERMEABLE DEPTH.
!     * JAN 02/95 - D.VERSEGHY. CLASS - VERSION 2.5.
!     *                         COMPLETION OF ENERGY BALANCE
!     *                         DIAGNOSTICS.
!     * AUG 18/95 - D.VERSEGHY. CLASS - VERSION 2.4.
!     *                         REVISIONS TO ALLOW FOR INHOMOGENEITY
!     *                         BETWEEN SOIL LAYERS.
!     * JUL 30/93 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.2.
!     *                                  NEW DIAGNOSTIC FIELDS.
!     * APR 24/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.
!     *                                  REVISED AND VECTORIZED CODE
!     *                                  FOR MODEL VERSION GCM7.
!     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -
!     *                         CLASS VERSION 2.0 (WITH CANOPY).
!     * APR 11/89 - D.VERSEGHY. MELTING OF SNOWPACK.
!
      IMPLICIT NONE
!
!     * INTEGER CONSTANTS.
!
      INTEGER IG,ILG,IL1,IL2,JL,I
!
!     * INPUT/OUTPUT ARRAYS.
!
      REAL HTC   (ILG,IG)

      REAL ZSNOW (ILG),   TSNOW (ILG),   QMELT (ILG),   R     (ILG),  &
           TR    (ILG),   GZERO (ILG),   RALB  (ILG),   HMFN  (ILG),  &
           HTCS  (ILG)
!
!     * INPUT ARRAYS.
!
      REAL FI    (ILG),   HCPSNO(ILG),   RHOSNO(ILG),   WSNOW (ILG)
!
      INTEGER             ISAND (ILG,IG)
!                                                                                 
!     * TEMPORARY VARIABLES.
!
      REAL HADD,HCONV,ZMELT,RMELT,RMELTS,TRMELT
!
!     * COMMON BLOCK PARAMETERS.
!
      REAL DELT,TFREZ,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,  &
           SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
!
      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,       &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,        &
                      TCGLAC,CLHMLT,CLHVAP
!-----------------------------------------------------------------------
      DO 100 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              IF(QMELT(I).GT.0. .AND. ZSNOW(I).GT.0.)           THEN
                  HTCS(I)=HTCS(I)-FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*    &
                          ZSNOW(I)/DELT
                  HADD =QMELT(I)*DELT                         
                  HCONV=(0.0-TSNOW(I))*HCPSNO(I)*ZSNOW(I) +            &
                                CLHMLT*RHOSNO(I)*ZSNOW(I)      
                  IF(HADD.LE.HCONV)                     THEN
                      ZMELT=HADD/((0.0-TSNOW(I))*HCPSNO(I)+            &
                            CLHMLT*RHOSNO(I))                           
                      RMELTS=ZMELT*RHOSNO(I)/(RHOW*DELT)         
                      RMELT=RMELTS
                      TRMELT=0.0                            
                      ZSNOW(I)=ZSNOW(I)-ZMELT                 
                      HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/ &
                          (RHOW*ZSNOW(I))
                      HTCS (I)=HTCS(I)-FI(I)*(QMELT(I)-CLHMLT*RMELT*   &
                               RHOW)
                  ELSE                                             
                      RMELTS=ZSNOW(I)*RHOSNO(I)/RHOW                  
                      RMELT=RMELTS+WSNOW(I)/RHOW
                      HADD=HADD-HCONV                                
                      TRMELT=HADD/(HCPW*RMELT)                      
                      RMELT=RMELT/DELT                             
                      RMELTS=RMELTS/DELT
                      ZSNOW (I)=0.0    
                      HCPSNO(I)=0.0
                      TSNOW (I)=0.0   
                      WSNOW (I)=0.0
                      HTCS (I)=HTCS(I)-FI(I)*(QMELT(I)-CLHMLT*RMELTS*  &
                               RHOW-HADD/DELT)
                  ENDIF                                            
                  HMFN (I)=HMFN(I)+FI(I)*CLHMLT*RMELTS*RHOW
                  TR   (I)=(R(I)*TR(I)+RMELT*TRMELT)/(R(I)+RMELT)
                  R    (I)=R(I)+RMELT
                  QMELT(I)=0.0
                  HTCS(I)=HTCS(I)+FI(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*    &
                          ZSNOW(I)/DELT
              ENDIF
              RALB(I)=R(I)
          ENDIF
  100 CONTINUE                                            
!
      DO 200 I=IL1,IL2
          IF(FI(I).GT.0.)                                          THEN
              IF(QMELT(I).GT.0. .AND. ISAND(I,1).GT.-4)      THEN
                  GZERO(I)=GZERO(I)+QMELT(I)
                  HTCS (I)=HTCS(I)-FI(I)*QMELT(I)
                  HTC(I,1)=HTC(I,1)+FI(I)*QMELT(I)
              ENDIF
              RALB(I)=R(I)
          ENDIF
  200 CONTINUE                                               
!                                                                                  
      RETURN                            
      END        
