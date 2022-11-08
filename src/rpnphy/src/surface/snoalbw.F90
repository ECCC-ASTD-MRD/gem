      SUBROUTINE SNOALBW(ALBSNO,RHOSNO,ZSNOW,HCPSNO,TSNOW,           &
                         FI,S,RMELT,WSNOW,RHOMAX,ISAND,              &
                         ILG,IG,IL1,IL2,JL)                             
!                                                                       
!     * APR 17/14 - D.VERSEGHY. MAKE SNOW ALBEDO REFRESHMENT VALUE
!     *                         CONSISTENT WITH SNOADD.
!     * MAR 07/07 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS.           
!     * MAR 24/06 - D.VERSEGHY. ALLOW FOR PRESENCE OF WATER IN SNOW.    
!     * SEP 23/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
!     * SEP 10/04 - R.HARVEY/D.VERSEGHY. INCREASE SNOW ALBEDO           
!     *                         REFRESHMENT AND WETTING THRESHOLDS.     
!     * AUG 04/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE              
!     *                         CALCULATIONS AGAIST ROUNDOFF            
!     *                         ERRORS.                                 
!     * APR 21/04 - F.SEGLENIEKS/D.VERSEGHY. BUG FIX IN SNOW            
!     *                         TEMPERATURE COMPARISONS.                
!     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.          
!     * OCT 20/00 - R.BROWN/D.VERSEGHY. MODIFIED SNOW DENSITY           
!     *                                 CALCULATIONS, ACCOUNTING        
!     *                                 FOR SETTLING IN WARM AND        
!     *                                 COLD SNOW.                      
!     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
!     *                         SPECIFY LOCATION OF ICE SHEETS          
!     *                         BY SOIL TEXTURE ARRAY RATHER            
!     *                         THAN BY SOIL COLOUR INDEX.              
!     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.                    
!     *                         COMPLETION OF ENERGY BALANCE            
!     *                         DIAGNOSTICS.                            
!     * MAR 13/92 - M.LAZARE.   CLASS - VERSION 2.1.                    
!     *                         CODE FOR MODEL VERSION GCM7 -           
!     *                         DIVIDE PREVIOUS SUBROUTINE              
!     *                         "SNOALB" INTO "SNOALBA" AND             
!     *                         "SNOALBW" AND VECTORIZE.                
!     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
!     *                         CLASS VERSION 2.0 (WITH CANOPY).        
!     * APR 11/89 - D.VERSEGHY. CALCULATE DECREASE IN SNOW ALBEDO       
!     *                         AND INCREASE IN DENSITY DUE TO          
!     *                         AGING. (ASSIGN DIFFERENT LOWER          
!     *                         SNOW ALBEDO LIMITS FOR DRY AND          
!     *                         MELTING SNOW.)                          
!                                                                       
      IMPLICIT NONE                                                     
!                                                                       
!     * INTEGER CONSTANTS.                                              
!                                                                       
      INTEGER ILG,IG,IL1,IL2,JL,I,IPTBAD                                
!                                                                       
!     * OUTPUT ARRAYS.                                                  
!                                                                       
      REAL ALBSNO(ILG),   RHOSNO(ILG),   ZSNOW (ILG),   HCPSNO(ILG)     
!                                                                       
!     * INPUT ARRAYS.                                                   
!                                                                       
      REAL TSNOW (ILG),   FI    (ILG),   S     (ILG),   RMELT (ILG),   &
           WSNOW (ILG)                                                  
!                                                                       
      INTEGER             ISAND (ILG,IG)                                
!                                                                       
!     * WORK ARRAY.                                                     
!                                                                       
      REAL RHOMAX(ILG)                                                  
!                                                                       
!     * TEMPORARY VARIABLES.                                            
!                                                                       
      REAL TIMFAC,RHOOLD                                                
!                                                                       
!     * COMMON BLOCK PARAMETERS.                                        
!                                                                       
      REAL DELT,TFREZ,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,     &
           SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP        
!                                                                       
      COMMON /CLASS1/ DELT,TFREZ                                        
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,          &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,           &
                      TCGLAC,CLHMLT,CLHVAP                              
!---------------------------------------------------------------------- 
      IPTBAD=0                                                          
      DO 100 I=IL1,IL2                                                  
          IF(ZSNOW(I).GT.0. .AND. FI(I).GT.0. .AND. S(I)*DELT.LT.1.0E-4)  THEN   
            IF(ALBSNO(I).GT.0.5001 .AND. (RMELT(I).GT.1.0E-7 .OR. TSNOW(I).GE.-0.01)) THEN
                  ALBSNO(I)=(ALBSNO(I)-0.50)*EXP(-0.01*DELT/3600.0)+ 0.50                                              
            ELSE IF(ALBSNO(I).GT.0.7001 .AND. RMELT(I).LE.1.0E-7)   THEN  
                  ALBSNO(I)=(ALBSNO(I)-0.70)*EXP(-0.01*DELT/3600.0)+ 0.70                                              
            ENDIF                                                     
          ENDIF                                                         
!                                                                       
          IF(FI(I).GT.0. .AND. ZSNOW(I).GT.0.0001) THEN  
              IF(TSNOW(I).LT.-0.01)     THEN              
                  RHOMAX(I)=450.0-(204.7/ZSNOW(I))*(1.0-EXP(-ZSNOW(I)/0.673))                        
              ELSE                                                      
                  RHOMAX(I)=700.0-(204.7/ZSNOW(I))*(1.0-EXP(-ZSNOW(I)/0.673))                        
              ENDIF                                                     
          ENDIF                                                         
!                                                                       
          IF(FI(I).GT.0. .AND. ZSNOW(I).GT.0.0001 .AND. RHOSNO(I).LT.(RHOMAX(I)-0.01)) THEN  
              RHOOLD=RHOSNO(I)                                          
              RHOSNO(I)=(RHOSNO(I)-RHOMAX(I))*EXP(-0.01*DELT/3600.0)+ RHOMAX(I)
              ZSNOW(I)=ZSNOW(I)*RHOOLD/RHOSNO(I)                        
              HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/(RHOW*ZSNOW(I))
          ENDIF                                                      
          IF((ALBSNO(I).LT.0.49 .OR. ALBSNO(I).GT.1.0) .AND. ZSNOW (I).GT.0. .AND. FI(I).GT.0.) IPTBAD=I
  100 CONTINUE                                                          
!                                                                       
      IF(IPTBAD.NE.0) THEN                                              
         WRITE(6,6100) IPTBAD,JL,ALBSNO(IPTBAD)                         
 6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALBSNO = ',F10.5)          
         CALL XIT('SNOALBW',-1)                                         
      ENDIF                                                             
!                                                                       
      RETURN                                                            
      END                                                               
