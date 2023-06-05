      SUBROUTINE SNOALBA(ALVSSN,ALIRSN,ALVSSC,ALIRSC,ALBSNO,           &
                         TRSNOWC, ALSNO, TRSNOWG, FSDB, FSFB, RHOSNO,  &
                         REFSN,BCSN,SNO,CSZ,ZSNOW,FSNOW,ASVDAT,ASIDAT, &
                         ALVSG, ALIRG,                                 &
                         ILG,IG,IL1,IL2,JL,IALS,NBS,ISNOALB)            
!                                                                       
!     * MAR 08/16 - M.Mackay    4-BAND SOLAR RADIATION DISABLED
!     * NOV 16/13 - J.COLE.     Final version for gcm17:                
!     *                         - Fixes to get the proper BC mixing ratio in 
!     *                           snow, which required passing in and using  
!     *                           the snow density RHON.                
!     * JUN 22/13 - J.COLE.     ADD CODE FOR "ISNOALB" OPTION,          
!     *                         WHICH IS BASED ON 4-BAND SOLAR.         
!     * FEB 05/07 - D.VERSEGHY. STREAMLINE CALCULATIONS OF              
!     *                         ALVSSN AND ALIRSN.                      
!     * APR 13/06 - D.VERSEGHY. SEPARATE ALBEDOS FOR OPEN AND           
!     *                         CANOPY-COVERED SNOW.                    
!     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
!     * MAR 18/02 - D.VERSEGHY. UPDATES TO ALLOW ASSIGNMENT OF          
!     *                         USER-SPECIFIED VALUES TO SNOW           
!     *                         ALBEDO.                                 
!     * JUN 05/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
!     *                         SPECIFY LOCATION OF ICE SHEETS          
!     *                         BY SOIL TEXTURE ARRAY RATHER            
!     *                         THAN BY SOIL COLOUR INDEX.              
!     * NOV 29/94 - M.LAZARE.   CLASS - VERSION 2.3.                    
!     *                         CALL ABORT CHANGED TO CALL XIT TO       
!     *                         ENABLE RUNNING ON PC'S.                 
!     * MAR 13/92 - M.LAZARE.   CODE FOR MODEL VERSION GCM7 -           
!     *                         DIVIDE PREVIOUS SUBROUTINE              
!     *                         "SNOALB" INTO "SNOALBA" AND             
!     *                         "SNOALBW" AND VECTORIZE.                
!     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
!     *                         CLASS VERSION 2.0 (WITH CANOPY).        
!     * APR 11/89 - D.VERSEGHY. DISAGGREGATE SNOW ALBEDO INTO           
!     *                         VISIBLE AND NEAR-IR PORTIONS;           
!     *                         CALCULATE TRANSMISSIVITY TO             
!     *                         SHORTWAVE RADIATION.                    
!                                                                       
      IMPLICIT NONE                                                     
!                                                                       
!     * INTEGER CONSTANTS.                                              
!                                                                       
      INTEGER ILG,IG,IL1,IL2,JL,IALS,IPTBAD,I,IB,NBS,ISNOALB            
!                                                                       
!     * OUTPUT ARRAYS.                                                  
!
      REAL   ALSNO(ILG,NBS), TRSNOWG(ILG,NBS)                           
      REAL   ALVSSN(ILG),  ALIRSN(ILG),  ALVSSC(ILG),  ALIRSC(ILG),    &
             ALVSG (ILG),  ALIRG (ILG),  TRSNOWC(ILG)                   
!                                                                       
!     * INPUT ARRAYS.                                                   
!                                                                       
      REAL   FSDB(ILG,NBS), FSFB(ILG,NBS)                               
      REAL   ALBSNO(ILG),  ZSNOW (ILG),  FSNOW (ILG),                  &
             ASVDAT(ILG),  ASIDAT(ILG),  REFSN (ILG),  BCSN  (ILG),    &
             CSZ   (ILG),  SNO   (ILG),  RHOSNO(ILG)                    
!                                                                       
!     * LOCAL ARRAYS                                                    
!                                                                       
      REAL SALBG(ILG,NBS), ALDIR(ILG,NBS), ALDIF(ILG,NBS),             &
                           TRDIR(ILG,NBS), TRDIF(ILG,NBS)               
      REAL REFSNO(ILG), BCSNO(ILG)                                      
      INTEGER C_FLAG(ILG)                                               
!                                                                       
!     * CONSTANTS.                                                      
!                                                                       
      REAL WDIRCT, WDIFF                                                
      INTEGER SUM_C_FLAG                                                
!------------------------------------------------------------------     
      IF (ISNOALB .NE. 0) THEN   !4-BAND SW DISABLED 
         CALL XIT('SNOALBA',-9)                                         
      ENDIF                                                             

      IPTBAD=0                                                          
      DO 100 I=IL1,IL2                                                  
         IF(ALBSNO(I).LT.0.50.AND.ALBSNO(I).GT.0.499) ALBSNO(I)=0.50    
         IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.0)              THEN          
             IF(ALBSNO(I).GT.0.70)                    THEN              
                 ALVSSN(I)=0.79*(ALBSNO(I)-0.70)+0.84                   
                 ALIRSN(I)=1.21*(ALBSNO(I)-0.70)+0.56                   
             ELSE                                                       
                 ALVSSN(I)=0.97*(ALBSNO(I)-0.50)+0.62                   
                 ALIRSN(I)=1.03*(ALBSNO(I)-0.50)+0.38                   
             ENDIF                                                      
             IF(ALVSSN(I).GT.0.999.OR.ALVSSN(I).LT.0.001) IPTBAD=I      
             IF(ALIRSN(I).GT.0.999.OR.ALIRSN(I).LT.0.001) IPTBAD=I      
         ELSE IF(FSNOW(I).GT.0.0 .AND. IALS.EQ.1)         THEN          
             ALVSSN(I)=ASVDAT(I)                                        
             ALIRSN(I)=ASIDAT(I)                                        
         ENDIF                                                          
         ALVSSC(I)=ALVSSN(I)                                            
         ALIRSC(I)=ALIRSN(I)                                            
         TRSNOWC(I)=EXP(-25.0*ZSNOW(I))                                 
  100 CONTINUE                                                          
!                                                                       
      IF(IPTBAD.NE.0) THEN                                              
         WRITE(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)          
 6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSSN,ALIRSN = ',2F10.5)  
         CALL XIT('SNOALBA',-1)                                         
      ENDIF                                                             
!                                                                       
!     IF (ISNOALB .EQ. 0) THEN                                          
         DO I = IL1, IL2                                                
            ALSNO(I,1) = ALVSSN(I)                                      
            ALSNO(I,2) = ALIRSN(I)                                      
            ALSNO(I,3) = ALIRSN(I)                                      
            ALSNO(I,4) = ALIRSN(I)                                      
                                                                        
            TRSNOWG(I,1:NBS) = TRSNOWC(I)                               
         END DO
!     END IF ! ISNOALB                                                  
                                                                        
      RETURN                                                            
      END                                                               
