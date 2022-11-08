      SUBROUTINE TSOLVE(ISNOW,FI,                                      &
                        QSWNET,QLWOUT,QTRANS,QSENS,QEVAP,EVAP,         &
                        TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,     &
                        FTEMP,FVAP,ILMO,UE,H,                          &
                        QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,                &
                        ALVISG,ALNIRG,CRIB,CPHCH,CEVAP,TVIRTA,         &
                        ZOSCLH,ZOSCLM,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,      &
                        GCONST,GCOEFF,TSTART,PCPR,TRSNOWG,FSSB,ALSNO,  &
                        THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPOND,         &
                        IWATER,IEVAP,ITERCT,ISAND,                     &
                        ISLFD,ITG,ILG,IG,IL1,IL2,JL,NBS,ISNOALB,       &
                        TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,               &
                        DCFLXM,CFLUXM,WZERO,TRTOP,A,B,                 &
                        LZZ0,LZZ0T,FM,FH,ITER,NITER,JEVAP,KF)           
!                                                                       
!     * OCT 26/16 - D.VERSEGHY. ADD ZPOND TO CALCULATION OF EVPMAX.
!     * JUL 22/15 - D.VERSEGHY. LIMIT CALCULATED EVAPORATION RATE
!     *                         ACCORDING TO WATER AVAILABILITY.
!     * JAN 09/15 - D.VERSEGHY. FIX TO SUPPRESS EVAPORATION FROM ROCK.
!     * JUN 27/14 - D.VERSEGHY. CHANGE ITERATION LIMIT BACK TO 50 FOR
!     *                         BISECTION SCHEME.
!     * NOV 16/13 - J.COLE/     FINAL VERSION FOR GCM17:                
!     *             M.LAZARE.   - FIX COMPUTATION OF QSWNI OVER SNOW FREE 
!     *                           BARE SOIL for ISNOW=0 and ISNOALB=1 (NEED 
!     *                           TO SUM OVER THE 3 NEAR-IR BANDS).     
!     * JUN 22/13 - J.COLE/     - ADD "ISNOALB" OPTION (4-BAND SOLAR).  
!     *             M.LAZARE.   - MODIFY ABORT CONDITION FOR TOO COLD   
!     *                           TEMPS FROM 173 TO 123, SO WON'T       
!     *                           BLOW UP OVER ANTARCTICA.              
!     * OCT 14/11 - D.VERSEGHY. FOR POST-ITERATION CLEANUP WITH N-R SCHEME,
!     *                         REMOVE CONDITION INVOLVING LAST ITERATION 
!     *                         TEMPERATURE.                             
!     * DEC 07/09 - D.VERSEGHY. RESTORE EVAPORATION WHEN PRECIPITATION   
!     *                         IS OCCURRING.                            
!     * MAR 13/09 - D.VERSEGHY. REPLACE SURFCON COMMON BLOCK WITH CLASSD2;
!     *                         REVISED CALL TO FLXSURFZ.               
!     * JAN 06/09 - D.VERSEGHY/M.LAZARE. SPLIT IF CONDITIONS FRAMING    
!     *                         300 LOOP.                               
!     * FEB 25/08 - D.VERSEGHY. STREAMLINE SOME CALCULATIONS; REMOVE    
!     *                         "ILW" SWITCH; SUPPRESS WATER VAPOUR FLUX  
!     *                         IF PRECIPITATION IS OCCURRING.            
!     * MAY 17/06 - D.VERSEGHY. SUPPRESS EVAPORATION WHEN PONDED WATER    
!     *                         IS FREEZING; ADD IL1 AND IL2 TO CALL TO   
!     *                         FLXSURFZ; REMOVE JL FROM CALL TO DRCOEF.  
!     * APR 13/05 - R.BROWN. ADD WINDLESS TRANFER COEFFICIENT TO QSENS    
!     *                         CALCULATION FOR SNOW PACKS.               
!     * DEC 17/04 - Y.DELAGE/D.VERSEGHY. ADD SWITCH TO USE EITHER SECANT/ 
!     *                         BISECTION OR NEWTON-RAPHSON ITERATION     
!     *                         SCHEME (WITH NUMBER OF ITERATIONS LIMITED 
!     *                         TO FIVE AND CORRECTION FOR REMAINING      
!     *                         RESIDUAL).                              
!     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.            
!     * AUG 06/04 - Y.DELAGE/D.VERSEGHY. PROTECT SENSITIVE CALCULATIONS 
!     *                         FROM ROUNDOFF ERRORS.                   
!     * NOV 07/02 - Y.DELAGE/D.VERSEGHY. NEW CALL TO FLXSURFZ.          
!     * JUL 26/02 - D.VERSEGHY. SHORTENED CLASS4 COMMON BLOCK.          
!     * MAR 28/02 - D.VERSEGHY. STREAMLINED SUBROUTINE CALL.            
!     *                         BYPASS EVAPORATION EFFICIENCY PARAMETER 
!     *                         IN CASES OF CONDENSATION.               
!     * JAN 18/02 - P.BARTLETT/D.VERSEGHY. NEW "BETA" FORMULATION FOR   
!     *                         BARE SOIL EVAPORATION BASED ON LEE AND  
!     *                         PIELKE.                                 
!     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.        
!     * OCT 06/00 - D.VERSEGHY. CONDITIONAL "IF" IN ITERATION SEQUENCE  
!     *                         TO AVOID DIVIDE BY ZERO.                
!     * DEC 07/99 - A.WU/D.VERSEGHY. NEW SOIL EVAPORATION FORMULATION.  
!     * JUL 24/97 - D.VERSEGHY. CLASS - VERSION 2.7.                    
!     *                         REPLACE BISECTION METHOD IN SURFACE     
!     *                         TEMPERATURE ITERATION SCHEME WITH       
!     *                         SECANT METHOD FOR FIRST TEN ITERATIONS.   
!     *                         PASS QZERO,QA,ZOMS,ZOHS TO REVISED        
!     *                         DRCOEF (ZOMS AND ZOHS ALSO NEW WORK ARRAYS
!     *                         PASSED TO THIS ROUTINE).                 
!     * JUN 20/97 - D.VERSEGHY. PASS IN NEW "CLASS4" COMMON BLOCK.      
!     * JAN 02/96 - D.VERSEGHY. CLASS - VERSION 2.5.                    
!     *                         COMPLETION OF ENERGY BALANCE            
!     *                         DIAGNOSTICS.  ALSO, PASS SWITCH "ILW"   
!     *                         THROUGH SUBROUTINE CALL, SPECIFYING     
!     *                         WHETHER QLWIN REPRESENTS INCOMING       
!     *                         (ILW=1) OR NET (ILW=2) LONGWAVE         
!     *                         RADIATION ABOVE THE GROUND.             
!     * NOV 30/94 - M.LAZARE.   CLASS - VERSION 2.3.                    
!     *                         NEW DRAG COEFFICIENT AND RELATED FIELDS,
!     *                         NOW DETERMINED IN ROUTINE "DRCOEF"      
!     *                         "CFLUX" NOW WORK FIELD INSTEAD OF "CLIMIT". 
!     * OCT 04/94 - D.VERSEGHY. CHANGE "CALL ABORT" TO "CALL XIT" TO       
!     *                         ENABLE RUNNING ON PCS.                    
!     * JAN 24/94 - M.LAZARE.   UNFORMATTED I/O COMMENTED OUT IN LOOP 200.
!     * JUL 29/93 - D.VERSEGHY. CLASS - VERSION 2.2.                     
!     *                         REMOVE RE-DEFINITION OF QMELT NEAR END   
!     *                         (SINCE DONE ELSEWHERE ALREADY) AND      
!     *                         REDEFINE QSWNET FOR DIAGNOSTIC PURPOSES 
!     *                         TO INCLUDE TRANSMISSION THROUGH         
!     *                         SNOWPACK.                               
!     * OCT 15/92 - D.VERSEGHY/M.LAZARE. CLASS - VERSION 2.1.           
!     *                                  REVISED AND VECTORIZED CODE    
!     *                                  FOR MODEL VERSION GCM7.        
!     * AUG 12/91 - D.VERSEGHY. CODE FOR MODEL VERSION GCM7U -          
!     *                         CLASS VERSION 2.0 (WITH CANOPY).        
!     * APR 11/89 - D.VERSEGHY. ITERATIVE SURFACE TEMPERATURE           
!     *                         CALCULATIONS FOR SNOW/SOIL.             
!                                                                       
      IMPLICIT NONE                                                     
                                                                        
!     * INTEGER CONSTANTS.                                              
!                                                                       
      INTEGER ISNOW,ISLFD,ITG,ILG,IG,IL1,IL2,JL,I,IB,NBS,ISNOALB        
!                                                                       
      INTEGER NUMIT,NIT,IBAD,ITERMX                                     
!                                                                       
!     * OUTPUT ARRAYS.                                                  
!                                                                       
      REAL QSWNET(ILG),    QLWOUT(ILG),    QTRANS(ILG),    QSENS (ILG), &
           QEVAP (ILG),    EVAP  (ILG),    TZERO (ILG),    QZERO (ILG), &
           GZERO (ILG),    QMELT (ILG),    CDH   (ILG),    CDM   (ILG), &
           RIB   (ILG),    CFLUX (ILG),    FTEMP (ILG),    FVAP  (ILG), &
           ILMO  (ILG),    UE    (ILG),    H     (ILG)                  
!                                                                       
!     * INPUT ARRAYS.                                                   
!                                                                       
      REAL FI    (ILG),    QLWIN (ILG),                                 &
           TPOTA (ILG),    QA    (ILG),    VA    (ILG),    PADRY (ILG), &
           RHOAIR(ILG),    ALVISG(ILG),    ALNIRG(ILG),    CRIB  (ILG), &
           CPHCH (ILG),    CEVAP (ILG),    TVIRTA(ILG),                 &
           ZOSCLH(ILG),    ZOSCLM(ILG),    ZRSLFH(ILG),    ZRSLFM(ILG), &
           ZOH   (ILG),    ZOM   (ILG),    GCONST(ILG),    GCOEFF(ILG), &
           TSTART(ILG),    FCOR  (ILG),    PCPR  (ILG),                 &
           RHOSNO(ILG),    ZSNOW (ILG),    ZPOND (ILG) 
!
      REAL THLIQ (ILG,IG), THLMIN(ILG,IG), DELZW (ILG,IG)
!                                                                       
      INTEGER          IWATER(ILG),        IEVAP (ILG)                  
      INTEGER          ITERCT(ILG,6,50),   ISAND(ILG,IG)                
!                                                                       
!     * BAND-DEPENDANT ARRAYS.                                          
!                                                                       
      REAL TRSNOWG(ILG,NBS), ALSNO(ILG,NBS), FSSB(ILG,NBS),             &
           TRTOP  (ILG,NBS)                                             
!                                                                       
!     * INTERNAL WORK ARRAYS.                                           
!                                                                       
      REAL TSTEP (ILG),    TVIRTS(ILG),    EVBETA(ILG),    Q0SAT (ILG), &
           RESID (ILG),    DCFLXM(ILG),    CFLUXM(ILG),                 &
           A     (ILG),    B     (ILG),                                 &
           LZZ0  (ILG),    LZZ0T (ILG),    FM    (ILG),    FH    (ILG), &
           WZERO (ILG),    EVPMAX(ILG)
!                                                                       
      INTEGER              ITER  (ILG),    NITER (ILG),    JEVAP (ILG), &
                           KF    (ILG)                                  
!                                                                       
!     * TEMPORARY VARIABLES.                                            
!                                                                       
      REAL QSWNV,QSWNI,DCFLUX,DRDT0,TZEROT,QEVAPT,BOWEN,EZERO           
!                                                                       
!     * COMMON BLOCK PARAMETERS.                                        
!                                                                       
      REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,HCPW,HCPICE,      &
           HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,        &
           RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP,DELTA,CGRAV,CKARM,CPD,      &
           AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX                          
!                                                                       
      COMMON /CLASS1/ DELT,TFREZ                                        
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN                   
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,           &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,            &
                      TCGLAC,CLHMLT,CLHVAP                              
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD                             
      COMMON /CLASSD2/ AS,ASX,CI,BS,BETA,FACTN,HMIN,ANGMAX              
!-----------------------------------------------------------------------
!     * INITIALIZATION AND PRE-ITERATION SEQUENCE.                      
!                                                                       
      IF(ITG.LT.2) THEN                                                 
          ITERMX=50                                                     
      ELSE                                                              
          ITERMX=5                                                      
      ENDIF                                                             
!                                                                       
!      IF(ISNOW.EQ.0) THEN                                              
!          EZERO=0.0                                                    
!      ELSE                                                             
!          EZERO=2.0                                                    
!      ENDIF                                                            
       EZERO=0.0                                                        
!                                                                       
      DO I=IL1,IL2                                                      
         QSWNET(I)=0.0                                                  
         QTRANS(I)=0.0                                                  
      END DO                                                            
!                                                                       
      IF(ISNOW.EQ.0)    THEN ! Use usual snow-free bare soil formulation
         DO I=IL1,IL2                                                   
            IF(FI(I).GT.0.) THEN                                        
               TRTOP(I,1)=0.                                            
               QSWNV=FSSB(I,1)*(1.0-ALVISG(I))                          
               IF (ISNOALB .EQ. 0) THEN                                 
                  QSWNI=FSSB(I,2)*(1.0-ALNIRG(I))                       
               ELSE IF (ISNOALB .EQ. 1) THEN                            
                  QSWNI=0.0                                             
                  DO IB = 2, NBS                                        
                     QSWNI=QSWNI+FSSB(I,IB)*(1.0-ALNIRG(I))             
                  END DO ! IB                                           
               ENDIF                                                    
               QSWNET(I)=QSWNV+QSWNI                                    
               QTRANS(I)=QSWNET(I)*TRTOP(I,1)                           
               QSWNET(I)=QSWNET(I)-QTRANS(I)                            
            END IF                                                      
         END DO ! I                                                     
      ELSE                                                              
         IF (ISNOALB .EQ. 0) THEN ! Use the existing snow albedo and transmission 
            DO I=IL1,IL2                                                
               IF(FI(I).GT.0.) THEN                                     
                  TRTOP(I,1)=TRSNOWG(I,1)                               
                  QSWNV=FSSB(I,1)*(1.0-ALSNO(I,1))                      
                  QSWNI=FSSB(I,2)*(1.0-ALSNO(I,2))                      
                  QSWNET(I)=QSWNV+QSWNI                                 
                  QTRANS(I)=QSWNET(I)*TRTOP(I,1)                        
                  QSWNET(I)=QSWNET(I)-QTRANS(I)                         
               END IF                                                   
            END DO ! I                                                  
         ELSE IF(ISNOALB .EQ. 1) THEN ! Use the band-by-band snow albedo and transmission
            DO I=IL1,IL2                                                
               QTRANS(I) = 0.0                                          
               QSWNET(I) = 0.0                                          
            END DO ! I                                                  
            DO IB = 1, NBS                                              
               DO I=IL1,IL2                                             
                  IF(FI(I).GT.0.) THEN                                  
                     TRTOP(I,IB)=TRSNOWG(I,IB)                          
                     QSWNV=FSSB(I,IB)*(1.0-ALSNO(I,IB))                 
                     QSWNET(I)=QSWNET(I)+FSSB(I,IB)*(1.0-ALSNO(I,IB))   
                     QTRANS(I)=QTRANS(I)+QSWNV*TRTOP(I,IB)              
                  END IF                                                
               END DO ! I                                               
            END DO ! IB                                                 
            DO I=IL1,IL2                                                
               IF(FI(I).GT.0.) THEN                                     
                  QSWNET(I)=QSWNET(I)-QTRANS(I)                         
               END IF                                                   
            END DO ! I                                                  
         END IF ! ISNOALB                                               
      END IF ! ISNOW                                                    
!                                                                       
      DO 50 I=IL1,IL2                                                   
          IF(FI(I).GT.0.)                                          THEN 
              TZERO(I)=TSTART(I)                                        
              TSTEP(I)=1.0                                              
              ITER(I)=1                                                 
              NITER(I)=1                                                
!                                                                       
              QMELT(I)=0.0                                              
              RESID(I)=999999.                                          
              DCFLXM(I)=0.0                                             
              CFLUX(I)=0.0                                              
              IF(ISNOW.EQ.1)                      THEN                  
                  KF(I)=3                                               
                  EVPMAX(I)=RHOSNO(I)*ZSNOW(I)/DELT
              ELSE                                                      
                  KF(I)=6                                               
                  EVPMAX(I)=RHOW*(ZPOND(I)+(THLIQ(I,1)-THLMIN(I,1))* &
                            DELZW(I,1))/DELT
              ENDIF                                                     
          ENDIF                                                         
   50 CONTINUE                                                          
!                                                                       
!     * ITERATION SECTION.                                              
!     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS 
!     * ON THE CURRENT LATITUDE CIRCLE(S).                              
!                                                                       
  100 CONTINUE                                                          
!                                                                       
      NUMIT=0                                                           
      NIT=0                                                             
      DO 150 I=IL1,IL2                                                  
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN 
              NIT=NIT+1                                                 
              CFLUXM(I)=CFLUX(I)                                        
              IF(TZERO(I).GE.TFREZ)                        THEN         
                  A(I)=17.269                                           
                  B(I)=35.86                                            
              ELSE                                                      
                  A(I)=21.874                                           
                  B(I)=7.66                                             
              ENDIF                                                     
              WZERO(I)=0.622*611.0*EXP(A(I)*(TZERO(I)-TFREZ)/           &
                    (TZERO(I)-B(I)))/PADRY(I)                           
              Q0SAT(I)=WZERO(I)/(1.0+WZERO(I))                          
              IF(IWATER(I).GT.0)                              THEN      
                  EVBETA(I)=1.0                                         
                  QZERO(I)=Q0SAT(I)                                     
              ELSE                                                      
                  EVBETA(I)=CEVAP(I)                                    
                  QZERO(I)=EVBETA(I)*Q0SAT(I)+(1.0-EVBETA(I))*QA(I)     
                  IF(QZERO(I).GT.QA(I) .AND. IEVAP(I).EQ.0) THEN        
                      EVBETA(I)=0.0                                     
                      QZERO(I)=QA(I)                                    
                  ENDIF                                                 
              ENDIF                                                     
              TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I))                    
          ENDIF                                                         
  150 CONTINUE                                                          
!                                                                       
      IF(NIT.GT.0)                                                  THEN
!                                                                       
!     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND   
!     * OTHER RELATED QUANTITIES.                                       
!                                                                       
        IF(ISLFD.LT.2) THEN                                             
            CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      &
                         CRIB,TVIRTS,TVIRTA,VA,FI,ITER,                 &
                         ILG,IL1,IL2)                                   
        ELSE                                                            
            CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            &
                          UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            &
                          TZERO,QZERO,H,ZOM,ZOH,                        &
                          LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )     
        ENDIF                                                           
!                                                                       
!     * REMAINING CALCULATIONS.                                         
!                                                                       
        DO 175 I=IL1,IL2                                                
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN 
              QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I)         
              IF(TZERO(I).LT.TPOTA(I))                        THEN      
                  QSENS(I)=(RHOAIR(I)*SPHAIR*CFLUX(I)+EZERO)*(TZERO(I)- &
                      TPOTA(I))                                         
              ELSE                                                      
                  QSENS(I)=RHOAIR(I)*SPHAIR*CFLUX(I)*(TZERO(I)-         &
                      TPOTA(I))                                         
              ENDIF                                                     
              EVAP(I)=RHOAIR(I)*CFLUX(I)*(QZERO(I)-QA(I))               
              IF(EVAP(I).GT.EVPMAX(I)) EVAP(I)=EVPMAX(I)
              QEVAP(I)=CPHCH(I)*EVAP(I)                                 
              GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I)                     
              RESID(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-QSENS(I)-QEVAP(I)-  &
                       GZERO(I)                                         
              IF(ABS(RESID(I)).LT.5.0)                       ITER(I)=0  
              IF(ABS(TSTEP(I)).LT. 1.0E-2)                   ITER(I)=0  
              IF(NITER(I).EQ.ITERMX .AND. ITER(I).EQ.1)      ITER(I)=-1 
          ENDIF                                                         
175     CONTINUE                                                        
      ENDIF                                                             
!                                                                       
      IF(ITG.LT.2) THEN                                                 
!                                                                       
!     * OPTION #1: BISECTION ITERATION METHOD.                          
!                                                                       
      IF(NIT.GT.0)                                                  THEN
        DO 180 I=IL1,IL2                                                
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN 
              IF(NITER(I).EQ.1) THEN                                    
                  IF(RESID(I).GT.0.0) THEN                              
                      TZERO(I)=TZERO(I)+1.0                             
                  ELSE                                                  
                      TZERO(I)=TZERO(I)-1.0                             
                  ENDIF                                                 
              ELSE                                                      
                  IF((RESID(I).GT.0. .AND. TSTEP(I).LT.0.) .OR.         &
                      (RESID(I).LT.0. .AND. TSTEP(I).GT.0.))    THEN    
                      TSTEP(I)=-TSTEP(I)/2.0                            
                  ENDIF                                                 
                  TZERO(I)=TZERO(I)+TSTEP(I)                            
              ENDIF                                                     
          ENDIF                                                         
!                                                                       
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN 
              NITER(I)=NITER(I)+1                                       
              NUMIT=NUMIT+1                                             
          ENDIF                                                         
  180   CONTINUE                                                        
      ENDIF                                                             
!                                                                       
!     DO 185 I=IL1,IL2                                                  
!         IF(FI(I).GT.0. .AND. ITER(I).EQ.-1)                      THEN 
!             WRITE(6,6250) I,JL,RESID(I),TZERO(I),RIB(I)               
!6250         FORMAT('0GROUND ITERATION LIMIT',3X,2I3,3(F8.2,E12.4))    
!         ENDIF                                                         
! 185 CONTINUE                                                          
!                                                                       
      IF(NUMIT.GT.0)                                    GO TO 100       
!                                                                       
      ELSE                                                              
!                                                                       
!     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.                     
!                                                                       
      IF(NIT.GT.0)                                                  THEN
        DO 190 I=IL1,IL2                                                
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                      THEN  
              IF(NITER(I).GT.1)                                 THEN    
                  DCFLUX=(CFLUX(I)-CFLUXM(I))/                          &
                      SIGN(MAX(.001,ABS(TSTEP(I))),TSTEP(I))            
                  IF(ABS(TVIRTA(I)-TVIRTS(I)).LT.0.4)                   &
                      DCFLUX=MAX(DCFLUX,0.8*DCFLXM(I))                  
                  DCFLXM(I)=DCFLUX                                      
              ELSE                                                      
                  DCFLUX=0.                                             
              ENDIF                                                     
              DRDT0= -4.0*SBC*TZERO(I)**3                               &
                 -RHOAIR(I)*SPHAIR*(CFLUX(I)+MAX(0.,TZERO(I)-TPOTA(I))  &
                 *DCFLUX) -GCOEFF(I)                                    &
                 +CPHCH(I)*RHOAIR(I)*(CFLUX(I)*WZERO(I)*A(I)            &
                 *EVBETA(I)*(B(I)-TFREZ)/((TZERO(I)-B(I))*              &
                 (1.0+WZERO(I)))**2-(QZERO(I)-QA(I))*DCFLUX)            
              TSTEP(I)=-RESID(I)/DRDT0                                  
              TSTEP(I)=MAX(-10.,MIN(5.,TSTEP(I)))                       
              TZERO(I)=TZERO(I)+TSTEP(I)                                
              NITER(I)=NITER(I)+1                                       
              NUMIT=NUMIT+1                                             
          ENDIF                                                         
  190   CONTINUE                                                        
      ENDIF                                                             
!                                                                       
      IF(NUMIT.GT.0)                                    GO TO 100       
!                                                                       
!     * IF CONVERGENCE HAS NOT BEEN REACHED, CALCULATE TEMPERATURE AND  
!     * FLUXES ASSUMING NEUTRAL STABILITY.                              
!                                                                       
      DO 195 I=IL1,IL2                                                  
          NUMIT=0                                                       
          JEVAP(I)=0                                                    
          IF(FI(I).GT.0. .AND.ITER(I).EQ.-1)                       THEN 
              TZEROT=TVIRTA(I)/(1.0+0.61*QZERO(I))                      
              IF(ABS(RESID(I)).GT.50.) THEN                             
                  TZERO(I)=TZEROT                                       
                  IF(TZERO(I).GE.TFREZ)                        THEN     
                      A(I)=17.269                                       
                      B(I)=35.86                                        
                  ELSE                                                  
                      A(I)=21.874                                       
                      B(I)=7.66                                         
                  ENDIF                                                 
                  WZERO(I)=0.622*611.0*EXP(A(I)*(TZERO(I)-TFREZ)/       &
                      (TZERO(I)-B(I)))/PADRY(I)                         
                  Q0SAT(I)=WZERO(I)/(1.0+WZERO(I))                      
                  QZERO(I)=EVBETA(I)*Q0SAT(I)+(1.0-EVBETA(I))*QA(I)     
                  QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I)     
                  GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I)                 
                  RESID(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-GZERO(I)        
                  IF(RESID(I).GT.0.)                 THEN               
                      QEVAP(I)=RESID(I)                                 
                  ELSE                                                  
                      QEVAP(I)=RESID(I)*0.5                             
                  ENDIF                                                 
                  IF(IEVAP(I).EQ.0) QEVAP(I)=0.0
                  QSENS(I)=RESID(I)-QEVAP(I)                            
                  RESID(I)=0.                                           
                  EVAP(I)=QEVAP(I)/CPHCH(I)                             
                  TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I))                
                  JEVAP(I)=1                                            
                  NUMIT=NUMIT+1                                         
              ENDIF                                                     
          ENDIF                                                         
  195 CONTINUE                                                          
!                                                                       
      IF(NUMIT.GT.0)                   THEN                             
        IF(ISLFD.LT.2) THEN                                             
            CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      &
                         CRIB,TVIRTS,TVIRTA,VA,FI,JEVAP,                &
                         ILG,IL1,IL2)                                   
        ELSE                                                            
            CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            &
                          UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            &
                          TZERO,QZERO,H,ZOM,ZOH,                        &
                          LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,JEVAP,JL )    
        ENDIF                                                           
      ENDIF                                                             
!                                                                       
      ENDIF                                                             
!                                                                       
!     * CHECK FOR BAD ITERATION TEMPERATURES.                           
!                                                                       
      IBAD=0                                                            
      DO 200 I=IL1,IL2                                                  
          IF(FI(I).GT.0. .AND. (TZERO(I).LT.123.16 .OR.                  &
                                 TZERO(I).GT.373.16))               THEN
              IBAD=I                                                    
          ENDIF                                                         
 200  CONTINUE                                                          
!                                                                       
      IF(IBAD.NE.0)                                                 THEN
          WRITE(6,6275) IBAD,JL,TZERO(IBAD),NITER(IBAD),ISNOW           
 6275     FORMAT('0BAD ITERATION TEMPERATURE',3X,2I3,F16.2,2I4)         
          WRITE(6,6280) QSWNET(IBAD),QLWIN(IBAD),QSENS(IBAD),           &
              QEVAP(IBAD),GZERO(IBAD),CFLUX(IBAD),RIB(IBAD)             
 6280     FORMAT(2X,7F12.4)                                             
          CALL XIT('TSOLVE',-1)                                         
      ENDIF                                                             
!                                                                       
!     * POST-ITERATION CLEAN-UP.                                        
!                                                                       
      NIT=0                                                             
      DO 300 I=IL1,IL2                                                  
          IF(FI(I).GT.0.)                                          THEN 
              IF(((IWATER(I).EQ.1 .AND. TZERO(I).LT.TFREZ) .OR.         &
                  (IWATER(I).EQ.2 .AND. TZERO(I).GT.TFREZ)) .OR.        &
                  (ISAND(I,1).EQ.-4 .AND. TZERO(I).GT.TFREZ))   THEN    
                  TZERO(I)=TFREZ                                        
                  WZERO(I)=0.622*611.0/PADRY(I)                         
                  QZERO(I)=WZERO(I)/(1.0+WZERO(I))                      
                  TVIRTS(I)=TZERO(I)*(1.0+0.61*QZERO(I))                
                  ITER(I)=1                                             
                  NIT=NIT+1                                             
              ELSE                                                      
                  ITER(I)=0                                             
              ENDIF                                                     
          ENDIF                                                         
  300 CONTINUE                                                          
!                                                                       
      IF(NIT.GT.0)                                                  THEN
!                                                                       
!       * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND 
!       * OTHER RELATED QUANTITIES.                                     
!                                                                       
        IF(ISLFD.LT.2) THEN                                             
            CALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      &
                         CRIB,TVIRTS,TVIRTA,VA,FI,ITER,                 &
                         ILG,IL1,IL2)                                   
        ELSE                                                            
            CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            &
                          UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            &
                          TZERO,QZERO,H,ZOM,ZOH,                        &
                          LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )     
        ENDIF                                                           
      ENDIF                                                             
!                                                                       
!     * REMAINING CALCULATIONS.                                         
!                                                                       
      DO 350 I=IL1,IL2                                                  
          IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                       THEN 
              QLWOUT(I)=SBC*TZERO(I)*TZERO(I)*TZERO(I)*TZERO(I)         
              IF(TZERO(I).LT.TPOTA(I))                        THEN      
                  QSENS(I)=(RHOAIR(I)*SPHAIR*CFLUX(I)+EZERO)*(TZERO(I)- &
                      TPOTA(I))                                         
              ELSE                                                      
                  QSENS(I)=RHOAIR(I)*SPHAIR*CFLUX(I)*(TZERO(I)-         &
                      TPOTA(I))                                         
              ENDIF                                                     
              EVAP(I)=RHOAIR(I)*CFLUX(I)*(QZERO(I)-QA(I))               
              IF(EVAP(I).GT.EVPMAX(I)) EVAP(I)=EVPMAX(I)
              QEVAP(I)=CPHCH(I)*EVAP(I)                                 
              GZERO(I)=GCOEFF(I)*TZERO(I)+GCONST(I)                     
              QMELT(I)=QSWNET(I)+QLWIN(I)-QLWOUT(I)-QSENS(I)-QEVAP(I)-  &
                       GZERO(I)                                         
              RESID(I)=0.0                                              
              IF(QMELT(I).LT.0.0) THEN                                  
                  QMELT(I)=QMELT(I)+QEVAP(I)                            
                  QEVAP(I)=0.0                                          
                  EVAP(I) =0.0                                          
              ENDIF                                                     
          ENDIF                                                         
!                                                                       
          IF(FI(I).GT.0.)                                 THEN          
              IF(ABS(EVAP(I)).LT.1.0E-8) THEN                           
                  RESID(I)=RESID(I)+QEVAP(I)                            
                  EVAP(I)=0.0                                           
                  QEVAP(I)=0.0                                          
              ENDIF                                                     
              IF((ISNOW.EQ.1 .AND. QMELT(I).LT.0.0) .OR.                &
                  (ISNOW.EQ.0 .AND. QMELT(I).GT.0.0))     THEN          
                  GZERO(I)=GZERO(I)+QMELT(I)                            
                  QMELT(I)=0.0                                          
              ENDIF                                                     
!              QSENS(I)=QSENS(I)+0.5*RESID(I)                           
!              GZERO(I)=GZERO(I)+0.5*RESID(I)                           
              QSENS(I)=QSENS(I)+RESID(I)                                
              QSWNET(I)=QSWNET(I)+QTRANS(I)                             
              EVAP(I)=EVAP(I)/RHOW                                      
              ITERCT(I,KF(I),NITER(I))=ITERCT(I,KF(I),NITER(I))+1       
          ENDIF                                                         
  350 CONTINUE                                                          
!                                                                       
      RETURN                                                            
      END                                                               
