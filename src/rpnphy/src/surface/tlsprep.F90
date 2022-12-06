      SUBROUTINE TLSPREP(GCOEFFS,GCONSTS,CPHCHS,TCSNOW,HCPSNO,IWATER,  &
                        ZRSLDM,ZRSLDH,ZRSLFM,ZRSLFH,ZDSLM,ZDSLH,       &
                        ZOSCLM,ZOSCLH,ZOMLNS,ZOELNS,ZOM,ZOH,           &
                        TVIRTA,TPOTA,CRIB,DRAGS,CEVAP,IEVAP,ISAND,     &
                        FLS,ZSNOW,TSNOW,RHOSNO,WSNOW,ZREFM,ZREFH,      &
                        ZDIAGM,ZDIAGH,TA,QA,VA,IZREF,ILG,IL1,IL2,IG,JL)
!
!     * MAY 13/15 - D.VERSEGHY. PREPARATION FOR LAKE SNOW TEMPERATURE
!     *                         CALCULATIONS (BASED ON CLASS 
!     *                         SUBROUTINES CLASST, TPREP AND TSPREP).
!
      IMPLICIT NONE
!                                                                                 
!     * INTEGER CONSTANTS.
!
      INTEGER IZREF,ILG,IL1,IL2,IG,JL,I,J
!
!     * OUTPUT ARRAYS.
!
      REAL GCOEFFS(ILG),   GCONSTS(ILG),   CPHCHS (ILG),   TCSNOW (ILG), &
           ZRSLDM (ILG),   ZRSLDH (ILG),   ZRSLFM (ILG),   ZRSLFH (ILG), &
           ZDSLM  (ILG),   ZDSLH  (ILG),   ZOSCLM (ILG),   ZOSCLH (ILG), &
           ZOMLNS (ILG),   ZOELNS (ILG),   ZOM    (ILG),   ZOH    (ILG), &
           TVIRTA (ILG),   TPOTA  (ILG),   CRIB   (ILG),   DRAGS  (ILG), &
           CEVAP  (ILG),   HCPSNO (ILG)
!
      INTEGER              IWATER(ILG),    IEVAP  (ILG)     
      INTEGER              ISAND (ILG,IG)
!
!     * INPUT ARRAYS.
!
      REAL FLS   (ILG),    ZSNOW (ILG),    TSNOW (ILG),    RHOSNO(ILG),  &
           WSNOW (ILG),    ZREFM (ILG),    ZREFH (ILG),    ZDIAGM(ILG),  &
           ZDIAGH(ILG),    TA    (ILG),    QA    (ILG),    VA    (ILG) 
!
!     * COMMON BLOCK PARAMETERS.
!
      REAL TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,RHOSOL,RHOOM,      &
           HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,    &
           SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP,        &
           DELTA,CGRAV,CKARM,CPD 
!
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,        &
                      RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,     &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,      &
                      TCGLAC,CLHMLT,CLHVAP
      COMMON /PHYCON/ DELTA,CGRAV,CKARM,CPD
!-----------------------------------------------------------------------
!     * CALCULATIONS FOR SNOW-COVERED GROUND.                           
!                                                                       
      DO 100 I=IL1,IL2                                              
          IF(FLS(I).GT.0.)                                      THEN
              ZOM(I)=0.001
              ZOMLNS(I)=LOG(ZOM(I))
              ZOH(I)=0.0003
              ZOELNS(I)=LOG(ZOH(I))
              IF(IZREF.EQ.1) THEN                                   
                  ZRSLDM(I)=ZREFM(I)                                
                  ZRSLDH(I)=ZREFH(I)                                
                  ZRSLFM(I)=ZREFM(I)-ZOM(I)                         
                  ZRSLFH(I)=ZREFH(I)-ZOM(I)                         
                  ZDSLM(I)=ZDIAGM(I)-ZOM(I)                         
                  ZDSLH(I)=ZDIAGH(I)-ZOM(I)                         
                  TPOTA(I)=TA(I)+ZRSLFH(I)*CGRAV/CPD                 
              ELSE                                                  
                  ZRSLDM(I)=ZREFM(I)+ZOM(I)                         
                  ZRSLDH(I)=ZREFH(I)+ZOM(I)                         
                  ZRSLFM(I)=ZREFM(I)                                
                  ZRSLFH(I)=ZREFH(I)                                
                  ZDSLM(I)=ZDIAGM(I)                                
                  ZDSLH(I)=ZDIAGH(I)                                
                  TPOTA(I)=TA(I)                                    
              ENDIF                                                 
              ZOSCLM(I)=ZOM(I)/ZRSLDM(I)                            
              ZOSCLH(I)=ZOH(I)/ZRSLDH(I)                            
              TVIRTA(I)=TPOTA(I)*(1.0+0.61*QA(I))                   
              CRIB(I)=-CGRAV*ZRSLDM(I)/(TVIRTA(I)*VA(I)**2)          
              DRAGS(I)=(CKARM/(LOG(ZRSLDM(I))-ZOMLNS(I)))**2  
          ENDIF                                                     
  100     CONTINUE                                                      
!                                                                       
!     * THERMAL PROPERTIES OF SNOW.                                     
!                                                                       
      DO 200 I=IL1,IL2                                                  
          IF(ZSNOW(I).GT.0.)                                        THEN
              HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/          &
                  (RHOW*ZSNOW(I))                                       
              IF(RHOSNO(I).LT.156.0) THEN                               
                  TCSNOW(I)=0.234E-3*RHOSNO(I)+0.023                    
              ELSE                                                      
                  TCSNOW(I)=3.233E-6*RHOSNO(I)*RHOSNO(I)-1.01E-3*       &
                      RHOSNO(I)+0.138                                   
              ENDIF                                                     
!mdm          TCSNOW(I)=2.576E-6*RHOSNO(I)*RHOSNO(I)+0.074       !test: Mellor
!mdm          TCSNOW(I)=0.021+4.2E-4*RHOSNO(I)+
!mdm 1            2.2E-9*RHOSNO(I)*RHOSNO(I)*RHOSNO(I)      !test: Rogers et al 1995
!mdm          TCSNOW(I)=2.22362*(RHOSNO(I)/1000.0)**1.885   !test: Yen 1981
!mdm          TCSNOW(I)=0.0688*EXP(0.0088*
!mcm 1                      (TSNOW(I)-273.16)+4.6682*RHOSNO(I)/1000.0) !Pitman&Zuckerman 1967
          ENDIF
  200 CONTINUE
!                                                                       
!     * CALCULATE COEFFICIENTS.
!
      DO 300 I=IL1,IL2
          IF(FLS(I).GT.0.)                                          THEN
              GCOEFFS(I)=3.0*TCSNOW(I)/ZSNOW(I)
              GCONSTS(I)=-3.0*TCSNOW(I)*TSNOW(I)/ZSNOW(I)
              CPHCHS(I)=CLHVAP+CLHMLT
              IWATER(I)=2             
          ELSE
              IWATER(I)=1
          ENDIF
          CEVAP(I)=1.0
          IEVAP(I)=1
          DO 250 J=1,IG
              ISAND(I,J)=-4
  250     CONTINUE
  300 CONTINUE
!
      RETURN                                                                      
      END
