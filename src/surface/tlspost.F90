      SUBROUTINE TLSPOST(GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZERO,TSNBOT,  &
                        HTCS,HMFN,QFN,EVAPS,RPCN,TRPCN,SPCN,TSPCN,      &
                        GCONSTS,GCOEFFS,T0,ZSNOW,TCSNOW,HCPSNO,QTRANS,  &
                        RPCP,TRPCP,SPCP,TSPCP,TZEROS,RHOSNI,            &
                        FLS,DELSKIN,ILG,IL1,IL2,JL          )
!
!     * JAN 30/18 - M.LAZARE.   LAST ELEMENT IN CALL, "N", REMOVED
!     *                         BECAUSE NOT IN CALL STATEMENT FROM 
!     *                         ROUTINE CLASSL, AND IS NOT USED.
!     * SEP 01/15 - D.VERSEGHY. LAKE SNOW TEMPERATURE AND HEAT FLUX
!     *                         CALCULATIONS (BASED ON CLASS SUBROUTINE
!     *                         TSPOST); NET SURFACE WATER FLUX TERMS
!     *                         (BASED ON CLASS SUBROUTINE WPREP).
!
      IMPLICIT NONE
!                                                                                 
!     * INTEGER CONSTANTS.
!
      INTEGER ILG,IL1,IL2,JL,I,J
!
!     * OUTPUT ARRAYS.
!
      REAL GZERO (ILG),    TSNBOT(ILG),                                 &
           RPCN  (ILG),    TRPCN (ILG),    SPCN  (ILG),    TSPCN (ILG)
!
!     * INPUT/OUTPUT ARRAYS.
!
      REAL GSNOW (ILG),    TSNOW (ILG),    WSNOW (ILG),    RHOSNO(ILG), &
           QMELTS(ILG),    HTCS  (ILG),    HMFN  (ILG),    QFN   (ILG), &
           EVAPS (ILG)
!
!     * INPUT ARRAYS.
!
      REAL T0    (ILG),    ZSNOW (ILG),    TCSNOW (ILG),   HCPSNO(ILG), &
           QTRANS(ILG),    GCONSTS(ILG),   GCOEFFS(ILG),                &
           RPCP  (ILG),    TRPCP  (ILG),   SPCP   (ILG),   TSPCP (ILG), &
           TZEROS(ILG),    RHOSNI (ILG),   FLS    (ILG)
!
      REAL DELSKIN
!
!     * TEMPORARY VARIABLES.
!
      REAL RADD  (ILG),    SADD  (ILG)
!
      REAL HADD,HCONV,WFREZ
!
!     * COMMON BLOCK PARAMETERS.
!
      REAL DELT,TFREZ,TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,            &
           RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,       &
           SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,TCGLAC,CLHMLT,CLHVAP
!
      COMMON /CLASS1/ DELT,TFREZ                                       
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,            &
                      RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,         &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,          &
                      TCGLAC,CLHMLT,CLHVAP
!-----------------------------------------------------------------------
!
      DO 100 I=IL1,IL2
          IF(FLS(I).GT.0.)                                          THEN
              TSNBOT(I)=(ZSNOW(I)*TSNOW(I)+DELSKIN*T0(I))/                &
                   (ZSNOW(I)+DELSKIN)
              GZERO(I)=-2.0*TCSNOW(I)*(TSNBOT(I)-TSNOW(I))/ZSNOW(I)
!     1                 +TCICE*(T0(I)-TSNBOT(I))/DELSKIN)
              IF(QMELTS(I).LT.0.)                               THEN
                  GSNOW(I)=GSNOW(I)+QMELTS(I)                           
                  QMELTS(I)=0.                                          
              ENDIF                                                     
              TSNOW(I)=TSNOW(I)+(GSNOW(I)-GZERO(I))*DELT/              &
                                (HCPSNO(I)*ZSNOW(I))-TFREZ             
              IF(TSNOW(I).GT.0.)                                THEN
                  QMELTS(I)=QMELTS(I)+TSNOW(I)*HCPSNO(I)*ZSNOW(I)/DELT
                  GSNOW(I)=GSNOW(I)-TSNOW(I)*HCPSNO(I)*ZSNOW(I)/DELT
                  TSNOW(I)=0.                                         
              ENDIF                                                  
!              GZERO(I)=GZERO(I)+QMELTS(I)
!              QMELTS(I)=0.0
          ENDIF
  100 CONTINUE
! 
      DO 200 I=IL1,IL2
           IF(FLS(I).GT.0. .AND. TSNOW(I).LT.0. .AND. WSNOW(I).GT.0.)    &
                                                                    THEN
             HTCS(I)=HTCS(I)-FLS(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/ &
                     DELT
             HADD=-TSNOW(I)*HCPSNO(I)*ZSNOW(I)
             HCONV=CLHMLT*WSNOW(I)
             IF(HADD.LE.HCONV)                           THEN      
                 WFREZ=HADD/CLHMLT
                 HADD=0.0
                 WSNOW(I)=MAX(0.0,WSNOW(I)-WFREZ)
                 TSNOW(I)=0.0
                 RHOSNO(I)=RHOSNO(I)+WFREZ/ZSNOW(I)
                 HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE+HCPW*WSNOW(I)/     &
                     (RHOW*ZSNOW(I))
             ELSE                 
                 HADD=HADD-HCONV 
                 WFREZ=WSNOW(I)
                 WSNOW(I)=0.0
                 RHOSNO(I)=RHOSNO(I)+WFREZ/ZSNOW(I)
                 HCPSNO(I)=HCPICE*RHOSNO(I)/RHOICE
                 TSNOW(I)=-HADD/(HCPSNO(I)*ZSNOW(I))
             ENDIF
             HMFN(I)=HMFN(I)-FLS(I)*CLHMLT*WFREZ/DELT
             HTCS(I)=HTCS(I)-FLS(I)*CLHMLT*WFREZ/DELT
             HTCS(I)=HTCS(I)+FLS(I)*HCPSNO(I)*(TSNOW(I)+TFREZ)*ZSNOW(I)/   &
                     DELT
          ENDIF
  200 CONTINUE
!
      DO 300 I=IL1,IL2
        QFN(I)=FLS(I)*EVAPS(I)*RHOW
        IF(SPCP(I).GT.0. .OR. EVAPS(I).LT.0.) THEN           
            SADD(I)=SPCP(I)-EVAPS(I)*RHOW/RHOSNI(I)
            IF(ABS(SADD(I)).LT.1.0E-12) SADD(I)=0.0
            IF(SADD(I).GT.0.0) THEN                          
                SPCN (I)=SADD(I)                            
                IF(SPCP(I).GT.0.0) THEN
                    TSPCN(I)=TSPCP(I)
                ELSE
                    TSPCN(I)=MIN((TZEROS(I)-TFREZ),0.0)
                ENDIF
                EVAPS(I)=0.0                               
            ELSE                                            
                EVAPS(I)=-SADD(I)*RHOSNI(I)/RHOW           
                SPCN (I)=0.0                               
                TSPCN(I)=0.0                               
            ENDIF                                           
        ELSE                                                
            SPCN (I)=0.0                                   
            TSPCN(I)=0.0                                   
        ENDIF

        IF(RPCP(I).GT.0.)                         THEN     
            RADD(I)=RPCP(I)-EVAPS(I)                      
            IF(ABS(RADD(I)).LT.1.0E-12) RADD(I)=0.0
            IF(RADD(I).GT.0.)   THEN                      
                RPCN (I)=RADD(I)                         
                TRPCN(I)=TRPCP(I)
                EVAPS(I)=0.0                            
            ELSE                                         
                EVAPS(I)=-RADD(I)                       
                RPCN (I)=0.0                            
                TRPCN(I)=0.0                           
            ENDIF                                       
        ELSE                                            
            RPCN (I)=0.0                               
            TRPCN(I)=0.0                               
        ENDIF                                             
  300 CONTINUE
!
      RETURN                    
      END
