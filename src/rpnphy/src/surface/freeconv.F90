      SUBROUTINE FREECONV(LKICEH,T0,TLAK,RHOIW,NLAK,NLAKMAX,ILG,IL1,IL2)
!=======================================================================
      IMPLICIT NONE
!
! ----* LAKE MODEL VARIABLES *----------------------------------------
!
      INTEGER NLAKMAX
      INTEGER,DIMENSION(ILG) :: NLAK
      REAL,DIMENSION(ILG) :: T0, LKICEH
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
!
! ----* INPUT *-------------------------------------------
!
      INTEGER ILG,IL1,IL2
      REAL RHOIW
!
! ----* CLASS COMMON BLOCKS *------------------------------------------
!
      REAL DELT,TFREZ
      REAL TKECN,TKECF,TKECE,TKECS,TKECL,HDPTHMIN,                 &
           TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,DUMAX
      COMMON /CLASS1/ DELT,TFREZ                                       
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,           &
                       TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,   &
                       DHMAX,TKECL,DUMAX
!
! ----* LOCAL VARIABLES *---------------------------------------------
!
      INTEGER I,J,K,NMIX
      REAL ZTOP,ZBOT,TTEST,RHO1,RHO2,TC1,TC2,TBAR,XXX,ICEBOT
!=======================================================================
!
      DO 100 I=IL1,IL2
         IF (LKICEH(I) .LE. 0.0) THEN
           TC1=T0(I)-TFREZ
           TC2=TLAK(I,1)-TFREZ
           CALL EQNST(XXX,RHO1,TC1,0.05)
           CALL EQNST(XXX,RHO2,TC2,0.5)
           IF (RHO1 .GT. RHO2) THEN
             TBAR=((DELSKIN*RHO1*T0(I))+(DELZLK*RHO2*TLAK(I,1)))/  &
                   ((DELSKIN*RHO1)+(DELZLK*RHO2))
             T0(I)=TBAR
             TLAK(I,1)=TBAR
           ENDIF
         ENDIF
         ICEBOT=RHOIW*LKICEH(I)
        
        NMIX=1
        DO 420, J=1,NLAK(I)-1
          ZTOP=DELSKIN + (J-1)*DELZLK
          ZBOT=ZTOP+DELZLK
          IF (ICEBOT .LE. ZTOP) THEN
           TC1=TLAK(I,J)-TFREZ
           TC2=TLAK(I,J+1)-TFREZ
!mdm       TTEST=(TC1-3.9816)*(TC2-3.9816)
           TTEST=(TC1-3.98275)*(TC2-3.98275)
           CALL EQNST(XXX,RHO1,TC1,ZBOT)
           CALL EQNST(XXX,RHO2,TC2,ZBOT+DELZLK)
!--------- MIX LAYERS IF RHO1>RHO2 OR TEMPERATURES SPAN 
!--------- T_MAXDENSITY=3.9816 C.
          IF ((RHO1 .GT. RHO2) .OR. (TTEST .LT. 0.0)) THEN
            TBAR=((NMIX*RHO1*TLAK(I,J))+(RHO2*TLAK(I,J+1)))/(NMIX*RHO1+RHO2)
            DO 430, K=J-NMIX+1,J+1
              TLAK(I,K)=TBAR
430         CONTINUE
            NMIX=NMIX+1
          ELSE
            NMIX=1
          ENDIF
         ENDIF
420     CONTINUE
100   CONTINUE
6666  FORMAT(A37,4I5,F5.1)
      RETURN
      END
