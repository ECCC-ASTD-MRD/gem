      SUBROUTINE TSOLVL(TLAK,T0,LKICEH,KLAK,GLAK,Q0SAT,            &
                        KSTAR,LSTAR,QSENS,QEVAP,EVAP,ALVS,ALIR,    &
                        QSWIN,QLWIN,CSZ,TA,QA,VA,PRES,RHOAIR,CDH,  &
                        GZEROL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,    &
                        CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,             &
                        CQ1BI,CQ2BI,CQ3BI,G0,                      &
                        NLAKMAX,ILG,IL1,IL2)
!======================================================================
!     * DEC 19/16 - M.LAZARE/   - NSTEP REMOVED (WASN'T USED).
!     *             D.VERSEGHY. - BUGFIX IN ACCOUNTING FOR EFFECT
!     *                           OF FICE ON ALBEDO AND LATENT HEAT.
!     *                         - {ALBW,ALBI} NOW DEFINED IN CLASSL
!     *                           AND PASSED IN, FOR CONSISTENCY.
!     * APR 11/16 - M.MACKAY.   THERMAL EXPANSION OF ICE BUG CORRECTED
!     *                         ICE DRAFT,FREEBOARD NOW INCLUDED
!     *                         SW ATTENUATION THROUGH LEADS INCLUDED
!     * SEP 02/15 - D.VERSEGHY. ADD EFFECTS OF SNOW COVERAGE ON SURFACE
!     *                         FLUXES; ADDITIONAL DIAGNOSTIC VARIABLES; 
!     *                         COSMETIC CHANGES TO CODE AND NAMING.
!     * SEP  2/11 - M.MACKAY.  	ICE COMPUTED IN SKIN LAYER
!     * SEP 28/07 - M.MACKAY.  	SFC ENERGY BALANCE NOW COMPUTED OVER
!     *                         SKIN OF FINITE WIDTH
!     * MAR 15/07 - M.MACKAY.   COMPUTES LAKE SURFACE ENERGY BALANCE
!     *  
!
      IMPLICIT NONE
!
! ----* GLOBAL LAKE VARIABLES *---------------------------------------
!
      INTEGER NLAKMAX
      REAL,DIMENSION(ILG) :: T0,LKICEH,KLAK,GLAK,Q0SAT
      REAL,DIMENSION(ILG,NLAKMAX) :: TLAK
!
! ----* OUTPUT FIELDS *------------------------------------------------
!
      REAL,DIMENSION(ILG) :: KSTAR,LSTAR,QSENS,QEVAP,EVAP,ALVS,ALIR
!
! ----* INPUT FIELDS *------------------------------------------------
!
      REAL,DIMENSION(ILG) :: QSWIN,QLWIN,CSZ,TA,QA,VA,PRES,RHOAIR,CDH, &
                             GZEROL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,   &
                             CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,            &
                             CQ1BI,CQ2BI,CQ3BI
      INTEGER ILG,IL1,IL2
!
! ----* INTERNAL ARRAYS *----------------------------------------------
!
      REAL,DIMENSION(ILG) :: G0,XICE
!
! ----* COMMON BLOCKS *------------------------------------------------
!
      REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,TCW,TCICE,     &
           TCSAND,TCCLAY,TCOM,TCDRYS,RHOSOL,RHOOM,HCPW,HCPICE,HCPSOL,&
           HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,&
           TCGLAC,CLHMLT,CLHVAP,CGRAV,CKARM,CPD,AS,ASX,CI,BS,        &
           BETA,FACTN,HMIN
      REAL TKECN,TKECF,TKECE,TKECS,HDPTHMIN,                         &
           TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,           &
           TKECL,DUMAX
      COMMON /CLASS1/ DELT,TFREZ                                       
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS3/ TCW,TCICE,TCSAND,TCCLAY,TCOM,TCDRYS,           &
                      RHOSOL,RHOOM
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,        &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,         &
                      TCGLAC,CLHMLT,CLHVAP
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,              &
                       TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,&
                       TKECL,DUMAX
!
! ----* LOCAL VARIABLES *--------------------------------------------
!
      INTEGER I,J,ITMAX
      REAL DZ,CA,CB,E0SAT,DS
      REAL ALBTOT,T0old,NEWICE,ESKIN,ECOOL
      REAL EAVAIL,ATTEN1,ATTEN2,ATTEN3,EHEAT,TC,CPHCH
      REAL RHOIW,ICETOP,ICEBOT
!
! ----* LOCAL PARAMETER DEFINITIONS *-------------------------------
!
      DZ=DELZLK
      DS=DELSKIN
      RHOIW=RHOICE/RHOW
!----------------------------------------------------------------------
      DO 100 I=IL1,IL2

!======================================================================
! COMPUTE SURFACE ENERGY BALANCE FOR CURRENT TIMESTEP
!======================================================================
! COMPUTE SW FLUXES --------------------------------------------------
!   KSTAR is net SW at surface of skin
!   KLAK is penetrating SW at base of skin
!
        XICE(I)=MAX(FICE(I)-FLS(I),0.0)
        IF(FICE(I).GT.(FLS(I)+0.001))      THEN
            ALBTOT=(XICE(I)*ALBI(I)+(1.-FICE(I))*ALBW(I))/(1.0-FLS(I))
        ELSE
            ALBTOT=ALBW(I)
        ENDIF

        KSTAR(I)=(1.-FLS(I))*(1.-ALBTOT)*QSWIN(I)+FLS(I)*QTRANSL(I)
        IF ( KSTAR(I) .LT. 0.0)  KSTAR(I)=0.0
        ALVS(I)=ALBTOT
        ALIR(I)=ALBTOT

!--- Attenuation through ice (now includes attenuation through leads)
!---
          ICEBOT=RHOIW*LKICEH(I)                !ICE DRAFT
          ICETOP=LKICEH(I)-ICEBOT               !ICE FREEBOARD
          IF ( LKICEH(I) .LE. 0.0) THEN         !NO ICE 
            ATTEN1=CQ1B(I)*DS
            ATTEN2=CQ2B(I)*DS
            ATTEN3=CQ3B(I)*DS
          ELSE 
             IF (ICEBOT .GT. DS) THEN        !Z inside ice
               ATTEN1=FICE(I)*(CQ1BI(I)*(DS+ICETOP)) + &
                      (1.-FICE(I))*CQ1B(I)*DS
               ATTEN2=FICE(I)*(CQ2BI(I)*(DS+ICETOP)) + &
                      (1.-FICE(I))*CQ2B(I)*DS
               ATTEN3=FICE(I)*(CQ3BI(I)*(DS+ICETOP)) + &
                      (1.-FICE(I))*CQ3B(I)*DS
             ELSE
               ATTEN1=FICE(I)*(CQ1BI(I)*LKICEH(I) + CQ1B(I)*(DS-ICEBOT))  &
                    + (1.-FICE(I))*CQ1B(I)*DS
               ATTEN2=FICE(I)*(CQ2BI(I)*LKICEH(I) + CQ2B(I)*(DS-ICEBOT))  &
                    + (1.-FICE(I))*CQ2B(I)*DS
               ATTEN3=FICE(I)*(CQ3BI(I)*LKICEH(I) + CQ3B(I)*(DS-ICEBOT))  &
                    + (1.-FICE(I))*CQ3B(I)*DS
             ENDIF
          ENDIF
          KLAK(I)=KSTAR(I)*(CQ1A(I)*EXP(-ATTEN1) +  &
                            CQ2A(I)*EXP(-ATTEN2) +  &
                            CQ3A(I)*EXP(-ATTEN3) )
!
! COMPUTE TURBULENT FLUXES -------------------------------------------
!     * CALCULATION OF E0SAT CONSISTENT WITH CLASSI
!     * BUT CONSTANTS DIFFER FROM ROGERS&YAU
!     * Rogers and Yau values
!         CA=17.67
!         CB=29.65
!
          IF(T0(I).GE.TFREZ) THEN                              
              CA=17.269                                       
              CB=35.86                                       
          ELSE                                              
              CA=21.874                                    
              CB=7.66                                     
          ENDIF                                          
      
        IF(FICE(I).GT.(FLS(I)+0.001))      THEN
            CPHCH=(XICE(I)*(CLHMLT+CLHVAP)+(1.-FICE(I))*CLHVAP)/(1.0-FLS(I))
        ELSE
            CPHCH=CLHVAP
        ENDIF
        QSENS(I)=(1.-FLS(I))*RHOAIR(I)*SPHAIR*CDH(I)*VA(I)*(T0(I)-TA(I))
        E0SAT=611.0*EXP(CA*(T0(I)-TFREZ)/(T0(I)-CB))
        Q0SAT(I)=0.622*E0SAT/(PRES(I)-0.378*E0SAT)
        EVAP(I)=(1.-FLS(I))*RHOAIR(I)*CDH(I)*VA(I)*(Q0SAT(I)-QA(I))
!mdm    QEVAP(I)=CPHCH*EVAP(I)
        QEVAP(I)=(1.-FLS(I))*RHOAIR(I)*CPHCH*CDH(I)*VA(I)*(Q0SAT(I)-QA(I))
!
! COMPUTE NET LW AND NET SFC ENERGY -------------------------------
! GLAK IS THERMAL FLUX AT BASE OF SKIN.  THERMAL CONDUCTIVITY BASED
! ON WEIGHTED AVERAGE FOR WATER AND ICE IF ICE PRESENT IN LAYER
!
        LSTAR(I)=(1.-FLS(I))*EMSW*(QLWIN(I)-SBC*T0(I)*T0(I)*T0(I)*T0(I))
        G0(I)=KSTAR(I)-KLAK(I)+LSTAR(I)-QSENS(I)-QEVAP(I)+HTCL(I)+FLS(I)*GZEROL(I)
        IF (ICEBOT .GE. DS) THEN 
          GLAK(I)=(-2.0*TCICE/(DZ+DS))*(TLAK(I,1)-T0(I))
        ELSE IF (ICEBOT .LT. DS .AND. LKICEH(I) .GT. 0.0) THEN
          TC=(ICEBOT*TCICE + (DS-ICEBOT)*TCW)/DS
!mdm      TC=20.0          !mdm test
          GLAK(I)=(-2.0*TC/(DZ+DS))*(TLAK(I,1)-T0(I))
        ELSE
          GLAK(I)=(-2.0*TCW/(DZ+DS))*(TLAK(I,1)-T0(I))
        ENDIF
!-----NET ENERGY FLUX INTO SKIN (W/M2)
        ESKIN= G0(I) - GLAK(I)
!
! STEP FORWARD SKIN TEMP T0
!
        T0old=T0(I)
        IF (LKICEH(I) .LE. 0.0) THEN
          T0(I) = T0(I) + (DELT/(DS*HCPW))*ESKIN
        ELSE IF (LKICEH(I) .GT. 0.0 .AND. ICEBOT .LE. DS) THEN
          T0(I) = T0(I) + (DELT/((LKICEH(I)*HCPICE)+(DS-ICEBOT)*HCPW))*ESKIN
        ELSE
          T0(I) = T0(I) + (DELT/((DS+ICETOP)*HCPICE))*ESKIN
        ENDIF
!
! ICE GROWTH OR DECAY
!
        IF (ESKIN .LT. 0.0 .AND. ICEBOT .LT. DS) THEN
!-----NET ENERGY FLUX USED TO LOWER T0 TO TFREZ 
          IF (T0old .GT. TFREZ) THEN
            ECOOL=(DS-ICEBOT)*HCPW*(T0old-TFREZ)/DELT
          ELSE
            ECOOL=0.0
          ENDIF
!-----REMAINING ENERGY FLUX (IF ANY) USED TO FREEZE ICE
          EAVAIL=ESKIN+ECOOL
          IF (EAVAIL .LT. 0.0) THEN
           NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL
           LKICEH(I)=LKICEH(I)+NEWICE
           ICEBOT=RHOIW*LKICEH(I) 
           T0(I)=TFREZ
!-----LIMIT ICE GROWTH TO THE CURRENT LAYER
           IF (ICEBOT .GT. DS) THEN
             EHEAT=(RHOICE*CLHMLT*(ICEBOT-DS))/DELT
             T0(I)=TFREZ - (EHEAT*DELT)/(DS*HCPICE)
             LKICEH(I)=DS/RHOIW
           ENDIF
          ENDIF
        ENDIF

        IF (ESKIN .GT. 0.0 .AND. LKICEH(I) .GT. 0.0) THEN
!-----NET ENERGY FLUX USED FIRST TO RAISE T TO ZERO
            IF (ICEBOT .LE. DS) THEN 
             EHEAT=LKICEH(I)*HCPICE*(TFREZ-T0old)/DELT
            ELSE
             EHEAT=(DS+ICETOP)*HCPICE*(TFREZ-T0old)/DELT
            ENDIF

!-----NET ENERGY FLUX USED TO MELT ICE
          IF (ESKIN .GT. EHEAT) THEN
           NEWICE=-(DELT/(RHOICE*CLHMLT))*(ESKIN-EHEAT)
           LKICEH(I)=LKICEH(I)+NEWICE
           T0(I)=TFREZ
!-----LIMIT ICE MELT TO THE CURRENT LAYER
           IF (LKICEH(I) .LT. 0.0) THEN
             EHEAT=-(RHOICE*CLHMLT*LKICEH(I))/DELT
             T0(I)=TFREZ + (EHEAT*DELT)/(DS*HCPW)
             LKICEH(I)=0.0
           ENDIF
          ENDIF
        ENDIF
! -----------------
100     CONTINUE

      RETURN
      END
