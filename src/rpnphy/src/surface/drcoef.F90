      SUBROUTINE DRCOEF(CDM,CDH,RIB,CFLUX,QG,QA,ZOMIN,ZOHIN,   &
                        CRIB,TVIRTG,TVIRTA,VA,FI,ITER,         &
                        ILG,IL1,IL2)
!
!     * NOV 04/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
!     * SEP 10/02 - K.ABDELLA.  BUGFIX IN CALCULATION OF "OLS" (2 PLACES).
!     * MAY 22/02 - N.MCFARLANE.USE THE ASYMPTOTIC VALUE FOR CDH 
!     *                         GENERALLY TO LIMIT FLUXES AND PREVENT
!     *                         OVERSHOOT. (COMMMENTED OUT FOR
!     *                         OFF-LINE RUNS.)
!     * APR 11/01 - M.LAZARE.   SHORTENED "CLASS2" COMMON BLOCK.
!     * OCT 26/99 - E. CHAN.    COMMENT OUT ARTIFICIAL DAMPING OF 
!     *                         TURBULENT FLUXES FOR STAND-ALONE TESTING.
!     * JUL 18/97 - M. LAZARE,  CLASS 2.7. PASS IN ADDITIONAL WORK ARRAYS
!     *             D.VERSEGHY. ZOMIN AND ZOHIN TO USE INTERNALLY,
!     *                         SO THAT INPUT ZOMIN AND ZOHIN DO NOT
!     *                         CHANGE WHEN PASSED BACK TO THE
!     *                         ITERATION (BUG FROM PREVIOUS
!     *                         VERSIONS). PREVIOUS ZOM AND ZOH
!     *                         BECOME WORK ARRAYS. 
!     * MAR 10/97 - M. LAZARE.  PASS IN QG AND QA AND USE TO ONLY
!     *                         DAMP TURBULENT FLUXES USING CDHMOD
!     *                         UNDER STABLE CONDITIONS IF MOISTURE
!     *                         FLUX IS UPWARDS. ALSO, BETTER 
!     *                         DEFINITION OF SURFACE LAYER TOP
!     *                         FROM K.ABDELLA.
!     * MAY 21/96 - K. ABDELLA. MODIFICATION FOR FREE-CONVECTIVE
!     *                         LIMIT ON UNSTABLE SIDE ADDED.
!     * JAN 10/96 - K. ABDELLA. CORRECT ERROR IN AU1 (UNSTABLE
!     *                         SIDE) AND PUT IN PRANDTL NUMBER
!     *                         RANGE (0.74->1)/
!     *                         "CDHMOD" USED ON BOTH STABLE AND
!     *                         UNSTABLE SIDE, TO LIMIT FLUXES
!     *                         OVER LONG TIMESTEP.     
!     * M. LAZARE - FEB 14/95.  USE VARIABLES "ZOLN" AND "ZMLN"
!     *                         ON UNSTABLE SIDE, FOR OPTIMIZATION.
!     *                         THIS IS PREVIOUS VERSION "DRCOEFX".
!     * K. ABDELLA/M. LAZARE. - NOV 30/94.
!
      IMPLICIT NONE

!     * CALCULATES DRAG COEFFICIENTS AND RELATED VARIABLES FOR CLASS.

!     * OUTPUT FIELDS ARE:
!     *    CDM    : STABILITY-DEPENDENT DRAG COEFFICIENT FOR MOMENTUM.
!     *    CDH    : STABILITY-DEPENDENT DRAG COEFFICIENT FOR HEAT.
!     *    RIB    : BULK RICHARDSON NUMBER.
!     *    CFLUX  : CD * MOD(V), BOUNDED BY FREE-CONVECTIVE LIMIT.
!
!     * INPUT FIELDS ARE:
!     *    ZOMIN/: ROUGHNESS HEIGHTS FOR MOMENTUM/HEAT NORMALIZED BY
!     *    ZOHIN   REFERENCE HEIGHT.
!     *    CRIB   : -RGAS*SLTHKEF/(VA**2), WHERE
!     *             SLTHKEF=-LOG(MAX(SGJ(ILEV),SHJ(ILEV)))
!     *    TVIRTG : "SURFACE" VIRTUAL TEMPERATURE.
!     *    TVIRTA : LOWEST LEVEL VIRTUAL TEMPERATURE.
!     *    VA     : AMPLITUDE OF LOWEST LEVEL WIND.
!     *    FI     : FRACTION OF SURFACE TYPE BEING STUDIED.
!     *    QG     : SATURATION SPECIFIC HUMIDITY AT GROUND TEMPERATURE.
!     *    QA     : LOWEST LEVEL SPECIFIC HUMIDITY.
!     *    ITER   : INDEX ARRAY INDICATING IF POINT IS UNDERGOING
!     *             FURTHER ITERATION OR NOT. 

!     *    ZOM/  : WORK ARRAYS USED FOR SCALING ZOMIN/ZOHIN
!     *    ZOH     ON STABLE SIDE, AS PART OF CALCULATION.

!     * INTEGER CONSTANTS.

      INTEGER ILG,IL1,IL2,JL,I

!     * OUTPUT ARRAYS.

      REAL CDM    (ILG),   CDH    (ILG),   RIB    (ILG),   CFLUX  (ILG) 

!     * INPUT ARRAYS.

      REAL ZOMIN  (ILG),   ZOHIN  (ILG),   CRIB   (ILG),   TVIRTG (ILG)
      REAL TVIRTA (ILG),   VA     (ILG),   FI     (ILG)
      REAL QG     (ILG),   QA     (ILG)

      INTEGER ITER(ILG)

!    * WORK ARRAYS.

      REAL ZOM    (ILG),   ZOH    (ILG)

!     * TEMPORARY VARIABLES.

      REAL AA,AA1,BETA,PR,ZLEV,ZS,ZOLN,ZMLN,CPR,ZI,OLSF,OLFACT,  &
           ZL,ZMOL,ZHOL,XM,XH,BH1,BH2,BH,WB,WSTAR,RIB0,WSPEED,   &
           AU1,OLS,PSIM1,PSIM0,PSIH1,PSIH0,USTAR,TSTAR,WTS,AS1,  &
           AS2,AS3,CLIMIT

!     * COMMON BLOCK PARAMETERS.

      REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN

      COMMON /CLASS1/ DELT,TFREZ
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
!-------------------------------------------------------------
      AA=9.5285714
      AA1=14.285714
      BETA=1.2
      PR = 1.
!      
      DO 100 I=IL1,IL2
        IF(FI(I).GT.0. .AND. ITER(I).EQ.1)                         THEN
          RIB(I)=CRIB(I)*(TVIRTG(I)-TVIRTA(I))
          IF(RIB(I).GE.0.0) THEN
            ZLEV=-CRIB(I)*TVIRTA(I)*(VA(I)**2)/GRAV  
            ZS=MAX(10.,5.*MAX(ZOMIN(I)*ZLEV, ZOHIN(I)*ZLEV))
            ZS=ZLEV*(1.+RIB(I))/(1.+(ZLEV/ZS)*RIB(I))
            ZOM(I)=ZOMIN(I)*ZLEV/ZS
            ZOH(I)=ZOHIN(I)*ZLEV/ZS
            RIB(I)=RIB(I)*ZS/ZLEV
          ELSE
            ZOM(I)=ZOMIN(I)
            ZOH(I)=ZOHIN(I)
          ENDIF
          ZOLN=LOG(ZOH(I))
          ZMLN=LOG(ZOM(I))
          IF(RIB(I).LT.0.0) THEN
            CPR=MAX(ZOLN/ZMLN,0.74)
            CPR=MIN(CPR,1.0)
            ZI=1000.0
            OLSF=BETA**3*ZI*VKC**2/ZMLN**3
            OLFACT=1.7*(LOG(1.+ZOM(I)/ZOH(I)))**0.5+0.9
            OLSF=OLSF*OLFACT
            ZL = -CRIB(I)*TVIRTA(I)*(VA(I)**2)/GRAV
            ZMOL=ZOM(I)*ZL/OLSF
            ZHOL=ZOH(I)*ZL/OLSF
            XM=(1.00-15.0*ZMOL)**(0.250)
            XH=(1.00-9.0*ZHOL)**0.25
            BH1=-LOG(-2.41*ZMOL)+LOG(((1.+XM)/2.)**2*(1.+XM**2)/2.)
            BH1=BH1-2.*ATAN(XM)+ATAN(1.)*2.
            BH1=BH1**1.5
            BH2=-LOG(-0.25*ZHOL)+2.*LOG(((1.00+XH**2)/2.00))
            BH=VKC**3.*BETA**1.5/(BH1*(BH2)**1.5)
            WB=SQRT(GRAV*(TVIRTG(I)-TVIRTA(I))*ZI/TVIRTG(I))
            WSTAR=BH**(0.333333)*WB
            RIB0=RIB(I)

            WSPEED=SQRT(VA(I)**2+(BETA*WSTAR)**2)
            RIB(I)=RIB0*VA(I)**2/WSPEED**2
            AU1=1.+5.0*(ZOLN-ZMLN)*RIB(I)*(ZOH(I)/ZOM(I))**0.25
            OLS=-RIB(I)*ZMLN**2/(CPR*ZOLN)*(1.0+AU1/                 &
                 (1.0-RIB(I)/(ZOM(I)*ZOH(I))**0.25))
            PSIM1=LOG(((1.00+(1.00-15.0*OLS)**0.250)/2.00)**2*       &
                  (1.0+(1.00-15.0*OLS)**0.5)/2.0)-2.0*ATAN(          &
                  (1.00-15.0*OLS)**0.250)+ATAN(1.00)*2.00
            PSIM0=LOG(((1.00+(1.00-15.0*OLS*ZOM(I))**0.250)/2.00)**2 &
                  *(1.0+(1.00-15.0*OLS*ZOM(I))**0.5)/2.0)-2.0*       &
                  ATAN((1.00-15.0*OLS*ZOM(I))**0.250)+ATAN(1.00)*2.0
            PSIH1=LOG(((1.00+(1.00-9.0*OLS)**0.50)/2.00)**2)
            PSIH0=LOG(((1.00+(1.00-9.0*OLS*ZOH(I))**0.50)/2.00)**2)

            USTAR=VKC/(-ZMLN-PSIM1+PSIM0)
            TSTAR=VKC/(-ZOLN-PSIH1+PSIH0)
            CDH(I)=USTAR*TSTAR/PR
            WTS=CDH(I)*WSPEED*(TVIRTG(I)-TVIRTA(I))
            WSTAR=(GRAV*ZI/TVIRTG(I)*WTS)**(0.333333)

            WSPEED=SQRT(VA(I)**2+(BETA*WSTAR)**2)
            RIB(I)=RIB0*VA(I)**2/WSPEED**2
            AU1=1.+5.0*(ZOLN-ZMLN)*RIB(I)*(ZOH(I)/ZOM(I))**0.25
            OLS=-RIB(I)*ZMLN**2/(CPR*ZOLN)*(1.0+AU1/                   &
                 (1.0-RIB(I)/(ZOM(I)*ZOH(I))**0.25))
            PSIM1=LOG(((1.00+(1.00-15.0*OLS)**0.250)/2.00)**2*         &
                  (1.0+(1.00-15.0*OLS)**0.5)/2.0)-2.0*ATAN(            &
                  (1.00-15.0*OLS)**0.250)+ATAN(1.00)*2.00
            PSIM0=LOG(((1.00+(1.00-15.0*OLS*ZOM(I))**0.250)/2.00)**2   &
                  *(1.0+(1.00-15.0*OLS*ZOM(I))**0.5)/2.0)-2.0*         &
                  ATAN((1.00-15.0*OLS*ZOM(I))**0.250)+ATAN(1.00)*2.0
            PSIH1=LOG(((1.00+(1.00-9.0*OLS)**0.50)/2.00)**2)
            PSIH0=LOG(((1.00+(1.00-9.0*OLS*ZOH(I))**0.50)/2.00)**2)
 
          ELSE

            WSPEED=VA(I)
            AS1=10.0*ZMLN*(ZOM(I)-1.0)
            AS2=5.00/(2.0-8.53*RIB(I)*EXP(-3.35*RIB(I))+0.05*RIB(I)**2)
!<<<
            AS2=AS2*PR*SQRT(-ZMLN)/2.
            AS3=27./(8.*PR*PR)
!>>>
            OLS=RIB(I)*(ZMLN**2+AS3*AS1*(RIB(I)**2+AS2*RIB(I)))/(AS1*RIB(I)-PR*ZOLN)
            PSIM1=-0.667*(OLS-AA1)*EXP(-0.35*OLS)-AA-OLS
            PSIM0=-0.667*(OLS*ZOM(I)-AA1)*EXP(-0.35*OLS*ZOM(I))-AA-OLS*ZOM(I)
            PSIH1=-(1.0+2.0*OLS/3.0)**1.5-0.667*(OLS-AA1)*EXP(-0.35*OLS)-AA+1.0
            PSIH0=-(1.0+2.0*OLS*ZOH(I)/3.0)**1.5-0.667*(OLS*ZOH(I)-AA1) &
                  *EXP(-0.35*OLS*ZOH(I))-AA+1.0

          ENDIF  
      
          USTAR=VKC/(-ZMLN-PSIM1+PSIM0)
          TSTAR=VKC/(-ZOLN-PSIH1+PSIH0)

          CDM(I)=USTAR**2.0
          CDH(I)=USTAR*TSTAR/PR    
!
!         * CALCULATE CD*MOD(V) UNDER FREE-CONVECTIVE LIMIT.
!
          IF(TVIRTG(I).GT.TVIRTA(I))    THEN
            CLIMIT=1.9E-3*(TVIRTG(I)-TVIRTA(I))**0.333333 
          ELSE
            CLIMIT=0.          
          ENDIF             
          CFLUX(I)=MAX(CDH(I)*WSPEED,CLIMIT)
        ENDIF
 100  CONTINUE  
!-------------------------------------------------------------
      RETURN
      END                                                        
