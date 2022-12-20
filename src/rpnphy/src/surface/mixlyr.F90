        SUBROUTINE MIXLYR(DTEMP,Q0,NLAK,USTAR,IL1,IL2,ILG, &
                    HDPTH,TKE,DELU,FQU,BFLX,DISS,EXPW,QSTAR,       &
                    FSHEAR,FENTRA,HLAK,LLAK,GRED,TRAN,             &
                    CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,RHOMIX,          &
                    LSTAR,QSENS,QEVAP,LKICEH)
!=======================================================================
!     * FEB 3/12  - M.MACKAY.   SUPPRESS BUOYANCY PRODUCTION UNDER ICE
!     *
!     * MAR 9/09  - M.MACKAY.   COMPUTE BFLX USING RHO INSTEAD OF RHO0 
!     *                         IN HCPW
!     * NOV 15/07 - M.MACKAY.   THERMOCLINE TKE LEAKAGE ADDED
!     * 
!     * OCT  5/07 - M.MACKAY.   LIMIT SET ON MAXIMUM DEEPENING ALLOWED 
!     *                         IN ONE TIMESTEP
!     * OCT  1/07 - M.MACKAY.   BFLX MODIFIED TO INCLUDE SKIN THICKNESS
!     *  
!     * MAR 15/07 - M.MACKAY.   COMPUTES LAKE MIXED LAYER DEPTH, TKE
!     *  
!
      IMPLICIT NONE
!
! ----* GLOBAL LAKE VARIABLES *---------------------------------------
!
      INTEGER,DIMENSION(ILG) :: NLAK
      REAL,DIMENSION(ILG) :: QSTAR,Q0,EXPW,DTEMP,HDPTH,TKE,DELU,DISS, &
                             BFLX,FQU,FSHEAR,FENTRA,HLAK,             &
                             LLAK,GRED,TRAN,RHOMIX,LKICEH
!
! ----* INPUT FIELDS *------------------------------------------------
!
      INTEGER ILG,IL1,IL2
      REAL,DIMENSION(ILG) :: USTAR, CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,  &
                             LSTAR,QSENS,QEVAP
!
! ----* COMMON BLOCKS *------------------------------------------------
!
      REAL DELT,TFREZ,RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN,               &
           HCPW,HCPICE,HCPSOL,                                       &
           HCPOM,HCPSND,HCPCLY,SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,&
           TCGLAC,CLHMLT,CLHVAP,DELTA,CGRAV,CKARM,CPD,AS,ASX,CI,BS,  &
           BETA,FACTN,HMIN
      REAL TKECN,TKECF,TKECE,TKECS,HDPTHMIN,                         &
           TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,           &
           TKECL,DUMAX
      COMMON /CLASS1/ DELT,TFREZ                                       
      COMMON /CLASS2/ RGAS,RGASV,GRAV,SBC,VKC,CT,VMIN
      COMMON /CLASS4/ HCPW,HCPICE,HCPSOL,HCPOM,HCPSND,HCPCLY,        &
                      SPHW,SPHICE,SPHVEG,SPHAIR,RHOW,RHOICE,         &
                      TCGLAC,CLHMLT,CLHVAP
      COMMON /LAKECON/ TKECN,TKECF,TKECE,TKECS,HDPTHMIN,             &
                      TKEMIN,DELMAX,DELMIN,EMSW,DELZLK,DELSKIN,DHMAX,&
                      TKECL,DUMAX
!
! ----* LOCAL VARIABLES *--------------------------------------------
!
      INTEGER I,J
      REAL CN,CF,CE,CS,C1,QH,QINT,QTERM,TKE1,TKE2,DENOM,DPTHMAX,MOL, &
           H2,TI,UMAX,HDPTH_OLD,TKE_OLD,DELU_OLD,G,DEPTH,CL,DELTADU, &
           LTERM
      REAL,DIMENSION(ILG) :: DHDT
!
! ----* LOCAL PARAMETERS *--------------------------------------------
!
      CN=TKECN
      CF=TKECF
      CE=TKECE
      CS=TKECS
      CL=TKECL
       G=GRAV
!
!=======================================================================
      DO 100 I=IL1,IL2
!-----------------------------------------------------------------------
       HDPTH_OLD=HDPTH(I)
       TKE_OLD=TKE(I)
       DELU_OLD=DELU(I)
!
!=======================================================================
! MEAN MIXING LAYER TKE 
!-----------------------------------------------------------------------
! (1) Buoyancy Flux 	
!
       LTERM=LSTAR(I)-QSENS(I)-QEVAP(I)
       DEPTH = HDPTH_OLD + DELSKIN
       QH = QSTAR(I)*( CQ1A(I)*EXP(-CQ1B(I)*DEPTH) + &
                       CQ2A(I)*EXP(-CQ2B(I)*DEPTH) + &
                       CQ3A(I)*EXP(-CQ3B(I)*DEPTH) )
       QINT=(2.*QSTAR(I)/DEPTH)*                                 &
                  ( (CQ1A(I)/CQ1B(I))*(EXP(-CQ1B(I)*DEPTH)-1.)   &
               + (CQ2A(I)/CQ2B(I))*(EXP(-CQ2B(I)*DEPTH)-1.)      &
               + (CQ3A(I)/CQ3B(I))*(EXP(-CQ3B(I)*DEPTH)-1.) )
       QTERM=QSTAR(I)+QH+QINT
       IF (QTERM .LT. 0.)    CALL XIT('MIXLYR',-1)

       BFLX(I)=0.5*DEPTH*(G*EXPW(I)/(SPHW*RHOMIX(I)))*(-LTERM - QTERM)  !m3/s3
!----------
! Suppress buoyancy production under ice
!
       IF (LKICEH(I) .GT. 0.0) THEN
         BFLX(I)=0.0
       ENDIF

!-----------------------------------------------------------------------
! (2) Mechanical Forcing 
!
       FQU(I)= 0.5*CN*CN*CN*USTAR(I)*USTAR(I)*USTAR(I)     !m3/s3

!-----------------------------------------------------------------------
! (3) Dissipation and transport of TKE to thermocline
!       (tendency in TKE due to dissipation and transport)
!
       DISS(I)=0.5*CE*SQRT(TKE_OLD*TKE_OLD*TKE_OLD)         !m3/s3
       TRAN(I)=0.5*CF*SQRT(TKE_OLD*TKE_OLD*TKE_OLD)         !m3/s3

!-----------------------------------------------------------------------
! (4) TKE (m2/s2)
!
!-- Forced tendency
       TKE1= (2.0*DELT/HDPTH_OLD)*( FQU(I)  + BFLX(I) )

!-- Dissipation and transport tendency
!-- Solved analytically using HDPTH_OLD for MLD  
       C1=(CE+CF)/HDPTH_OLD
       TKE2= (1.0/( (0.5*C1*DELT)+(1.0/SQRT(TKE_OLD)) ))  &
           * (1.0/( (0.5*C1*DELT)+(1.0/SQRT(TKE_OLD)) )) - TKE_OLD

       TKE(I)= TKE_OLD + TKE1 + TKE2
!
!=======================================================================
! MIXING LAYER DEPTH (HDPTH) AND TENDENCY (DHDT)
!
       DENOM=TKE_OLD-CS*DELU_OLD*DELU_OLD+EXPW(I)*G*HDPTH_OLD*DTEMP(I)
       IF (ABS(DENOM) .LE. 1.E-10) THEN      
        IF (ABS(DTEMP(I)) .GT. 1.E-10) THEN   
! *** set H to shear penetration depth (Pollard et al 1973) ****
! --need to add correction from Spigel at al 1986
          HDPTH(I)=CS*DELU_OLD*DELU_OLD/(EXPW(I)*G*DTEMP(I)) 
!         print*, "reset to shear penetration depth=========: ",HDPTH(I)
        ELSE
          HDPTH(I)=HDPTH_OLD
        ENDIF
       ELSE
        HDPTH(I)=HDPTH_OLD + DELT*((CF-CL)*SQRT(TKE_OLD*TKE_OLD*TKE_OLD))/DENOM
       ENDIF
!
! *** limit deepening to a maximum of DHMAX in 1 timestep
!
       IF ( (HDPTH(I)-HDPTH_OLD) .GT. DHMAX ) THEN
         HDPTH(I) = HDPTH_OLD + DHMAX
       ENDIF

       DPTHMAX=NLAK(I)*DELZLK
       HDPTH(I)=MIN(MAX(HDPTH(I),HDPTHMIN),DPTHMAX)
       DHDT(I)=(HDPTH(I)-HDPTH_OLD)/DELT

!-----------------------------------------------------------------------
! MIXED LAYER RETREAT FOLLOWING RAYNER (1980)
!
       IF (TKE(I) .LE. TKEMIN) THEN
        TKE(I)=TKEMIN
        IF (-BFLX(I) .GE. 1.0E-15) THEN
         MOL=-HDPTH_OLD*FQU(I)/BFLX(I)         !Monin-Obukhov length
!mdm     HDPTH(I)=MAX(HDPTHMIN,MOL)
         HDPTH(I)=MIN(MAX(MOL,HDPTHMIN),DPTHMAX)
        ELSE
         HDPTH(I)=HDPTHMIN
        ENDIF
        IF (HDPTH(I) .LE. 0.)  CALL XIT('MIXLYR',-2)
        DHDT(I)=0.0
        DELU(I)=0.0
       ENDIF
     
!=======================================================================
! SHEAR PRODUCTION AND ENTRAINMENT FLUXES DIAGNOSTIC
!
       FSHEAR(I)=0.5*CS*DELU_OLD*DELU_OLD*DHDT(I)
       FENTRA(I)=0.5*EXPW(I)*G*HDPTH_OLD*DTEMP(I)*DHDT(I)
!
!=======================================================================
! MEAN MIXING LAYER MOMENTUM (DELU)
!
       DELU(I)=DELU_OLD + DELT*(  (USTAR(I)*USTAR(I)/HDPTH_OLD)  &
                            - (DHDT(I)*DELU_OLD/HDPTH_OLD) )
       DELU(I)=MAX(DELU(I),0.)
       DELTADU=DELU(I)-DELU_OLD
       IF (DELTADU .GT. DUMAX) THEN
         DELU(I)=DELU_OLD + DUMAX
       ENDIF
!
!-----------------------------------------------------------------------
! RESET DELU IF MAXIMUM VALUE REACHED
! (EG Spigel and Imberger, 1980; Spigel et al 1986) 
!
       H2=HLAK(I)-HDPTH(I)
       IF (H2 .GE. 1.0E-12 .AND. GRED(I) .GT. 0.0) THEN
         TI=2.0*LLAK(I)/(SQRT(GRED(I)*HDPTH(I)*H2/HLAK(I)))
!        TI=0.0         !test mdm to turn off shear term
       ELSE
         TI=0.0
       ENDIF
       UMAX=USTAR(I)*USTAR(I)*TI/(4.0*HDPTH(I))
       IF (DELU(I) .GE. UMAX) THEN
         DELU(I)=0.0
       ENDIF
!-----------------------------------------------------------------------
100   CONTINUE

      RETURN
      END
