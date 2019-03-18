!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------
!**S/P KUO2
!
      SUBROUTINE KUO5 ( TE,QE,CRR,CSR,ILAB,CCK,CCB,OMEGAP,CLDW, &
                        TP1,TM1,QP1,QM1,GZM1,PSP1,PSM1, &
                        SIGMA, TAU, N, NI, NK, DBGKUO, SATUCO, CCCTIM, CMELT)
      use tdpack
      implicit none
#include <arch_specific.hf>
!
      LOGICAL SATUCO
!
      INTEGER N,NK,NI
      INTEGER ILAB(NI,NK)
      LOGICAL DBGKUO
      REAL TAU
      REAL TP1(N,NK),TM1(N,NK),QP1(N,NK),QM1(N,NK),GZM1(N,NK)
      REAL PSP1(N),PSM1(N),SIGMA(NI,NK)
      REAL TE(NI,NK),QE(NI,NK),CRR(NI),CSR(NI)
      REAL CCK(NI,NK),OMEGAP(N,NK),CLDW(N,NK)
      REAL CCB(NI,NK)
!
!Author
!          J.F.Geleyn E.C.M.W.F.     13/05/82.
!
!Revision
! 001      J. Mailhot   R.P.N. 08/02/85. Adaptation and
!          modification for EFR model.
! 002      J. Mailhot (Mar 1987) base of condensation
! 003      G.Pellerin (Oct87)Adaptation to revised code
! 004      J. Mailhot (Mar 1988) threshold of evaporation
! 005      G.Pellerin(August90)Adaptation to thermo functions
! 006      C. Girard(Nov 90)
!          Remove calculations related to the boundary layer
!          Introduction of the entrainment for the cloud profile
!          Modification of the heating/moistening partition:NFRAC**3
! 007      N. Brunet  (May91)
!          New version of thermodynamic functions
!          and file of constants
! 008      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 009      J. Mailhot  (Dec 1992) - Bug correction to the
!          Newton Method
! 010      C. Girard (Nov92) - Clean-up, and implicit
!          calculation of evaporation and precipitation.
!          Recycling of evaporated precipitation for
!          modelling the increase in humidity and calculation
!          of cloud fraction
! 011      A. Methot (Sept 93) -
!          -Bug correction : L/Cp missing in precip. evap.
!          -Evap. set to zero : logical variable EVAP skips unnecessary
!           calculations when no evaporation of precip.
!          -Change heating/moistening partition :
!           BETA=((1-RH)/(1-RHC))**3;(1-RHC)**-3=4;;RHC~.37
!          -SATUCO set to .FALSE. locally by the use of SATUC.
!           Change the value of CHLS depending on SATUC.
! 012      G. Lemay (Oct 93) - Dynamic memory allocation with stkmemw
! 013      R. Benoit (Dec 93) - Restore the ILAB output for use by
!          Sundqvist scheme. Also use icvmgt to handle integer ILAB.
! 014      A. Methot (Dec 93) - SATUCO is back as before revision 011
!             but the value of CHLS is determined by SATUCO
!             Also ZTSATCO determine whether CHLS or CHLC is used
!             -criteria on vertical velocity (OMEGAP)at SIGMAK1 SIGMAK2
! 015      B. Bilodeau (Feb 94) - Cleanup - Change name from KUO to KUO2
! 018      R. Sarrazin (June 95) evap, beta, liquid water
! 016      B. Bilodeau (June 94) - New physics interface
! 017      G. Pellerin (Jan 94) - Omitted bugfix for ZTSATCO
! 018      R. Sarrazin (Summer 95) - Corrections; add diagnostic cloud water.
! 019      B. Bilodeau (Jan 2001) - Automatic arrays
! 020      M. Lepine (March 2003) -  CVMG... Replacements
! 021      G. Pellerin (Mai 03) - Conversion IBM
!            - calls to optimized routine MFOQST
! 022      B. Bilodeau (Nov 2004) - Correct CRR bug
! 023      L. Spacek (Oct 2011)   - Add CCB
!
!Object
!          to calculate the tendencies of T and Q by moist convective
!          adjustment according to KUO equations
!
!Arguments
!
!          - Output -
! TE       temperature tendency
! QE       specific humidity tendency
! CRR      rate of liquid precipitation
! CSR      rate of solid precipitation (not available)
! ILAB     label array from Kuo scheme
!
!          - Output -
! CCK      cloud fraction due to the Kuo scheme
! CCB      cloud fraction due to the boundary layer scheme
!
!          - Input -
! OMEGAP   vertical velocity in pressure coordinates
!
!          - Output -
! CLDW     diagnostic cloud water
!
!          - Input/Output -
! TP1      temperature at (T+DT)
! TM1      temperature at (T-DT)
! QP1      specific humidity at (T+DT)
! QM1      specific humidity at (T-DT)
!
!          - Input -
! GZM1     geopotential (N,NK)
! PSP1     surface pressure at (T+DT)
! PSM1     surface pressure at (T-DT)
! SIGMA    sigma levels
! TAU      timestep
! N        dimension NI*NJ
! NI       1st horizontal dimension
! NK       vertical dimension
! DBGKUO   .TRUE. to debug KUO routine
!          .FALSE. to not debug
! SATUCO   .TRUE. if water/ice phase for saturation
!          .FALSE. if water phase only for saturation
!
!Notes
!          The process accounts for the moisture convergence for the
!          formation of cumulus clouds, the humidification of the
!          environment and the formation of precipitation. There is
!          an exact balance between these 3 features.  The process
!          adapts the T values at (T+DT) and Q values at (T-DT) (the
!          difference between T-DT and T+DT will give the moisture
!          convergence) and the proportional intensity to available
!          convergence.  The routine is divided into 5 different
!          steps:
!          1)allocation and position for work space
!          2)preliminary computations
!          3)cloud ascent and flagging
!          4)total moisture convergence and mean beta-parameter
!          5)moistening, condensation and evaporation of rain/snow
!
!
!REFERENCE
!     ECMWF RESEARCH MANUAL (VOLUME 3)
!     PHYSICAL PARAMETRIZATION  CHAPITRE 4
!*
!
      real wc0, wmr, pc1, pc2
!
      parameter (wc0=5.0E-4)
      parameter (wmr=5.0E-4)
      parameter (pc1=300.)
      parameter (pc2=0.5)
!
      LOGICAL LO,EVAP
      INTEGER IERR
      INTEGER NKP1,NKM1
      INTEGER IS,IKS,IKS2
      INTEGER JK,JL
      INTEGER HEAPERR
      REAL TEMP1, TEMP2, TEMP3, temp4, temp5, temp6
!
!
!
!***********************************************************************
!     AUTOMATIC ARRAYS
!***********************************************************************
!
      logical, dimension(ni) :: lo1
      integer, dimension(ni) :: iqcd, km1, km2, kp1, kp2
      real   , dimension(ni) :: zcdp   ,  zcpd   ,  zrfl   ,  zsfl   ,&
                                zcucov ,  zrfln  ,  zsfln  ,  sigmak1,&
                                sigmak2,  kww1   ,  kww2   ,  zqcd
!
      real, dimension(ni,nk) :: zqsatc ,  zldcpe ,  zpp1   ,  zdsg   ,  &
                                zdpp1  ,  zqac   ,  ztp1   ,  zqp1   ,  &
                                zqsate ,  zbeta  ,  zdqtot ,  ztc    ,  &
                                zqc    ,  ztvp1  ,  zdqloc ,  zdtloc ,  &
                                zdttot
!
!***********************************************************************
!
      REAL ZEVAP, ZMELT, ZCCTIM, ZRCYCL, ZEPFLM, ZEPFLS
      REAL ZEPCDP, ZEPCOV, ZTMST, ZCONS2, ZCONS3
      REAL ZCONS5, ZDPM1, ZLDCP0, ZCOR
      REAL ZDP, ZDQK, ZCUPRO, ZCPDL, ZDTK, ZSQFLN, ZNIMP
      REAL ZRFLN0, ZLVDCP, ENTRM
      real cevap, ccctim, cmelt
!
!*    PHYSICAL CONSTANTS.
!     -------- ----------
!
!
      REAL CHLS, DELTA2

!     *NSTAB* REFERS TO THE MINIMUM NUMBER OF STABLE LAYERS BETWEEN THE
!     TOP AND THE CLOUD BASE.
!     IF *EVAP* IS TRUE: THE EVAPORATION OF PRECIPITATION IS ACTIVATED
!                        OTHERWISE THERE IS NO EVAPORATION
!     *ZEVAP* IS A CONSTANT FOR THE EVAPORATION OF PRECIPITATION
!     *ZCCTIM* IS THE CUMULUS LIFE-TIME VALUE TO COMPUTE ITS CLOUD COVER
!     *ZRCYCL* IS THE EVAPORATED PRECIP RECYCLED FRACTION
!
      NKP1=NK+1
      NKM1=NK-1
!
      EVAP=.TRUE.
      IF ( EVAP ) THEN
       ZEVAP  = 6.0E-05
      ELSE
       ZEVAP  = 0.0
      ENDIF
      ZMELT  = CMELT
      ZCCTIM = CCCTIM
      ZRCYCL = 1.0
      ENTRM  = 5.E-6
!
!
!*    SECURITY PARAMETERS.
!     --------------------
!
!         *ZEPFLM* IS A MINIMUM FLUX TO AVOID DIVIDING BY ZERO IN THE IC
!     PROPORTION CALCULATIONS, *ZEPCDP* AVOIDS DIVIDING BY ZERO IN THE
!     ABSENCE OF CLOUD AND *ZEPCOV* IS A MINIMUM CLOUD COVER.
!
      ZEPFLM=1.E-24
      ZEPFLS=SQRT(ZEPFLM)
      ZEPCDP=1.E-12
      ZEPCOV=1.E-12
!
!*    COMPUTATIONAL CONSTANTS.
!     ------------- ----------
!
      ZTMST  = TAU
      DELTA2 = CPV/CPD - 1.

      IF (SATUCO) THEN
       CHLS   = CHLC + CHLF
      ELSE
       CHLS   = CHLC
      ENDIF
!
      ZCONS2 = 1./(ZTMST*GRAV)
      ZCONS3 = ZCCTIM/ZTMST
      ZCONS5 = 1./ZTMST
!
!
!     ------------------------------------------------------------------
!
!*         1.     ALLOCATION AND POSITION FOR WORK SPACE.
!                 ---------- --- -------- --- ---- ------
!
  100 CONTINUE
!
!
!
!
!***
!
!     METHOD.
!     -------
!
!          IN ORDER TO KEEP THE CODE LINEAR THE ROUTINE HAS BEEN DIVIDED
!     INTO THREE PARTS THAT WE SHALL CALL HERE (A) (B) AND (C).
!          IN (A) PRELIMINARY COMPUTATIONS ARE FIRST PERFORMED. THEN A
!     NEARLY ADIABATIC ASCENT IS ATTEMPTED FOR A CLOUD PARCEL STARTING
!     FROM THE LOWEST MODEL LAYER. THIS CLOUD ASCENT IS COMPUTED
!     IN TERMS OF TEMPERATURE AND SPECIFIC HUMIDITY.
!     ENTRAINMENT IS SIMULATED VIA THE ENTRAINMENT PARAMETER.
!     THE LEVELS ARE FLAGGED ACCORDING TO THE FOLLOWING CODE:
!     0 = STABLE,
!     1 = PART OF THE WELL MIXED BOUNDARY LAYER OR DRY UNSTABLE SLAB,
!     2 = MOIST UNSTABLE SLAB (I.E. CLOUD LAYER) AND
!     3 = LIFTING LEVEL IF DIFFERENT FROM THE GROUND.
!          IN (B) THE TOTAL MOISTURE CONVERGENCE FOR EACH SLAB OF
!     NON-0-FLAGS IS STORED INTO ALL THE CORRESPONDING LAYERS IF IT IS
!     POSITIVE AND THE SAME IS DONE FOR THE MEAN OF 1-RELATIVE HUMIDITY
!     FOR EACH CORRESPONDING SLAB OF 2-FLAGS.
!          IN (C) THE ACTUAL MODIFICATIONS OF TEMPERATURE AND SPECIFIC
!     HUMIDITY ARE COMPUTED. FIRST THE ENVIRONMENTAL MOISTENING IS
!     TAKEN PROPORTIONAL TO THE MOISTURE CONVERGENCE WEIGHTED BY "BETA"
!     AND TO THE SATURATION DEFICIT OF SPECIFIC HUMIDITY. AT THE SAME
!     TIME THE MOISTURE CONVERGENCE IS TAKEN AWAY FROM THE HUMIDITY FIELD.
!     SECOND THE FORMATION OF PRECIPITATION IS TAKEN PROPORTIONAL TO THE
!     MOISTURE CONVERGENCE WEIGHTED BY "1 - BETA"  AND TO THE
!     TEMPERATURE DIFFERENCE BETWEEN THE CLOUD AND THE ENVIRONMENT.
!     "BETA" , THE MOISTENING PARAMETER, IS A FUNCTION OF MEAN REL. HUM.
!     A CLOUD-COVER VALUE IS OBTAINED BY COMPARING THE TIME AT WHICH THE
!     ENVIRONMENT WOULD REACH EQUILIBRIUM WITH THE CLOUD TO A PRESCRIBED
!     LIFE-TIME VALUE FOR THE CLOUD ITSELF.
!
!
!     ------------------------------------------------------------------
!
!*         2.     PRELIMINARY COMPUTATIONS.
!                 ----------- -------------
!
!         Must initialize CRR to zero because it is stored
!         in the permanent bus
          DO JL=1,NI
             CRR(JL) = 0.0
          END DO
!
  200 CONTINUE
!
!*         2.1     NECESSARY SETTING FOR THE CASE THERE WOULD NOT BE
!*                 ANY CONVECTIVE POINT AROUND THE LATITUDE CIRCLE.
!
  210 CONTINUE
!
!*         2.2     MOISTURE ACCESSION, T+1 T,Q VARIABLES AND SATURATION
!*                 MIXING RATIO PLUS INITIAL VALUE FOR TEST FLAG, ETC.
!

!     The moisture accession is set to zero for all levels when the
!     vertical velocity in pressure coordinates OMEGAP is positive
!     (downward) at both sigma levels SIGMAK1 and SIGMAK2.
!
!     Since SIGMAK1 and 2 do not necessarly coincide with sigma levels
!     of the model, OMEGAP is linearly interpolated at levels SIGMAK1 and 2
!     using weighting factors KWW1 and KWW2. KP1 or KP2 and KM1 or KM2 are
!     the indices of the model's sigma levels from which the interpolation
!     takes place.

!
!     SIGMAK1 and SIGMAK2 are the sigma levels at which OMEGAP is tested

      DO 10 JL = 1,NI
         SIGMAK1(JL) = 0.9
         SIGMAK2(JL) = 0.7
         KM1    (JL) = NK
         KM2    (JL) = NK
10    CONTINUE
!
!     general case
      DO 20 JK = 1,NK
!
         DO 30 JL = 1,NI
!
            IF (SIGMA(JL,JK) .LE. SIGMAK1(JL)) KM1(JL) = JK
            IF (SIGMA(JL,JK) .LE. SIGMAK2(JL)) KM2(JL) = JK
!
            KP1(JL) = KM1(JL) + 1
            KP2(JL) = KM2(JL) + 1
!
            KWW1(JL)=(SIGMAK1(JL)-SIGMA(JL,KM1(JL)))/ &
                     (SIGMA(JL,KP1(JL))-SIGMA(JL,KM1(JL)))
!
            KWW2(JL)=(SIGMAK2(JL)-SIGMA(JL,KM2(JL)))/ &
                     (SIGMA(JL,KP2(JL))-SIGMA(JL,KM2(JL)))
!
30       CONTINUE
!
20    CONTINUE
!
!     loop 40 covers special cases
      DO 40 JL = 1,NI
!
         IF (SIGMA(JL,1).GT.SIGMAK1(JL)) THEN
            SIGMAK1(JL) = SIGMA(JL,1)
            KM1(JL)     = 1
            KP1(JL)     = 1
            KWW1(JL)    = 1.0
         ENDIF

         IF (SIGMA(JL,1).GT.SIGMAK2(JL)) THEN
            SIGMAK2(JL) = SIGMA(JL,1)
            KM2(JL)     = 1
            KP2(JL)     = 1
            KWW2(JL)    = 1.0
         ENDIF

         IF ( KM1(JL)  .EQ. NK )         THEN
            KP1(JL)     = KM1(JL)
            KWW1(JL)    = 1.0
         ENDIF

         IF ( KM2(JL)  .EQ. NK )         THEN
            KP2(JL)     = KM2(JL)
            KWW2(JL)    = 1.0
         ENDIF
!
40    CONTINUE
!
!
      do 220 jl=1,ni
      ZDSG(jl,1)=0.5*(SIGMA(jl,2)-SIGMA(jl,1))
      DO 225 JK=2,NKM1
      ZDSG(jl,JK)=0.5*(SIGMA(jl,JK+1)-SIGMA(jl,JK-1))
  225 CONTINUE
      ZDSG(jl,NK)=0.5*(1.-SIGMA(jl,NKM1))+0.5*(1.-SIGMA(jl,NK))
  220 continue
!
      DO 221 JK=1,NK
      DO 221 JL=1,NI
      ZPP1(JL,JK)=SIGMA(jl,JK)*PSP1(JL)
      ZTP1(JL,JK)=TP1(JL,JK)
      ZQP1(JL,JK)=QP1(JL,JK)
      ZDPP1(JL,JK)=ZDSG(jl,JK)*PSP1(JL)
      ZDPM1       =ZDSG(jl,JK)*PSM1(JL)
      ZQAC(JL,JK)=ZQP1(JL,JK)*ZDPP1(JL,JK)-QM1(JL,JK)*ZDPM1
      LO=(OMEGAP(JL,KP1(JL))*KWW1(JL) + &
          OMEGAP(JL,KM1(JL))*(1-KWW1(JL))).GT.0.
      LO=LO .AND. (OMEGAP(JL,KP2(JL))*KWW2(JL)+ &
                   OMEGAP(JL,KM2(JL))*(1-KWW2(JL))).GT.0.
      if (LO) ZQAC(JL,JK)=0.
      ILAB(JL,JK)=0
      ZBETA(JL,JK)=0.
      ZDQTOT(JL,JK)=-1.
  221 CONTINUE
!
      IF(SATUCO)THEN
       CALL mfoqst3(ZQSATE,ZTP1,ZPP1,NI,NK,NI)
      ELSE
       CALL mfoqsa3(ZQSATE,ZTP1,ZPP1,NI,NK,NI)
      ENDIF
!
!*         2.3     COMPUTATIONS WITHIN THE BOUNDARY LAYER.
!
  230 CONTINUE
!
!*         2.4     SPECIFY TC AND QC AT THE LOWEST LAYER TO START THE
!*                 CLOUD ASCENT. CHECK FOR POSITIVE MOISTURE ACCESSION
!*                 BETWEEN SURFACE AND CLOUD BASE.
!*                 ZQC=0 INDICATES STABLE CONDITIONS.
!
  240 CONTINUE
      DO 241 JL=1,NI
      LO=ZQAC(JL,NK).GT.0.
      ZTC(JL,NK)=ZTP1(JL,NK)
      if (LO) then
         ZQC(JL,NK) = QM1(JL,NK)
      else
         ZQC(JL,NK) = 0.
      endif
      IF (LO) ILAB(JL,NK) = 1
  241 CONTINUE
!
!
!     ------------------------------------------------------------------
!
!*         3.     CLOUD ASCENT AND FLAGGING.
!                 ----- ------ --- ---------
!
  300 CONTINUE
!
!*         3.1     CALCULATE TC AND QC AT NEXT LEVEL BY DRY ADIABATIC
!*                 LIFTING AND CONSIDERING LATENT HEAT RELEASE (THE
!*                 GEOPOTENTIAL DIFFERENCE IS CORRECTED TO TAKE INTO
!*                 ACCOUNT THE VIRTUAL TEMPERATURE DIFFERENCE BETWEEN
!*                 THE ASCENT AND THE ENVIRONMENT). THE CONDENSATION
!*                 CALCULATIONS ARE DONE WITH TWO ITERATIONS.
!
  310 CONTINUE
!
      DO 311 JL=1,NI
      ZCPD(JL)=CPD*(1.+DELTA2*ZQC(JL,NK))
      ZTVP1(JL,NK) = FOTVT(ZTP1(JL,NK),QM1(JL,NK))
  311 CONTINUE
!***
      DO 322 JK=NKM1,1,-1
!***
      DO 312 JL=1,NI
      ZTC(JL,JK)=ZTC(JL,JK+1)+(GZM1(JL,JK+1)-GZM1(JL,JK))* &
         (1./ZCPD(JL)+ENTRM*MAX(0.,ZTC(JL,JK+1)-ZTP1(JL,JK+1)))
      ZQC(JL,JK)=ZQC(JL,JK+1)+(GZM1(JL,JK+1)-GZM1(JL,JK))* &
         (            ENTRM*MAX(0.,ZQC(JL,JK+1)-QM1(JL,JK+1)))
      TEMP1 = ZTP1(JL,JK)
      ZTVP1(JL,JK) = FOTVT(ZTP1(JL,JK),QM1(JL,JK))
      LO=(FOTVT(ZTC(JL,JK),ZQC(JL,JK)).GT.ZTVP1(JL,JK)).AND. &
         (ZQC(JL,JK).NE.0.)
      IF (LO) ILAB(JL,JK) = 1
      LO = ZTC(JL,JK).LT.TRPL .and. SATUCO
!      ZLDCPE = CVMGT(CHLS,CHLC,LO)/ZCPD(JL)
      if (LO) then
         ZLDCPE(jl,jk) = CHLS/ZCPD(JL)
      else
         ZLDCPE(jl,jk) = CHLC/ZCPD(JL)
      endif
  312 CONTINUE

      IF(SATUCO)THEN
      CALL mfoqst3(ZQSATC(1,jk),ZTC(1,jk),ZPP1(1,jk),NI,1,NI)
      DO 313 JL=1,NI
!      ZQSATC(jl,jk)=FOQST(ZTC(JL,JK),ZPP1(JL,JK))
      ZCOR=ZLDCPE(jl,jk)*FODQS(ZQSATC(jl,jk),ZTC(JL,JK))
      ZQCD(jl)=AMAX1(0.,(ZQC(JL,JK)-ZQSATC(jl,jk))/(1.+ZCOR))
  313 CONTINUE
      ELSE
      CALL mfoqsa3(ZQSATC(1,jk),ZTC(1,jk),ZPP1(1,jk),NI,1,NI)
      DO 317 JL=1,NI
!      ZQSATC(jl,jk)=FOQSA(ZTC(JL,JK),ZPP1(JL,JK))
      ZCOR=ZLDCPE(jl,jk)*FODQA(ZQSATC(jl,jk),ZTC(JL,JK))
      ZQCD(jl)=AMAX1(0.,(ZQC(JL,JK)-ZQSATC(jl,jk))/(1.+ZCOR))
  317 CONTINUE
      ENDIF

       DO JL=1,NI
       LO=ZQCD(jl).EQ.0.
       IQCD(JL) = 1
       IF (LO) IQCD(JL) = 0
       ZQC(JL,JK)=ZQC(JL,JK)-ZQCD(jl)
       ZTC(JL,JK)=ZTC(JL,JK)+ZQCD(jl)*ZLDCPE(jl,jk)
       LO1(JL)=.FALSE.
       LO = ZTC(JL,JK).LT.TRPL
       if (LO) then
         ZLDCPE(jl,jk) = CHLS/ZCPD(JL)
       else
         ZLDCPE(jl,jk) = CHLC/ZCPD(JL)
       endif
       ENDDO

      IS=0
      DO 314 JL=1,NI
      IS=IS+IQCD(JL)
  314 CONTINUE
      IF (IS.NE.0) THEN
      IF(SATUCO)THEN
      CALL mfoqst3(ZQSATC(1,jk),ZTC(1,jk),ZPP1(1,jk),NI,1,NI)
      DO 315 JL=1,NI
!      ZQSATC(jl,jk)=FOQST(ZTC(JL,JK),ZPP1(JL,JK))
      ZCOR=ZLDCPE(jl,jk)*FODQS(ZQSATC(jl,jk),ZTC(JL,JK))
      ZQCD(jl)=(ZQC(JL,JK)-ZQSATC(jl,jk))/(1.+ZCOR)
  315 CONTINUE
      ELSE
      CALL mfoqsa3(ZQSATC(1,jk),ZTC(1,jk),ZPP1(1,jk),NI,1,NI)
      DO 318 JL=1,NI
!      ZQSATC(jl,jk)=FOQSA(ZTC(JL,JK),ZPP1(JL,JK))
      ZCOR=ZLDCPE(jl,jk)*FODQA(ZQSATC(jl,jk),ZTC(JL,JK))
      ZQCD(jl)=(ZQC(JL,JK)-ZQSATC(jl,jk))/(1.+ZCOR)
  318 CONTINUE
      ENDIF
       DO JL=1,NI
       LO1(JL)=IQCD(JL).NE.0
       if (.not. LO1(JL)) ZQCD(jl) = 0.
       ZQC(JL,JK)=ZQC(JL,JK)-ZQCD(jl)
       ZTC(JL,JK)=ZTC(JL,JK)+ZQCD(jl)*ZLDCPE(jl,jk)
       ENDDO
      ENDIF

      DO 316 JL=1,NI
      LO=(FOTVT(ZTC(JL,JK),ZQC(JL,JK)).GT.ZTVP1(JL,JK)) &
         .AND.LO1(JL)
      IF (LO) ILAB(JL,JK) = 2
      LO1(JL)=ILAB(JL,JK).EQ.0
      if (LO1(JL)) ZTC(JL,JK) = ZTP1(JL,JK)
      if (LO1(JL)) ZQC(JL,JK) = 0.
  316 CONTINUE
!
!
!*         3.2     IF NOT AT THE TOP CHECK FOR NEW LIFTING LEVEL, I.E.
!*                 MOISTURE CONVERGENCE IN A STABLE LAYER.
!
!***
      IF (JK.NE.1) THEN
!***
      DO 321 JL=1,NI
      LO=LO1(JL).AND.(ZQAC(JL,JK).GT.0.)
      if (LO) ZTC(JL,JK) = ZTP1(JL,JK)
      if (LO) ZQC(JL,JK) = QM1(JL,JK)
      ZCPD(JL)=CPD*(1.+DELTA2*ZQC(JL,JK))
  321 CONTINUE
!***
      ENDIF
  322 CONTINUE
!
!***
!
!*         3.3     ILAB=0 FOR DRY UNSTABLE LAYERS IF NO CLOUD IS ABOVE
!*                 ILAB=3 FOR LIFTING LEVEL IF LAYER ABOVE IS UNSTABLE.
!*                 IKS2 INDICATES THE HIGHEST TOP OF A CLOUD AROUND THE
!*                 LATITUDE CIRCLE (TO AVOID UNNECESSARY COMPUTATIONS
!*                 LATER).
!
  330 CONTINUE
      IKS2=NKP1
      DO 331 JL=1,NI
      LO=ILAB(JL,1).EQ.1
      IF (LO) ILAB(JL,1) = 0
  331 CONTINUE
      IS=0
      DO 332 JL=1,NI
      IS=IS+ILAB(JL,1)
  332 CONTINUE
      IF (IS.NE.0) IKS2=1
      DO 335 JK=2,NK
      DO 333 JL=1,NI
      LO=(ILAB(JL,JK).EQ.1).AND.(ILAB(JL,JK-1).EQ.0)
      IF (LO) ILAB(JL,JK) = 0
  333 CONTINUE
      IF (IKS2.EQ.NKP1) THEN
      IS=0
      DO 334 JL=1,NI
      IS=IS+ILAB(JL,JK)
  334 CONTINUE
      IF (IS.NE.0) IKS2=JK
      ENDIF
  335 CONTINUE
!***
      IF (IKS2.EQ.NKP1) GO TO 600
!***
      DO 337 JK=NK,2,-1
      DO 336 JL=1,NI
      LO=(ILAB(JL,JK).EQ.0).AND.(ILAB(JL,JK-1).NE.0)
      IF (LO) ILAB(JL,JK) = 3
  336 CONTINUE
  337 CONTINUE
!
!
!     ------------------------------------------------------------------
!
!*         4.     TOTAL MOISTURE CONVERGENCE AND MEAN BETA-PARAMETER.
!                 ----- -------- ----------- --- ---- ---------------
!
  400 CONTINUE
!
!*         4.1     CALCULATE TOTAL MOISTURE ACCESSION FOR UNSTABLE
!*                 LAYERS AND THE PARTITION PARAMETER BETA
!*                 AVERAGED OVER CLOUD LAYERS.
!
  410 CONTINUE
      DO 411 JL=1,NI
      LO=ILAB(JL,NK).GT.0
      if (.not. LO) ZQAC(JL,NK) = 0.
      ZQSATE(JL,NK) = AMAX1(ZQSATE(JL,NK),QM1(JL,NK))
      ZBETA(JL,NK)=0.
      ZCDP(JL)=0.
  411 CONTINUE
      DO 413 JK=NKM1,IKS2,-1
      DO 412 JL=1,NI
      LO=ILAB(JL,JK).GT.0
      if (.not. LO) ZQAC(JL,JK) = 0.
      LO=LO.AND.(ILAB(JL,JK).NE.3)
      if (LO) ZQAC(JL,JK) = ZQAC(JL,JK+1)+ZQAC(JL,JK)
      ZQSATE(JL,JK) = AMAX1(ZQSATE(JL,JK),QM1(JL,JK))
      ZBETA(JL,JK) = AMAX1(0.,QM1(JL,JK))/ZQSATE(JL,JK)
      LO=ILAB(JL,JK).EQ.2
      if (.not. LO) ZBETA(JL,JK) = 0.
      if (LO) then
         ZDP = ZDPP1(JL,JK)
      else
         ZDP = 0.
      endif
      LO=LO.AND.(ILAB(JL,JK+1).EQ.2)
      if (LO) then
         ZBETA(JL,JK) = (ZCDP(JL)*ZBETA(JL,JK+1)+ZDP*ZBETA(JL,JK)) &
                        /AMAX1(ZCDP(JL)+ZDP,ZEPCDP)
      endif
      if (LO) then
         ZCDP(JL) = ZCDP(JL)+ZDP
      else
         ZCDP(JL) = ZDP
      endif
  412 CONTINUE
  413 CONTINUE
!
!
!*         4.2     REPLACE THE MOISTURE ACCESSION AT CLOUD LAYERS BY THE
!*                 TOTAL MOISTURE ACCESSION OF THE WHOLE UNSTABLE LAYER
!*                 AND DO THE SAME FOR THE BETA-PARAMETER MEAN VALUE.
!*                 UPDATE IKS2 DURING THE PROCESS.
!
  420 CONTINUE
      IKS2=NKP1
      DO 421 JL=1,NI
      LO=ZQAC(JL,1).GT.0.
      IF (.NOT.LO) ILAB(JL,1) = 0
  421 CONTINUE
      IS=0
      DO 422 JL=1,NI
      IS=IS+ILAB(JL,1)
  422 CONTINUE
      IF (IS.NE.0) IKS2=1
      DO 425 JK=2,NK
      DO 423 JL=1,NI
      LO=(ILAB(JL,JK).EQ.2).AND.(ILAB(JL,JK-1).EQ.2)
      if (LO) ZQAC(JL,JK) = ZQAC(JL,JK-1)
       temp1=(1-ZBETA(JL,JK))*(1-ZBETA(JL,JK))
        temp1=temp1*(1-ZBETA(JL,JK))
      ZBETA(JL,JK) = MIN( 1.0, 6.0*temp1 )
      if (LO) ZBETA(JL,JK) = ZBETA(JL,JK-1)
      LO=(ZQAC(JL,JK).LE.0.).AND.(ILAB(JL,JK).EQ.2)
      IF (LO) ILAB(JL,JK) = 0
      LO=(ILAB(JL,JK).NE.2).AND.(ILAB(JL,JK-1).EQ.0)
      IF (LO) ILAB(JL,JK) = 0
  423 CONTINUE
      IF (IKS2.EQ.NKP1) THEN
      IS=0
      DO 424 JL=1,NI
      IS=IS+ILAB(JL,JK)
  424 CONTINUE
      IF (IS.NE.0) IKS2=JK
      ENDIF
  425 CONTINUE
!
!
!***
      IF (IKS2.EQ.NKP1) GO TO 600
!***
!
!     ------------------------------------------------------------------
!
!*         5.     MOISTENING, CONDENSATION AND EVAPORATION OF RAIN/SNOW.
!                 ----------- ------------ --- ----------- -- ----------
!
  500 CONTINUE
!
!*         5.1     COMPUTE THE TOTAL MOISTURE DEFICIT IN THE CLOUD
!*                 LAYERS.
!
  510 CONTINUE
      DO 511 JL=1,NI
      ZDQLOC(JL,NK)=0.
      ZDQTOT(JL,NK)=0.
  511 CONTINUE
      DO 513 JK=NKM1,IKS2,-1
      DO 512 JL=1,NI
      ZDQLOC(JL,JK) = ZQSATE(JL,JK)-QM1(JL,JK)
      ZDQK = ZDQLOC(JL,JK)*ZDPP1(JL,JK)
      LO=ILAB(JL,JK).EQ.2
      if (LO) then
         ZDQTOT(JL,JK) = ZDQK
      else
         ZDQTOT(JL,JK) = 0.
      endif
      LO=LO.AND.(ILAB(JL,JK+1).EQ.2)
      if (LO) ZDQTOT(JL,JK) = ZDQTOT(JL,JK+1)+ZDQK
  512 CONTINUE
  513 CONTINUE
      IF (IKS2.EQ.1) THEN
      DO 514 JL=1,NI
      LO=ZDQTOT(JL,1).GT.0.
      if (.not. LO) ZDQTOT(JL,1) = -1.
  514 CONTINUE
      ENDIF
      IKS=MAX0(2,IKS2)
      DO 516 JK=IKS,NK
      DO 515 JL=1,NI
      LO=(ILAB(JL,JK).EQ.2).AND.(ILAB(JL,JK-1).EQ.2)
      if (LO) ZDQTOT(JL,JK) = ZDQTOT(JL,JK-1)
      LO=ZDQTOT(JL,JK).GT.0.
      if (.not. LO) ZDQTOT(JL,JK) = -1.
  515 CONTINUE
  516 CONTINUE
!
!
!*         5.2     REMOVE THE MOISTURE ACCESSION IN THE RELEVANT LAYERS
!*                 (BY RETURNING TO THE PREVIOUS TIMESTEP VALUES)
!*                 AND DO THE ENVIRONMENTAL MOISTENING.
!
  520 CONTINUE
      DO 523 JK=IKS2,NK
      DO 522 JL=1,NI
      LO=ILAB(JL,JK).GT.0
      ZDPM1 = ZDSG(jl,JK)*PSM1(JL)
      if (LO) ZQP1(JL,JK) = QM1(JL,JK)*ZDPM1/ZDPP1(JL,JK)
      ZCUPRO = ZDQLOC(JL,JK)*ZBETA(JL,JK)*ZQAC(JL,JK)/ZDQTOT(JL,JK)
      LO=ILAB(JL,JK).EQ.2
      if (.not. LO) ZCUPRO = 0.
      ZQP1(JL,JK) = ZQP1(JL,JK)+ZCUPRO
  522 CONTINUE
  523 CONTINUE
!
!
!*         5.3     COMPUTATION OF THE TOTAL VIRTUAL TEMPERATURE DEFICIT
!*                 IN CLOUD LAYERS.
!
  530 CONTINUE
      DO 531 JL=1,NI
      ZTC(JL,NK)=ZTP1(JL,NK)
      ZQC(JL,NK)=QM1(JL,NK)
      ZDTLOC(JL,NK)=0.
      ZDTTOT(JL,NK)=0.
  531 CONTINUE
      DO 533 JK=NKM1,IKS2,-1
      DO 532 JL=1,NI
      LO = ZTC(JL,JK).LT.TRPL
!      ZCPDL=CPD*(1.+DELTA2*QM1(JL,JK))/(CVMGT(CHLS,CHLC,LO)*(1.+DELTA
!     *      *QM1(JL,JK)))
      if (LO) then
         ZCPDL=CPD*(1.+DELTA2*QM1(JL,JK))/(CHLS*(1.+DELTA*QM1(JL,JK)))
      else
         ZCPDL=CPD*(1.+DELTA2*QM1(JL,JK))/(CHLC*(1.+DELTA*QM1(JL,JK)))
      endif
      LO=ILAB(JL,JK).EQ.2
!      ZTC(JL,JK) = CVMGT(ZTC(JL,JK),ZTP1(JL,JK),LO)
      if (.not. LO) ZTC(JL,JK) = ZTP1(JL,JK)
!      ZQC(JL,JK) = CVMGT(ZQC(JL,JK),QM1(JL,JK),LO)
      if (.not. LO) ZQC(JL,JK) = QM1(JL,JK)
!      ZQSATE(JL,JK) = CVMGT(ZQC(JL,JK),ZQSATE(JL,JK),LO)
      if (LO) ZQSATE(JL,JK) = ZQC(JL,JK)
      ZDTLOC(JL,JK)=(FOTVT(ZTC(JL,JK),ZQC(JL,JK))-ZTVP1(JL,JK)) &
                    *ZCPDL
      ZDTK = ZDTLOC(JL,JK)*ZDPP1(JL,JK)
!      ZDTTOT(JL,JK) = CVMGT(ZDTK,0.,LO)
      if (LO) then
         ZDTTOT(JL,JK) = ZDTK
      else
         ZDTTOT(JL,JK) = 0.
      endif
      LO=LO.AND.(ILAB(JL,JK+1).EQ.2)
      if (LO) ZDTTOT(JL,JK) = ZDTTOT(JL,JK+1)+ZDTK
  532 CONTINUE
  533 CONTINUE
      IF (IKS2.EQ.1) THEN
      DO 534 JL=1,NI
      LO=ZDTTOT(JL,1).GT.0.
      if (.not. LO) ZDTTOT(JL,1) = -1.
  534 CONTINUE
      ENDIF
      DO 536 JK=IKS,NK
      DO 535 JL=1,NI
      LO=(ILAB(JL,JK).EQ.2).AND.(ILAB(JL,JK-1).EQ.2)
      if (LO) ZDTTOT(JL,JK) = ZDTTOT(JL,JK-1)
      LO=ZDTTOT(JL,JK).GT.0.
      if (.not. LO) ZDTTOT(JL,JK) = -1.
  535 CONTINUE
  536 CONTINUE
!
!
!*         5.4     FORMATION OF PRECIPITATIONS.
!
  540 CONTINUE
      DO 541 JL=1,NI
      ZRFL(JL)=0.
      ZSFL(JL)=0.
      ZCUCOV(JL)=ZEPCOV
  541 CONTINUE
!
!***
!
      DO 582 JK=IKS2,NK
!***
      DO 542 JL=1,NI
      LO=ILAB(JL,JK).EQ.2
      if (.not. LO) ZQAC(JL,JK) = 0.
      temp6 = ((1.-ZBETA(JL,JK))*ZQAC(JL,JK)/ZDTTOT(JL,JK))*ZCONS3
      ZCUCOV(JL) = AMIN1(AMAX1(ZCUCOV(JL),temp6),0.5)
      ZCUPRO = AMAX1(ZDTLOC(JL,JK)*((1.-ZBETA(JL,JK))*ZQAC(JL,JK) &
                 /ZDTTOT(JL,JK)),0.)
      ZRFLN(JL)= ZRFL(JL)+(ZCUPRO*ZDPP1(JL,JK)*ZCONS2)
      ZSFLN(JL)= ZSFL(JL)+(ZCUPRO*ZDPP1(JL,JK)*ZCONS2)
!
!
!     diagnostic d'eau liquide
!     ------------------------
!
      if( ilab(jl,jk) .eq. 2 ) then
!
      temp6 = amax1(1.E-12,amin1(temp6,0.5))
      temp4 = ((zrfln(jl)-zrfl(jl)) * (grav/zdpp1(jl,jk))) / &
                   (temp6*wc0*wmr)
      temp1 = temp4
!
      temp5 = temp1*temp1
      temp2 = temp1*(1.0-exp(-temp5)) - temp4
      temp3 = 1.0 + (2.0*temp5 - 1.0)*exp(-temp5)
      temp1 = amax1(temp1 - ( temp2/(temp3+1.E-12) ), 0.0)
!
      temp5 = temp1*temp1
      temp2 = temp1*(1.0-exp(-temp5)) - temp4
      temp3 = 1.0 + (2.0*temp5 - 1.0)*exp(-temp5)
      temp1 = amax1(temp1 - ( temp2/(temp3+1.E-12) ), 0.0)
!
      if(ztc(jl,jk) .lt. 268.) then
        temp2 = wmr / (1.0+pc2*(268.-ztc(jl,jk))**0.5)
      else
        temp2 = wmr
      endif
!
      temp2 = temp2 / (1.0 + pc1*(0.5*(zrfln(jl)+zrfl(jl)))**0.5)
!
      cldw(jl,jk) = temp1 * temp6 * temp2
!
      endif
!
  542 CONTINUE
!
!***
      IF (JK.GT.1 .AND. EVAP ) THEN
!***
      DO 561 JL=1,NI
!
!*         5.5     EVAPORATION OF PRECIPITATIONS.
!                  WITH RECYCLING.
!
      ZSQFLN = SQRT( ZRFLN(JL)/ZCUCOV(JL) )
      TEMP1 = ZQSATE(JL,JK)
      TEMP2 = ZTP1(JL,JK)
       ZCPD(JL)=CPD*(1.+DELTA2*ZQC(JL,JK))
      LO = TEMP2.LT.TRPL .AND. ZTP1(JL,NK).LT.TRPL
      if (LO) then
         ZLDCP0 = CHLS/ZCPD(JL)
      else
         ZLDCP0 = CHLC/ZCPD(JL)
      endif
      ZNIMP = 1. + 2.*(1.+ ZLDCP0*FODQS(TEMP1,TEMP2)) &
                    *ZEVAP*ZSQFLN/ZCONS2
!
      ZRFLN0 = ZRFLN(JL)
!      ZRFLN(JL) = ZCUCOV(JL)*(AMAX1(0.,ZSQFLN-ZEVAP*ZDPP1(JL,JK)
!     *           *AMAX1(0.,ZQSATE(JL,JK)-ZQC(JL,JK))/ZNIMP ))**2
       temp2=AMAX1(0.,ZQSATE(JL,JK)-ZQC(JL,JK))
        temp3=AMAX1(0.,ZSQFLN-ZEVAP*ZDPP1(JL,JK)*temp2/ZNIMP)
       ZRFLN(JL) = ZCUCOV(JL)*temp3*temp3
!
      ZQP1(JL,JK) = ZQP1(JL,JK)-(1.-ZRCYCL)* &
                     (ZRFLN(JL)-ZRFLN0)/(ZDPP1(JL,JK)*ZCONS2)
!
!*         5.6     MELTING/FREEZING OF PRECIPITATIONS.
!                  NO MORE CONSIDERED.
!
  561 CONTINUE
!***
      ENDIF
!***
!
!*         5.7     ADD CONVECTIVE TENDENCIES FOR T AND Q.
!
  570 CONTINUE
!
      DO 571 JL=1,NI
      LO = ZTC(JL,JK).LT.TRPL .AND. ZTP1(JL,NK).LT.TRPL
      if (LO) then
         ZLVDCP = CHLS/ZCPD(JL)
      else
         ZLVDCP = CHLC/ZCPD(JL)
      endif
      TE(JL,JK)=ZLVDCP*(ZRFLN(JL)-ZRFL(JL)) &
               *(GRAV/ZDPP1(JL,JK))
      QE(JL,JK)=(ZQP1(JL,JK)-QP1(JL,JK))*ZCONS5
  571 CONTINUE
!
!*         5.8     SWAP OF FLUXES, END OF VERTICAL LOOP
!
      DO 581 JL=1,NI
      ZRFL(JL)=ZRFLN(JL)
      ZSFL(JL)=ZSFLN(JL)
  581 CONTINUE
!***
  582 CONTINUE
!
!***
!*        5.9    EFFECTS OF RECYCLING,
!*               CONVECTIVE CLOUD FRACTION AND
!*               CONVECTIVE RAIN RATE.
!
  590 CONTINUE
!
      DO 591 JK=IKS2,NK
      DO 591 JL=1,NI
      TE(JL,JK) = TE(JL,JK)+MAX(TE(JL,JK),0.)*ZRCYCL &
                  *(ZSFL(JL)-ZRFL(JL))/MAX(1.E-12,ZSFL(JL))
      LO=ILAB(JL,JK).EQ.2
      if (lo) then
         cck(jl,jk) = zcucov(jl)
      else
         cck(jl,jk) = 0.
      endif
  591 CONTINUE
!
!
      DO 592 JL=1,NI
      CRR(JL)=ZRFL(JL)+(ZSFL(JL)-ZRFL(JL))*ZRCYCL
      CSR(JL)=0.0
  592 CONTINUE
!
!     ------------------------------------------------------------------
!
!*         6.     NECESSARY COMPUTATIONS IF SUBROUTINE IS BY-PASSED.
!                 --------- ------------ -- ---------- -- ----------
!
  600 CONTINUE
!***
!
!     ------------------------------------------------------------------
!
!*         7.     RETURN WORKSPACE.
!                 ------ ----------
!
  700 CONTINUE
!
!
      do jk=1,nk
         do jl=1,ni
            if(ilab(jl,jk).eq.2) then
!               nuages de convection profonde (kuo)
!               f(fbl +ik) = f(fdc+ik)
                ccb(jl,jk) = cck(jl,jk)
             else
!                nuages de convection restreinte
!                f(fdc+ik) = f(fbl+ik)*0.5
                 cck(jl,jk)= ccb(jl,jk)*0.5
             endif
         end do
      end do
!
      RETURN
      END
