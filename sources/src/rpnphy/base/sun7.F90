!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

subroutine sun7b(FDOWN,FUP,HEAT,VOZO,OZOTOIT, &
     DSIG,DSH,DSC, &
     DZ,PSOL,TM,WV,aSIG, &
     QOF,CNEB,TAUAE,RMUO,ALBS, &
     LMX,KMX,KMXP,NPTS,MM, &
     RJ,RK,UD,REFZ,TR, &
     UM,RL,RUEF,C1I,RMUE,RAY,TRA,TAUAZ,PIZAZ, &
     CGAZ,OMEGAT,CG,TAU,RE1,RE2,TR1,TE2,S,G,R23, &
     REF,RMU,R1,W,TO1,GG,RNEB,RMUZ,KK, &
     reduc, &
     flss, srd, afldsig, aflsig, R0R, &
     flkmx, flkmxp, RADFIX, RADFLTR)
   use tdpack_const
   implicit none
!!!#include <arch_specific.hf>

   external intchamps
   integer I,J,K,L,N
   integer KI,IK,K0,K1,JJ,MM
   integer KP1,LP1,KM1,ILG,LMX,JJP,KKI,N2J,IAE
   integer IABS,KMX,KMXP,KREF,NPTS,KFIN,LIND,IFIN,KKP4
   integer flkmx, flkmxp
   real CONS,CCAR,RE,WH2O
   real DSHH,DSCC,XNU,VA,XMUE
   real PT,XMPT,BPT,ZC,RE11,RKI,EE
   real A,AA,Y,R,R0R
   real REC_Y,REC_101325,REC_35
   real ZEPSCQ,ZEPSCT,ZEPSC,RAYL,XNET
   external WFLUX
   external SUN_RADFIX1
   logical RADFIX,RADFLTR

   real DSIG(LMX,KMX),DSH(LMX,KMX),DSC(LMX,KMX),DZ(LMX,KMX)
   real VOZO(LMX)

   real PSOL(LMX),WV(MM,KMX),TM(MM,KMX),QOF(LMX,KMX), &
        CNEB(LMX,KMX),TAUAE(LMX,KMX,5) ,RMUO(LMX),ALBS(LMX)
   real aSIG(LMX,KMX)
   real Z(LMX)

   real FDOWN(LMX,flkmxp),FUP(LMX,flkmxp),HEAT(LMX,flkmx)
   real OZOTOIT(LMX)

   real UD(LMX,3,KMXP),UM(LMX,KMXP),REFZ(LMX,2,KMXP), &
        RJ(LMX,6,KMXP),RK(LMX,6,KMXP),TR(LMX,2,KMXP), &
        C1I(LMX,KMXP),RL(LMX,8),RMUE(LMX,KMXP),RAY(LMX,KMXP), &
        RUEF(LMX,8),RE1(LMX),RE2(LMX),TR1(LMX),TE2(LMX),S(LMX),G(LMX), &
        R23(LMX),REF(LMX),RMU(LMX),R1(LMX),W(LMX),TO1(LMX), &
        GG(LMX),RNEB(LMX),RMUZ(LMX),TRA(LMX,KMXP), &
        TAUAZ(LMX,KMX),PIZAZ(LMX,KMX),CGAZ(LMX,KMX), &
        OMEGAT(LMX,KMX),CG(LMX,KMX),TAU(LMX,KMX), &
        UDI1(LMX),UDI2(LMX)
   integer KK(LMX)

   logical reduc
   real flss(lmx,flkmxp), srd(lmx,kmxp)
   real afldsig(lmx,flkmx), aflsig(lmx,flkmx)

   !@Author L.Garand (1989)

   !@Revision
   ! 001      G.Pellerin(Mar90)Standard documentation
   ! 002      N. Brunet  (May91)
   !                New version of thermodynamic functions
   !                and file of constants
   ! 003      B. Bilodeau  (August 1991)- Adaptation to UNIX
   ! 004      C. Girard (Nov 1992) - Adjustment on parameters
   !          controlling the effect of clouds
   ! 005      B. Bilodeau (Apr 1993) - Optimization
   ! 006      B. Bilodeau (May 1993) - R0 variable
   ! 007      R. Benoit (Aug 93) Local Sigma
   ! 008      B. Bilodeau (November 1993) - Total ozone
   !          (cm stp) above model roof in ozotoit
   !          Change name from SUN1 to SUN2
   ! 009      Wei Yu (Aug 94) - New option KUO2SUN2
   !          Change name from SUN2 to SUN3
   ! 010      M. Gagnon (June 1995) - Reduction mode
   ! 011      Louis Garand (April 1995) - Remove calculations of
   !          optical parameters now done in CLDOPTX; minor change
   !          to cloud overlapping calculation (as in IR code)
   ! 012      B. Dugas (Sep 96) - RADFIX to control 1) fixes and
   !          2) loss of solar flux due to ozone above model lid
   ! 013      V. Lee (Apr 99) - Correct bug (loop 904)
   ! 014      G. Lemay, A. Patoine and B. Bilodeau (Sep 99) - Correct
   !                    multitasking bug by removing comdeck "solfact"
   ! 015      B. Bilodeau (Jan 01) - Automatic arrays
   ! 016      M. Lepine (March 2003) -  CVMG... Replacements
   ! 017      D. Talbot (June 2003) - IBM conversion
   !               - calls to vsexp routine (from massvp4 library)
   !               - calls to exponen4 (to calculate power function '**')
   !               - divisions replaced by reciprocals
   ! 018      A-M. Leduc (Jun 2003)  - Add RADFLTR switch (sun6 to sun7)
   ! 019      F. Lemay (Sept 2003)   - Call to sun_radfix
   ! 020      B. Bilodeau (Aug 2003) - exponen4 replaced by vspown1
   ! 021      B. Bilodeau (Mar 2004) - add double loop INSIDE sun_radfix

   !@Object
   !          to calculate the solar heating rates and fluxes after
   !          Y.Fouquart and B. Bonnel, 1980, Beitr. Phys. Atmosph. 53,1,
   !          35-62

   !@Arguments
   !          - Output -
   ! FDOWN    downward flux at each flux level (W/m2) as output
   ! FUP      upward flux at each flux level (W/m2) as output
   ! HEAT     solar heat in Kelvin/second as output
   ! OZOTOIT  total ozone (cm stp) above model roof
   ! VOZO     transmissivity of ozone layer above model roof
   !          - Input -
   ! DSIG     sigma thickness of each level
   ! DSH      sigma thickness affected by exponent 1.9
   ! DSC      sigma thickness affected by exponent 1.75
   ! DZ       thickness of each layer in metres
   ! PSOL     surface pressure in Newtons/m2
   ! TM       temperature in middle of layer
   ! WV       mixing water vapour ratio in kg/kg (middle of layer)
   ! ASIG     reduced sigma levels
   ! QOF      amount of ozone in cm STP in each layer
   ! CNEB     cloud fraction in each layer
   ! ZLWC     amount of liquid water per layer in kg/m3
   ! TAUAE    optical thickness of aerosols per layer (no units)
   ! RMUO     cosines of the solar Zenith angle (1. at noon)
   ! ALBS     surface albedo
   ! LMX      maximum number of points to execute (identical to the 1st
   !          dimension of defined multi-dimensional tables)
   ! KMX      number of layers
   ! KMXP     number of flux levels (KMX+1)
   ! NPTS     number of points requested with ILG<=LMX
   ! MM       1st dimension of TM and WV
   ! RJ       downward fluxes
   ! RK       upward fluxes
   ! UD       downward radiation
   ! REFZ     reflectivities:
   !          REFZ(1,KM1) for reflectivity of layer KM1 for a non-
   !          reflecting underlying layer (R=0);
   !          REFZ(2,KM1) for reflectivity of layer KM1 for a reflecting
   !          underlying layer (R is not 0)
   ! TR       transmissivities
   ! UM       upward radiation
   ! RL       work space
   ! RUEF     work space
   ! C1I      total cloudiness above level K assuming a random
   !          overlapping of the cloud layers, taking into account the
   !          actual optical thickness of the cloud (via the forward
   !          scattering peak)
   ! RMUE     work space
   ! RAY      Rayleigh scattering
   ! TRA      work space
   ! TAUAZ    work space
   ! PIZAZ    work space
   ! CGAZ     work space
   ! OMEGAT   cloud layer single scattering albedo (0 to 1)
   ! CG       cloud layer asymmetry parameter (0 to 1)
   ! TAU      cloud layer optical thickness (no dimension, >=0)
   ! RE1      work space
   ! RE2      work space
   ! TR1      work space
   ! TE2      work space
   ! S        work space
   ! G        work space
   ! R23      work space
   ! REF      work space
   ! RMU      work space
   ! R1       work space
   ! W        work space
   ! TO1      work space
   ! GG       work space
   ! RNEB     work space
   ! RMUZ     work space
   ! KK       work space
   ! reduc    .true. to use interpolation
   !          .false. means we are working on full levels
   ! flss     full "flux" sigma levels
   ! srd      reduced "flux" sigma levels
   ! afldsig  thickness between full "flux" sigma levels
   ! aflsig   full sigma levels
   ! flkmx    full number of levels excluding ground
   ! flkmxp   full number of levels including ground
   ! RADFIX   .true. to use curve fit at the end
   ! RADFLTR  .true. to apply smoothing of net fluxes
   ! R0R      factor close to 1.0 (+/- 3%) that takes into account
   !          the variation in sun-earth distance

   !@Notes
   !          Modified for inclusion of aerosol effect by D. Tanre, Aug.
   !          1982.Vectorized and optimized by J.J. Morcrette , March
   !          1984.Provided by Jean-Pierre Blanchet, CCRN, AES,
   !          Toronto.This is the version with one spectral interval.
   !          Version with two intervals (O.25-0.68 and 0.68-4Microns)
   !          is also available.


   !-----------------------------------------------------------------------

   logical, dimension(    FLKMX ) :: TOIT

   real, dimension(LMX,FLKMX ) :: FLDSIG
   real, dimension(LMX,FLKMX ) :: FLSIG
   real, dimension(LMX,KMX   ) :: SIG
   real, dimension(LMX,FLKMXP) :: FLFDOWN
   real, dimension(LMX,FLKMXP) :: FLFUP

   real CH2O,CCO2
   real TAUA(5),PIZA(5),CGA(5)  !#,CAER(4,5)
   real APAD(3,6),BPAD(3,6),DQ(3),AKI(2)
   save TAUA,PIZA,CGA,APAD,BPAD,AKI,DQ !#,CAER
   save CH2O,CCO2
   logical LO1

   real SOMME,DIV,P1,P2,P3,P4,P5,P6,NUM,DENOM,ZY,NUME,DENOMI

   data DQ/ 0.61 , 0.93 , 0.00 /

   data AKI / 0.457, 0.00636 /

   data ((APAD(I,J),J=1,6),I=1,3)/ 6.653722092E-05, 1.228372707E-01, &
        1.434798127E+01, 1.028708866E+02, 3.309603691E+01, 0.0 , &
        6.074323894E+06, 2.379915061E+08, &
        1.978907838E+08, 7.973068584E+06, 1.221333677E+04, 0.0 , &
        4.875774044E-04, 6.493865154E-01, &
        8.363354452E+01, 6.305273645E+02, 7.929825858E+01, 0.0 /

   data ((BPAD(I,J),J=1,6),I=1,3) /6.653722092E-05, 1.242963636E-01, &
        1.534122504E+01, 1.289549216E+02, 6.193430339E+01, 1.0 , &
        6.074323894E+06, 2.392851737E+08, &
        2.070297278E+08, 9.314161596E+06, 1.966084935E+04, 1.0 , &
        4.875774044E-04, 6.499019064E-01, &
        8.435265650E+01, 6.443741405E+02, 9.864863951E+01, 1.0 /


   !** OPTICAL PARAMETERS FOR THE STANDARD AEROSOLS OF THE RADIATION
   !   COMMISSION

   !MODEL 1= CONTINENTAL, 2=MARITIME, 3=URBAN, 4=VOLCANIC, 5=STRATOSPHERIC


   !-- SHORTWAVE (ADAPTED TO THE L.O.A. SHORTWAVE SCHEME)

   data TAUA / .730719, .912819, .725059, .745405, .682188 /
   data PIZA / .872212, .982545, .623143, .944887, .997975 /
   data CGA  / .647596, .739002, .580845, .662657, .624246 /

   !-- LONGWAVE (ADAPTED TO THE L.O.A. LONGWAVE SCHEME SPECTRALS INTERVALS)

!!$      DATA CAER / .038520, .037196, .040532, .054934, &
!!$                  .12613 , .18313 , .10357 , .064106, &
!!$                  .012579, .013649, .018652, .025181, &
!!$                  .011890, .016142, .021105, .028908, &
!!$                  .013792, .026810, .052203, .066338 /

   data CH2O/5.3669274E-3/, CCO2/3.3E-4/

   !     FONCTIONS-FORMULES POUR LE CALCUL DES APPROXIMANTS DE PADE

   SOMME(P1,P2,P3,P4,P5,P6,ZY)= &
        (((((P6*ZY+P5)*ZY+P4)*ZY+P3)*ZY+P2)*ZY+P1)

   DIV(NUM,DENOM,ZY)=(NUM/DENOM)*(1.-ZY) + ZY

   !    CAUTION: ALL ENTRIES ARE INVERTED VERTICALLY

   ILG=NPTS
   if(ILG.gt.LMX)then
      write(6,9001)
9001  format(1X,' ILG DOIT ETRE < OU = A LMX DANS SUN2: FATAL')
   endif


   KFIN=flKMX*.5 + 1

   do K=1,KFIN
      KI=flKMX-K+1

      do J=1,LMX
         flsig(J,K)=aflsig(J,KI)
         flsig(J,KI)=aflsig(J,K)
         fldsig(J,K)=afldsig(J,KI)
         fldsig(J,KI)=afldsig(J,K)
      enddo
   enddo

   KFIN=KMX*.5
   do K=1,KFIN
      KI=KMX-K+1
      do J=1,LMX
         A=WV(J,K)
         WV(J,K)=WV(J,KI)
         WV(J,KI)=A
         A=DZ(J,K)
         DZ(J,K)=DZ(J,KI)
         DZ(J,KI)=A
         A=TM(J,K)
         TM(J,K)=TM(J,KI)
         TM(J,KI)=A
         A=QOF(J,K)
         QOF(J,K)=QOF(J,KI)
         QOF(J,KI)=A
         A=CNEB(J,K)
         CNEB(J,K)=CNEB(J,KI)
         CNEB(J,KI)=A
         A=cg(J,K)
         cg(J,K)=cg(J,KI)
         cg(J,KI)=A
         A=tau(J,K)
         tau(J,K)=max(tau(J,KI),1.e-10)
         tau(J,KI)=max(A,1.e-10)
         A=omegat(J,K)
         omegat(J,K)=omegat(J,KI)
         omegat(J,KI)=A

         A=TAUAE(J,K,1)
         TAUAE(J,K,1)=TAUAE(J,KI,1)
         TAUAE(J,KI,1)=A
         A=TAUAE(J,K,2)
         TAUAE(J,K,2)=TAUAE(J,KI,2)
         TAUAE(J,KI,2)=A
         A=TAUAE(J,K,3)
         TAUAE(J,K,3)=TAUAE(J,KI,3)
         TAUAE(J,KI,3)=A
         A=TAUAE(J,K,4)
         TAUAE(J,K,4)=TAUAE(J,KI,4)
         TAUAE(J,KI,4)=A
         A=TAUAE(J,K,5)
         TAUAE(J,K,5)=TAUAE(J,KI,5)
         TAUAE(J,KI,5)=A

      enddo
   enddo


   do I=1,KFIN
      do J=1,LMX
         IK=KMX-I+1
         A=DSIG(J,I)
         DSIG(J,I)=DSIG(J,IK)
         DSIG(J,IK)=A
         A=DSH(J,I)
         DSH(J,I)=DSH(J,IK)
         DSH(J,IK)=A
         A=DSC(J,I)
         DSC(J,I)=DSC(J,IK)
         DSC(J,IK)=A
         SIG(J,I)=aSIG(J,IK)
         SIG(J,IK)=aSIG(J,I)
      enddo
   enddo

   CCAR=CCO2*0.0088509
   CONS = GRAV/CPD


   !-----------------------------------------------------------------------


   !-- THRESHOLDS FOR CLOUDINESS, HUMIDITY, CLOUD OPTICAL THICKNESS

   ZEPSC =1.E-02
   ZEPSCQ=1.E-10
   ZEPSCT=1.E-03

   RAYL=0.06


   !         Y  : MAGNIFICATION FACTOR
   !         UD : DOWNWARD RADIATION
   !         UM : UPWARD RADIATION

   !     INTERACTIONS BETWEEN OZONE ABSORPTION AND SCATTERING ARE NEGLECTED

   !     ABSORBER AMOUNT IN THE K TH LAYER
   !     H2O : UD(1,K)
   !     CO2 : UD(2,K)

   IABS = 3

   REC_35=1./35.

   !     si on introduit ce vsqrt sun6 plante car on perd en precision
   !     CALL VSQRT(RMU,1224.*RMUO*RMUO+1.,ILG)
   !     DO I=1,ILG
   !     RMU(I)=SQRT(1224.*RMUO(I)*RMUO(I)+1.)*REC_35
   !     RMU(I)=RMU(I)*REC_35
   !     ENDDO

   do I=1,ILG
      C1I(I,KMXP)=0.
      UD(I,3,KMXP)=0.
      RMU(I)=sqrt(1224.*RMUO(I)*RMUO(I)+1.)*REC_35
      W(I) = OZOTOIT(I)/RMU(I)
      NUME  = SOMME(APAD(IABS,1),APAD(IABS,2),APAD(IABS,3),APAD(IABS,4), &
           APAD(IABS,5),APAD(IABS,6),W(I))
      DENOMI= SOMME(BPAD(IABS,1),BPAD(IABS,2),BPAD(IABS,3),BPAD(IABS,4), &
           BPAD(IABS,5),BPAD(IABS,6),W(I))

      !     POUR TENIR COMPTE DE L'ABSORPTION DU RAYONNEMENT SOLAIRE
      !     CAUSEE PAR L'OZONE AU-DESSUS DU TOIT DU MODELE, REACTIVER
      !     LA LIGNE SUIVANTE. ELLE L'EST DEJA SI RADFIX EST FAUX
      if (RADFIX) then
         vozo(i) = 1.
      else
         vozo(i) = DIV(NUME,DENOMI,DQ(IABS))
      endif
   enddo

   !     OZONE ABSORPTION ABOVE MODEL ROOF (NON-OPTIMIZED CODE)
   !     CALL TTTT(IABS,W,VOZO,WORKX,LMX,ILG,APAD,BPAD,DQ,AKI)


   !-----------------------------------------------------------------------
   !-- CALCULATES OZONE AMOUNTS FOR DOWNWARD LOOKING PATHS (SEC=1./RMU)

   do K=1,KMX
      KP1=K+1
      L=KMXP-K
      LP1=L+1
      do I=1,ILG
         UD(I,3,L)=UD(I,3,LP1)+QOF(I,L)/RMU(I)
      enddo
   enddo

   !-----------------------------------------------------------------------
   !  CALCULATES OZONE AMOUNTS FOR UPWARD LOOKING PATHS (SEC=1.66)

   Y=.6024
   REC_Y=1./Y
   do I=1,ILG
      UM(I,1)=UD(I,3,1)
   enddo

   do K=2,KMXP
      KM1=K-1
      do I=1,ILG
         UM(I,K)=UM(I,KM1)+QOF(I,KM1)*REC_Y
      enddo
   enddo

   !-----------------------------------------------------------------------
   !   CALCULATES AMOUNTS OF WATER VAPOR AND UNIFORMLY MIXED GASES (O2+CO2)
   !   FOR A VERTICAL PATH

   REC_101325=1./101325.
   do K=1,KMX
      KP1=K+1
      do I=1,ILG
         !     TZ=TM(I,K)
         !     Z = TCDK/TZ
         Z(I) = TCDK/TM(I,K)
      enddo

      call VSPOWN1  (UDI1,Z,0.45,ILG)
      call VSPOWN1  (UDI2,Z,1.375,ILG)
      do I=1,ILG
         WH2O=AMAX1(WV(I,K),1.E-8)
         DSHH= DSH(I,K) * CH2O*WH2O
         !     DSCC= DSC(I,K) * CCAR*REC_101325
         DSCC=DSC(I,K)*CCAR
         !     UD(I,1,K)= PSOL(I) * DSHH * Z**0.45
         !     UD(I,2,K)= PSOL(I) * DSCC * Z**1.375
         !     UD(I,1,K)= PSOL(I) * DSHH * (exp(0.45 *log(Z)))
         !     UD(I,2,K)= PSOL(I) * DSCC * (exp(1.375*log(Z)))
         UD(I,1,K)= PSOL(I) * DSHH * UDI1(I)
         UD(I,2,K)= PSOL(I) * DSCC * UDI2(I)
      enddo
   enddo


   !-----------------------------------------------------------------------
   !   C1I(*) TOTAL CLOUDINESS ABOVE LEVEL K ASSUMING A RANDOM OVERLAPPING
   !    OF THE CLOUD LAYERS, TAKING INTO ACCOUNT THE ACTUAL OPTICAL
   !    THICKNESS OF THE CLOUD (VIA THE FORWARD SCATTERING PEAK)


   do I=1,ILG
      R23(I)=0.
      C1I(I,KMXP)=0.
   enddo
   !     DO 11 K=1,KMX
   !     L=KMXP-K
   !     DO 10 I=1,ILG
   !-- 0.28 = 1.- PIZERO*G*G
   !     CORR=1.-EXP(-0.28*TAU(I,L)/RMU(I))
   !     R23(I)=1.-(1.-CNEB(I,L)*CORR)*(1.-R23(I))
   !     C1I(I,L)=R23(I)
   !     IF(I.EQ.1)WRITE(6,77)L,C1I(I,L),TAU(I,L),CNEB(I,L)
   ! 77  FORMAT(1X,'C1I TAU NEB:',I6,3E12.4)
   !10   CONTINUE
   !11   CONTINUE




   do i=1,ilg
      ! te2 is memory of level for new cloud layer
      ! where maximum overlap occurs
      ! tr1 records the mimimum transmittance in that cloud
      te2(i)=1.
      tr1(i)=1.
   end do

   do J=1,KMXP
      do I=1,ILG
         RAY(I,J)=1.
      enddo
   enddo
   do J=2,KMXP
      LIND=MAX0(1,KMX-(J-2)+1)
      LIND=MIN0(LIND,KMX)
      do I=1,ILG
         K=KMX-(J-1)+1
         XNU=1.-CNEB(I,K)
         VA=CNEB(I,LIND)
         if(va.lt.0.01)then
            te2(i)=ray(i,k)
            tr1(i)=xnu
         else
            tr1(i)=min(tr1(i),xnu)
         endif
         ray(i,j)=tr1(i)*te2(i)
      enddo
   enddo
   !  RAY EST LE VECTEUR NUAGES DU SOMMET AUX AUTRES NIVEAUX
   !  AVEC RANDOM OVERLAP POUR COUCHES SEPAREES ET FULL OVERLAP
   !  POUR COUCHES ADJACENTES

   do K=1,KMXP
      do I=1,ILG
         FUP(I,K)=1.-RAY(I,KMXP-K+1)
         !     IF(I.EQ.1)WRITE(6,9004)K,FUP(I,K)
         !9004 FORMAT(' NOUV C1I: ',I5,2X,E12.4)
      enddo
   enddo


   do K = 1, KMX
      do I=1,ILG
         RE2(I)=0.
         RE1(I)=0.
         C1I(I,K)=FUP(I,K)
         FDOWN(I,K)=CNEB(I,K)
         PIZAZ(I,K) = 0.
      enddo
   enddo

   !     IF(1.EQ.1)GO TO 630
   !     DO 412 K=1,KMX
   !        L=KMXP-K
   !        K P 1 = K + 1
   !        DO 411 I=1,ILG
   !           C NEB 1 =CNEB(I,KMX-K+1)
   ! . . . . COMPUTES A DOWNWARD REDUCED CLOUDINESS TO CANCELL RANDOM OVERL
   !           C NEB 2 = CVMGT (0.,(RAY(I,KP1)-RAY(I,K) )/
   !    1                      AMAX1(1. -RAY(I,K), 0.001)    ,
   !    2                        RAY(I,K) .GE. 0.999          )
   !           FDOWN(I,L) = (1. - RE1(I)) * C NEB 1 + RE1(I) * C NEB 2
   !           FDOWN(I,L) = AMIN1(FDOWN(I,L),CNEB1)
   !           RE2(I) = RE2(I) + TAU(I,L)
   !           PIZAZ(I,L) = PIZAZ(I,L) + TAU(I,L)*(CNEB1-FDOWN(I,L))
   !           FACOR = 1. - OMEGAT(I,L)*CG(I,L)*CG(I,L)
   !           RE1(I) = 1. - EXP(-FACOR * RE2(I) / RMU(I))
   !           C1I(I,L) = RAY(I,KP1) * RE1(I)
   !411     CONTINUE
   !412  CONTINUE

   !     DO 422 I=1,ILG
   !422    RE1(I)=0.
   !     DO 423 J=1,KMX
   !     DO 424 I=1,ILG
   !424  RE1(I)=RE1(I)+FDOWN(I,J)
   !423    CONTINUE

   !     KMX M 1 = KMX - 1
   !     DO 416 K = 1, KMX M 1
   !         L = KMXP - K
   !         DO 415 I = 1,ILG
   !           FUP(I,L) = TAU (I,L)
   !           IF (FDOWN(I,L) .GT. .01)THEN
   !             TAU (I,L) = (TAU(I,L) * FDOWN(I,L) + PIZAZ(I,L-1)) /
   !    1                                FDOWN(I,L)
   !           TAU(I,L)= FDOWN(I,L)*RE2(I)/RE1(I)
   !           tau(i,l) = min(tau(i,l),150.)
   !             OMEGAT(I,L) = 1. - (1.-OMEGAT(I,L))*FUP(I,L)/TAU(I,L)
   !             OMEGAT(I,L) = AMAX1(OMEGAT(I,L),0.9936)
   !             omegat(i,l) = amin1 ( omegat(i,l) , 0.999999 )
   !           END IF
   ! 415     CONTINUE
   ! 416 CONTINUE

   !630     CONTINUE



   !-----------------------------------------------------------------------
   !-- OPTICAL PARAMETERS FOR THE MIXING OF AEROSOLS

   do K=1,KMX
      do I=1,ILG
         CGAZ(I,K)=0.
         PIZAZ(I,K)=0.
         TAUAZ(I,K)=0.
      enddo
      do IAE=1,5
         do I=1,ILG
            TAUAZ(I,K)=TAUAZ(I,K)+TAUAE(I,K,IAE)*TAUA(IAE)
            PIZAZ(I,K)=PIZAZ(I,K)+TAUAE(I,K,IAE)*TAUA(IAE)*PIZA(IAE)
            CGAZ(I,K)=CGAZ(I,K)+TAUAE(I,K,IAE)*TAUA(IAE)*PIZA(IAE)*CGA(IAE)
         enddo
      enddo
      do I=1,ILG
         CGAZ(I,K)=CGAZ(I,K)/PIZAZ(I,K)
         PIZAZ(I,K)=PIZAZ(I,K)/TAUAZ(I,K)
      enddo
   enddo

   !------------------------------------------------------------------
   !  REFLECTIVITIES AND TRANSMITTIVITIES ARE FIRST CALCULATED WITHOUT
   !   GASEOUS ABSORPTION, I.E., FOR A PURELY SCATTERING ATMOSPHERE
   !   INCLUDING RAYLEIGH, CLOUDS, AEROSOLS

   !     SURFACE CONDITIONS    ALBS

   do I=1,ILG
      REFZ(I,2,1)=ALBS(I)
      REFZ(I,1,1)=ALBS(I)
   enddo

   Y=1.66
312 format( 5X,'REFZ1',10X,'REFZ2',12X,'TR1',13X,'TR2',/)

   do K=2,KMXP
      KM1=K-1
      do I=1,ILG
         RNEB(I)=FDOWN(I,KM1)
         !      RE1(I)=0.
         !      TR1(I)=0.
         !      RE2(I)=0.
         !      TE2(I)=0.

         !--   EQUIVALENT ZENITH ANGLE OBTAINED AS A MIX OF THE ZENITH ANGLE FOR
         !     DIRECT RADIATION AND OF THE DIFFUSE (SEC=1.66) ZENITH ANGLE
         !     DEPENDING ON THE OVERHEAD CLOUDINESS

         XMUE=(1.-C1I(I,K))/RMU(I)+C1I(I,K)*1.66

         RMUE(I,K)=1./XMUE

         !     REFLECTIVITY OF LAYER KM1 DUE TO RAYLEIGH SCATTERING

         R=RAYL*DSIG(I,KM1)

         !-- ORIGINAL FORMULA (FOUQUART)
         !     RAY(I,KM1)=0.5*R/(RMUE(I,K)+R)

         !-- TANRE'S FORMULA
         !     RAY(I,KM1)= R / (2.*RMUE(I,K) + R)

         !-- BURIEZ' FORMULA, JUIN 82
         !     R=DSIG(I,KM1)
         !     RAY(I,KM1)=0.0536*R/(RMUE(I,K)+0.063+0.055*R)

         !-- TANRE'S FORMULA MODIFYING THE CONTRIBUTION OF RAYLEIGH SCATTERING
         !    TO ACCOUNT FOR AEROSOLS

         PT=R+PIZAZ(I,KM1)*TAUAZ(I,KM1)
         XMPT=(1.-PIZAZ(I,KM1))*TAUAZ(I,KM1)
         BPT=0.5*(1.-TAUAZ(I,KM1)*CGAZ(I,KM1)/(R+TAUAZ(I,KM1)))*PT
         TRA(I,KM1)=1./(1.+XMPT*XMUE+BPT*XMUE)
         RAY(I,KM1)=TRA(I,KM1)*BPT*XMUE

         !-- INTRODUCING THE EFFECT OF A CLOUD LAYER IN THE SCATTERING PROCESS

         K0=0
         K1=1
         LO1=RNEB(I).gt.ZEPSC
         if (LO1) then
            ZC = 0.
         else
            ZC = 1.
         endif
         KK(I)=int(ZC)
         if (LO1) then
            TAU(I,KM1) = AMAX1(TAU(I,KM1),ZEPSCT)
         else
            TAU(I,KM1) = 0.
         endif

         W(I) = OMEGAT(I,KM1)
         TO1(I)=TAU(I,KM1)/W(I)
         REF(I)=REFZ(I,1,KM1)
         GG(I)=CG(I,KM1)
         RMUZ(I)=RMUE(I,K)
      enddo

      call WFLUX (REF, GG, W, TO1, RE1, TR1, RMUZ, RE2, TE2,LMX,ILG)

      !    REFZ(1,KM1): REFLECTIVITY OF LAYER KM1
      !       FOR A NON-REFLECTING UNDERLYING LAYER  R=0.
      !    REFZ(2,KM1): FOR A REFLECTING (R#0.) UNDERLYING LAYER

      do I=1,ILG
         REFZ(I,1,K)=(1.-RNEB(I))*(RAY(I,KM1)+REFZ(I,1,KM1)* &
              (  TRA(I,KM1) ))+ RNEB(I)*RE2(I)
         TR(I,1,KM1)=RNEB(I)*TE2(I)+(  TRA(I,KM1) )*(1.-RNEB(I))
         REFZ(I,2,K)=(1.-RNEB(I))*(RAY(I,KM1)+KK(I)*REFZ(I,2,KM1)* &
              (  TRA(I,KM1) ))+ RNEB(I)*RE1(I)
         TR(I,2,KM1)=RNEB(I)*TR1(I)+(1.-RNEB(I))*(  TRA(I,KM1) )
      enddo
   enddo

   do JJ=1,2

      !  RJ : DOWNWARD
      !  RK : UPWARD

      do I=1,ILG
         RJ(I,JJ,KMXP)=1.
         RK(I,JJ,KMXP)=REFZ(I,JJ,KMXP)
      enddo

      do K=1,KMX
         L=KMXP-K
         LP1=L+1
         do I=1,ILG
            RE11=RJ(I,JJ,LP1)*TR(I,JJ,L)
            RJ(I,JJ,L)=RE11
            RK(I,JJ,L)=RE11*REFZ(I,JJ,L)
         enddo
      enddo
   enddo

   !----------------------------------------------------------------
   !  REFLECTIVITIES AND TRANSMITTIVITIES ARE NOW CALCULATED WITH A
   !  SIMULATED GASEOUS ABSORPTION USING AN EMPIRICAL GREY ABSORPTION
   !  COEFFICIENTS  AKI(IABS)

   !  IABS=1  WATER VAPOR  ; IABS=2  UNIFORMLY MIXED GASES

   N=2

   DO28: do IABS=1,2
      RKI=AKI(IABS)

      do I=1,ILG
         REFZ(I,2,1)=ALBS(I)
         REFZ(I,1,1)=ALBS(I)
      enddo

      DO23: do K=2,KMXP
         KM1=K-1
         do i=1,ilg
            AA=UD(I,IABS,KM1)
            TO1(I) = RKI*aa
            s(i)=-RKI*AA*y
            g(i)=-RKI*AA/RMUE(I,K)
            RMUZ(I)=RMUE(I,K)
         enddo
         call vsexp(s,s,ilg)
         call vsexp(g,g,ilg)
         do I=1,ILG
            RNEB(I)=FDOWN(I,KM1)
            !      AA=UD(I,IABS,KM1)
            !      S(I)=EXP(-RKI*AA*Y)
            !      G(I)=EXP(-RKI*AA/RMUE(I,K))
            !      TR1(I)=0.
            !      RE1(I)=0.
            !      TE2(I)=0.
            !      RE2(I)=0.

            LO1=RNEB(I).gt.ZEPSC
            if (LO1) then
               ZC = 0.
            else
               ZC = 1.
            endif
            KK(I)=int(ZC)
            if (LO1) then
               TAU(I,KM1) = AMAX1(TAU(I,KM1),ZEPSCT)
            else
               TAU(I,KM1) = 0.
            endif


            !      TO1(I) = RKI*AA + TAU(I,KM1)/OMEGAT(I,KM1)
            TO1(I) = TO1(i) + TAU(I,KM1)/OMEGAT(I,KM1)
            W(I)=TAU(I,KM1)/TO1(I)
            REF(I)=REFZ(I,1,KM1)
            GG(I)=CG(I,KM1)
            !      RMUZ(I)=RMUE(I,K)
         enddo

         call WFLUX (REF, GG, W, TO1, RE1, TR1, RMUZ, RE2, TE2,LMX,ILG)

         do I=1,ILG
            REFZ(I,2,K)=(1.-RNEB(I))*(RAY(I,KM1)+KK(I)*REFZ(I,2,KM1)* &
                 (  TRA(I,KM1) ))*G(I)*S(I)+RNEB(I)*RE1(I)
            TR(I,2,KM1)=(1.-RNEB(I))*(  TRA(I,KM1) )*G(I)+RNEB(I)*TR1(I)
            REFZ(I,1,K)=(1.-RNEB(I))*(RAY(I,KM1)+REFZ(I,1,KM1)* &
                 (  TRA(I,KM1) ))*G(I)*S(I)+RNEB(I)*RE2(I)
            TR(I,1,KM1)=(1.-RNEB(I))*(  TRA(I,KM1) )*G(I)+RNEB(I)*TE2(I)
         enddo
      enddo DO23

      do KREF=1,2
         N=N+1

         do I=1,ILG
            RJ(I,N,KMXP)=1.
            RK(I,N,KMXP)=REFZ(I,KREF,KMXP)
         enddo

         do K=1,KMX
            L=KMXP-K
            LP1=L+1
            do I=1,ILG
               RE11=RJ(I,N,LP1)*TR(I,KREF,L)
               RJ(I,N,L)=RE11
               RK(I,N,L)=RE11*REFZ(I,KREF,L)
            enddo
         enddo
      enddo
   enddo DO28



   !-----------------------------------------------------------------------
   !   UPWARD (RK) AND DOWNWARD (RJ) FLUXES ARE NOW CALCULATED
   !    WITHOUT GASEOUS ABSORPTION     : N= 1 , 2
   !    WITH GASEOUS ABSORPTION BY H2O : N= 3 , 4
   !    WITH GASEOUS ABSORPTION BY CO2 : N= 5 , 6
   !    N= 1, 3, 5 : NON-REFLECTING UNDERLYING LAYER
   !    N= 2, 4, 6 :  REFLECTING UNDERLYING LAYER

   !  EFFECTIVE ABSORBER AMOUNTS    UEFF= - LN (F(K)/F(0))/K

   !----- EE IS AN ADJUSTABLE THRESHOLD DEPENDING ON THE COMPUTER PRECISION
   !       ADDED TO 'RJ' TO PREVENT LOGARITHM OF ZERO

   EE=1.E-10
   do K=1,KMXP
      do JJ=1,5,2
         JJP =JJ+1
         do I=1,ILG
            RJ(I,JJ,K)=        RJ(I,JJ,K) - RJ(I,JJP ,K)
            RK(I,JJ,K)=        RK(I,JJ,K) - RK(I,JJP ,K)
            RJ(I,JJ,K)= AMAX1( RJ(I,JJ,K) , EE )
            RK(I,JJ,K)= AMAX1( RK(I,JJ,K) , EE )
         enddo
      enddo
   enddo

   do K=1,KMXP
      do JJ=2,6,2
         do I=1,ILG
            !     RJ(I,JJ,K)= EE + RJ(I,JJ,K)
            !     RK(I,JJ,K)= EE + RK(I,JJ,K)
            RJ(I,JJ,K)= AMAX1( RJ(I,JJ,K) , EE )
            RK(I,JJ,K)= AMAX1( RK(I,JJ,K) , EE )
         enddo
      enddo
   enddo

   DO42: do K=1,KMXP
      KKI=1

      DO40: do JJ=1,2
         RKI=1.0/AKI(JJ)

         DO39: do N=1,2
            N2J=N+2*JJ
            KKP4=KKI+4

            !  EFFECTIVE ABSORBER AMOUNT
            do I=1,ILG
               !      W(I)=ALOG(RJ(I,N,K)/RJ(I,N2J,K))/RKI
               W(I)=RJ(I,N,K)/RJ(I,N2J,K)
            enddo
            call vslog(w,w,ilg)
            do i=1,ilg
               w(i)=w(i)*rki
               NUME  = SOMME(APAD(JJ,1),APAD(JJ,2),APAD(JJ,3),APAD(JJ,4), &
                    APAD(JJ,5),APAD(JJ,6),W(I))
               DENOMI= SOMME(BPAD(JJ,1),BPAD(JJ,2),BPAD(JJ,3),BPAD(JJ,4), &
                    BPAD(JJ,5),BPAD(JJ,6),W(I))
               R1(I) = DIV(NUME,DENOMI,DQ(JJ))
               RUEF(I,KKI)=W(I)
            enddo

            !  TRANSMISSION FUNCTION
            !     CET APPEL A ETE REMPLACE PAR LES FONCTIONS-FORMULE SOMME ET DIV.
            !     CALL TTTT(JJ,W,R1,WORKX,LMX,ILG,APAD,BPAD,DQ,AKI)

            do I=1,ILG
               RL(I,KKI)=R1(I)
               !     W(I)=ALOG(RK(I,N,K)/RK(I,N2J,K))/RKI
               W(I)=RK(I,N,K)/RK(I,N2J,K)
            enddo
            call vslog(w,w,ilg)
            do i=1,ilg
               w(i)=w(i)*rki
               NUME  = SOMME(APAD(JJ,1),APAD(JJ,2),APAD(JJ,3),APAD(JJ,4), &
                    APAD(JJ,5),APAD(JJ,6),W(I))
               DENOMI= SOMME(BPAD(JJ,1),BPAD(JJ,2),BPAD(JJ,3),BPAD(JJ,4), &
                    BPAD(JJ,5),BPAD(JJ,6),W(I))
               R1(I) = DIV(NUME,DENOMI,DQ(JJ))
            enddo

            !     CET APPEL A ETE REMPLACE PAR LES FONCTIONS-FORMULE SOMME ET DIV.
            !     CALL TTTT(JJ,W,R1,WORKX,LMX,ILG,APAD,BPAD,DQ,AKI)

            do I=1,ILG
               RL(I,KKP4)=R1(I)
               RUEF(I,KKP4)=W(I)
            enddo

            KKI=KKI+1
         enddo DO39
      enddo DO40


      ! UPWARD AND DOWNWARD FLUXES WITH H2O AND CO2 ABSORPTION

      do I=1,ILG
         FDOWN(I,K)=RJ(I,1,K)*RL(I,1)*RL(I,3) &
              +RJ(I,2,K)*RL(I,2)*RL(I,4)
         FUP(I,K)=RK(I,1,K)*RL(I,5)*RL(I,7) &
              +RK(I,2,K)*RL(I,6)*RL(I,8)
      enddo
   enddo DO42

   !-----------------------------------------------------------------------
   !    OZONE ABSORPTION

   IABS=3
   do K=1,KMXP
      do I=1,ILG
         W(I)=UD(I,3,K)
         NUME  = SOMME(APAD(IABS,1),APAD(IABS,2),APAD(IABS,3),APAD(IABS,4), &
              APAD(IABS,5),APAD(IABS,6),W(I))
         DENOMI= SOMME(BPAD(IABS,1),BPAD(IABS,2),BPAD(IABS,3),BPAD(IABS,4), &
              BPAD(IABS,5),BPAD(IABS,6),W(I))
         R1(I) = DIV(NUME,DENOMI,DQ(IABS))
      enddo

      !     CET APPEL A ETE REMPLACE PAR LES FONCTIONS-FORMULE SOMME ET DIV.
      !     CALL TTTT(IABS,W,R1,WORKX,LMX,ILG,APAD,BPAD,DQ,AKI)

      do I=1,ILG
         FDOWN(I,K) = R1(I)*CONSOL*R0R*VOZO(I)*RMUO(I)*FDOWN(I,K)
         W(I)=UM(I,K)
         NUME  = SOMME(APAD(IABS,1),APAD(IABS,2),APAD(IABS,3),APAD(IABS,4), &
              APAD(IABS,5),APAD(IABS,6),W(I))
         DENOMI= SOMME(BPAD(IABS,1),BPAD(IABS,2),BPAD(IABS,3),BPAD(IABS,4), &
              BPAD(IABS,5),BPAD(IABS,6),W(I))
         R1(I) = DIV(NUME,DENOMI,DQ(IABS))
      enddo

      !     CET APPEL A ETE REMPLACE PAR LES FONCTIONS-FORMULE SOMME ET DIV.
      !     CALL TTTT(IABS,W,R1,WORKX,LMX,ILG,APAD,BPAD,DQ,AKI)

      do I=1,ILG
         FUP(I,K) = R1(I)*CONSOL*R0R*RMUO(I)*FUP(I,K)
      enddo
   enddo

   !**** NET FLUX AND HEATING RATES (DEG.K / DAY)

   do I=1,ILG
      !     ALBP(I)=FUP(I,KMXP)/FDOWN(I,KMXP)
      !     FSOL(I)=FDOWN(I,1)-FUP(I,1)
   enddo

   do K=1,KMXP
      do I=1,ILG
         !     Q(I,K)=FUP(I,K)-FDOWN(I,K)
      enddo
   enddo

   !     En mode reduction, interpoler fup et fdown aux niveaux de flux complets

   if( reduc ) then

      !       La condition .gt. de intchamps exige qu'on passe
      !       des champs dans le sens normal (1=haut,kmx=bas).
      !       Il faut inverser les flux a interpoler.

      IFIN=KMXP*.5

      do I=1,IFIN
         KI=KMXP-I+1

         do J=1,LMX
            A=FUP(J,I)
            FUP(J,I)=FUP(J,KI)
            FUP(J,KI)=A
            A=FDOWN(J,I)
            FDOWN(J,I)=FDOWN(J,KI)
            FDOWN(J,KI)=A
         enddo
      enddo

      call intchamps(flfdown,fdown,flss,srd,lmx,flkmxp,kmxp)
      call intchamps(flfup,fup,flss,srd,lmx,flkmxp,kmxp)

      !       Inverser le resultat avant de poursuivre les calculs.

      IFIN=flkmxp*.5 + 1

      do I=1,IFIN
         KI=flkmxp-I+1

         do J=1,LMX
            FUP(J,I)=flFUP(J,KI)
            FUP(J,KI)=flFUP(J,I)
            FDOWN(J,I)=flFDOWN(J,KI)
            FDOWN(J,KI)=flFDOWN(J,I)
         enddo
      enddo

   endif

   !     SMOOTHING OF FLUXES If RADFLTR IS TRUE

   do K=1,flkmxp
      do I=1,ILG
         RMUE(I,K)=FUP(I,K)
         RAY(I,K)=FDOWN(I,K)
      enddo
   enddo


   if ( RADFLTR ) then
      do K=2,flKMX
         do I=1,ILG
            FUP(I,K)=0.25*(RMUE(I,K+1)+RMUE(I,K-1)) + 0.5*RMUE(I,K)
            FDOWN(I,K) = 0.25*(RAY(I,K+1)+RAY(I,K-1)) + 0.5*RAY(I,K)
         enddo
      enddo
   else
      do K=2,flKMX
         do I=1,ILG
            FUP(I,K)= RMUE(I,K)
            FDOWN(I,K) = RAY(I,K)
         enddo
      enddo
   endif


   !     Calcul des taux (HEAT)

   do K=1,flKMX
      L=flkmxp-K
      LP1=L+1

      TOIT(L) = L.eq.flkmxp-1

      do I=1,ILG

         XNET=FUP(I,L)-FDOWN(I,L) - (FUP(I,LP1)-FDOWN(I,LP1))
         RE=RMUE(I,L)-RAY(I,L)-(RMUE(I,LP1)-RAY(I,LP1))
         if ( L.eq.1) XNET= RE

         !          AU NIVEAU DE LA SURFACE ON UTILISE LES FLUX SANS SMOOTHING

         HEAT(I,L)=max(CONS*XNET/(PSOL(I)*flDSIG(I,L)),0.)

      enddo

   enddo

   if(RADFIX) then
      call SUN_RADFIX1(flSIG,RMUO,RMUO,PSOL,HEAT,TOIT,ILG,flkmx)
   endif

   !     INVERT BACK most INPUTS
   KFIN=KMX*.5
   do K=1,KFIN
      KI=KMX-K+1
      do J=1,LMX
         A=DZ(J,K)
         DZ(J,K)=DZ(J,KI)
         DZ(J,KI)=A
         A=WV(J,K)
         WV(J,K)=WV(J,KI)
         WV(J,KI)=A
         A=TM(J,K)
         TM(J,K)=TM(J,KI)
         TM(J,KI)=A
         A=QOF(J,K)
         QOF(J,K)=QOF(J,KI)
         QOF(J,KI)=A
         A=CNEB(J,K)
         CNEB(J,K)=CNEB(J,KI)
         CNEB(J,KI)=A

      enddo
   enddo

   do I=1,KFIN
      do J=1,LMX
         IK=KMX-I+1
         A=DSIG(J,I)
         DSIG(J,I)=DSIG(J,IK)
         DSIG(J,IK)=A
         A=DSH(J,I)
         DSH(J,I)=DSH(J,IK)
         DSH(J,IK)=A
         A=DSC(J,I)
         DSC(J,I)=DSC(J,IK)
         DSC(J,IK)=A
      enddo
   enddo

   IFIN=flKMXP*.5

   do I=1,IFIN
      KI=flKMXP-I+1

      do J=1,LMX
         A=FUP(J,I)
         FUP(J,I)=FUP(J,KI)
         FUP(J,KI)=A
         A=FDOWN(J,I)
         FDOWN(J,I)=FDOWN(J,KI)
         FDOWN(J,KI)=A
      enddo
   enddo

   ifin=flkmx*.5

   do i=1,ifin
      ki=flkmx-i+1

      do j=1,lmx
         a = heat(j,i)
         heat(j,i) = heat(j,ki)
         heat(j,ki) = a
      enddo
   enddo

   return
end subroutine sun7b
