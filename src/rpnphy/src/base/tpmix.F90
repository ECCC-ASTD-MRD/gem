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
!** S/P TPMIX
       SUBROUTINE TPMIX(P,THTU,TU,QU,QLIQ,QICE,QNEWLQ,QNEWIC,RATIO2,RL, &
                         XLV0,XLV1,XLS0,XLS1, &
                         ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!
!
      implicit none
!!!#include <arch_specific.hf>
!
      INTEGER ITCNT
!
      REAL ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE
      REAL CP,C5
      REAL DQ,DQICE,DQLIQ,DT
      REAL ES,ESICE,ESLIQ,F0,F1,F2,P,PPI
      REAL QICE,QLIQ,QNEW,QNEWIC,QNEWLQ,QS,QTOT
      REAL QU,RATIO2,RL,RLL,RV,THTGS,THTU,TU,T0,T1
      REAL XLS0,XLS1,XLV0,XLV1
!
!
!AUTHOR
!          Jack Kain and JM Fritcsh (Oct 14,1990)
!
!REVISION
! 001      Stephane Belair (1994)
! 002      Richard Moffet (1999) - Set limits to QS,F1-F0.
!                                  Add new thermodynamic functions.
! 003      N. Brunet and S. Belair (Oct 2000) - Remove new thermodynamic
!          functions and constants in order to have consistency with the
!          rest of Kain-Fritsch code.
! 004      Stephane Belair (2001) - Impose minimum values for humidity
!
!OBJECT
!          This subroutine iteratively extracts wet-bulb temperature
!          from equivalent potential temperature, then checks to see
!          if sufficient moisture is available to achieve saturation.
!          If not, temperature is adjusted accordingly, if so, the
!          residual liquid water/ice concentration is determined.
!
!Arguments
!
!          - Input -
! P        pressure
! THTU     equivalent potential temperature of the updraft
! TU       temperature of the updraft
! QU       specific humidity of the updraft
!
!          - Input/Output -
! QLIQ     liquid water in the updraft
! QICE     cloud ice in the updraft
!
!          - Output -
! QNEWLQ   liquid water after correction (evap/cond)
! QNEWIC   cloud water after correction
!
!          - Input -
! RATIO2   degree of glaciation
! RL       latent heat of vaporization
! XLV0     =3.147E+6 (constant for calculation of latent heating)
! XLV1     =2369. (constant for calculation of latent heating)
! XLS0     =2.905E+6(constant for calculation of latent heating)
! XLS1     =259.532(constant for calculation of latent heating)
! ALIQ     =613.3(constant for calcul. of saturation vapor pressure)
! BLIQ     =17.502(constant for calcul. of saturation vapor pressure)
! CLIQ     =4780.8(constant for calcul. of saturation vapor pressure)
! DLIQ     =32.19(constant for calcul. of saturation vapor pressure)
! AICE     =613.2(constant for calcul. of saturation vapor pressure)
! BICE     =22.452(constant for calcul. of saturation vapor pressure)
! CICE     =6133.0(constant for calcul. of saturation vapor pressure)
! DICE     =0.61(constant for calcul. of saturation vapor pressure)
!
!Notes
!
!*
!
!
      C5=1.0723E-3
      RV=461.51
!
!
!                         Iteration to find wet bulb temperature as a
!                         function of equivalent potential temperature
!                         and pressure, assuming saturation vapor pressure.
!                         RATIO2 is the degree of glaciation.
!
!                         First calculate the eq. pot. temperature based
!                         on the temperature of the updraft.
!
         IF(RATIO2.LT.1.E-6)THEN
!
!                         Below the glaciation zone...
!
           ES=ALIQ*EXP((BLIQ*TU-CLIQ)/(TU-DLIQ))
           ES=MIN( ES , 0.5*P )
           QS=0.622*ES/(P-ES)
           QS=MAX(QS,0.)
           QS=MIN(QS,0.050)
           PPI=(1.E5/P)**(0.2854*(1.-0.28*QS))
           THTGS=TU*PPI*EXP((3374.6525/TU-2.5403)*QS* &
           (1.+0.81*QS))
         ELSEIF(ABS(RATIO2-1.).LT.1.E-6)THEN
!
!                         Above the glaciation zone...
!
           ES=AICE*EXP((BICE*TU-CICE)/(TU-DICE))
           ES=MIN( ES , 0.5*P )
           QS=0.622*ES/(P-ES)
           QS=MAX(QS,0.)
           QS=MIN(QS,0.050)
           PPI=(1.E5/P)**(0.2854*(1.-0.28*QS))
           THTGS=TU*PPI*EXP((3114.834/TU-0.278296)*QS* &
           (1.+0.81*QS))
         ELSE
!
!                         In the glaciation zone...
!
           ESLIQ=ALIQ*EXP((BLIQ*TU-CLIQ)/(TU-DLIQ))
           ESICE=AICE*EXP((BICE*TU-CICE)/(TU-DICE))
           ES=(1.-RATIO2)*ESLIQ+RATIO2*ESICE
           ES=MIN( ES , 0.5*P )
           QS=0.622*ES/(P-ES)
           QS=MAX(QS,0.)
           QS=MIN(QS,0.050)
           PPI=(1.E5/P)**(0.2854*(1.-0.28*QS))
           THTGS=TU*PPI*EXP(RL*QS*C5/TU*(1.+0.81*QS))
         ENDIF
!
!
!                           First adjustment towards the updraft
!                           temperature.
!
!
      F0=THTGS-THTU
      T1=TU-0.5*F0
      T1=MAX(T1,150.)
      T0=TU
      ITCNT=0
!
!
!
!
90    IF(RATIO2.LT.1.E-6)THEN
!
!                          Below the glaciation zone...
!
         ES=ALIQ*EXP((BLIQ*T1-CLIQ)/(T1-DLIQ))
         ES=MIN( ES , 0.5*P )
         QS=0.622*ES/(P-ES)
         QS=MAX(QS,0.)
         QS=MIN(QS,0.050)
         PPI=(1.E5/P)**(0.2854*(1.-0.28*QS))
         THTGS=T1*PPI*EXP((3374.6525/T1-2.5403)*QS* &
         (1.+0.81*QS))
      ELSEIF(ABS(RATIO2-1.).LT.1.E-6)THEN
!
!                          Above the glaciation zone...
!
         ES=AICE*EXP((BICE*T1-CICE)/(T1-DICE))
         ES=MIN( ES , 0.5*P )
         QS=0.622*ES/(P-ES)
         QS=MAX(QS,0.)
         QS=MIN(QS,0.050)
         PPI=(1.E5/P)**(0.2854*(1.-0.28*QS))
         THTGS=T1*PPI*EXP((3114.834/T1-0.278296)*QS* &
         (1.+0.81*QS))
      ELSE
!
!                          In the glaciation zone...
!
         ESLIQ=ALIQ*EXP((BLIQ*T1-CLIQ)/(T1-DLIQ))
         ESICE=AICE*EXP((BICE*T1-CICE)/(T1-DICE))
         ES=(1.-RATIO2)*ESLIQ+RATIO2*ESICE
         ES=MIN( ES , 0.5*P )
         QS=0.622*ES/(P-ES)
         QS=MAX(QS,0.)
         QS=MIN(QS,0.050)
         PPI=(1.E5/P)**(0.2854*(1.-0.28*QS))
         THTGS=T1*PPI*EXP(RL*QS*C5/T1*(1.+0.81*QS))
      ENDIF
!
!
!                          Departure of the "guess" form the updraft
!                          eq. pot. temperature.
!
      F1=THTGS-THTU
!
!
!                          Iterative process has converged.
!
      IF(ABS(F1).LT.0.01)GOTO 50
!
!
!                          OR we need one more iteration.
!
      ITCNT=ITCNT+1
!
!
!                          Too many iterations already ... get out.
!
      IF(ITCNT.GT.10)GOTO 50
!
!
!                          New adjustment for the next iteration.
      F2=F1-F0
      if ( abs(F2).lt.0.05) GO TO 50
      DT=F1*(T1-T0)/(F1-F0)
      T0=T1
      F0=F1
      T1=T1-DT
      T1=MAX(T1,150.)
      GOTO 90
!
!
!                          ITERATION PROCESS IS OVER.
!
!
!
!                          If the parcel is supersaturated, calculate
!                          the fresh condensation...
!
50    IF(QS.LE.QU)THEN
        QNEW=QU-QS
        QU=QS
        GOTO 96
      ENDIF
!
!
!
!                          If the humidity of the parcel is less than
!                          saturation, temperature and mixing ratio must
!                          be adjusted.  If liquid water or ice is
!                          present, must evaporate or sublimate.
!
      QNEW=0.
      DQ=QS-QU
      QTOT=QLIQ+QICE
!
!
!
!                          If there is enough liquid or ice to saturate
!                          the parcel, the temperature stays at the wet
!                          bulb value, the vapor mixing ratio is at the
!                          saturated level, and the mixing ratios of
!                          liquid water and ice are adjusted to make up
!                          for the original saturation deficit.
!                          Otherwise, any available liquid water or ice
!                          vaporizes and appropriate adjustments to parcel
!                          temperature, vapor, liquid and ice mixing
!                          ratios are done.
!
!                          Note that the liquid water and ice may be
!                          present in proportions slightly different than
!                          suggested by the value of RATIO2.  Check to
!                          make sure that liquid water and ice concentrations
!                          are not reduced to values below zero when
!                          evaporation/sublimation occurs.
!
!                          First case:  enough liquid water annd ice to
!                          saturate the updraft.
!
      IF(QTOT.GE.DQ)THEN
        DQICE=0.0
        DQLIQ=0.0
        QLIQ=QLIQ-(1.-RATIO2)*DQ
        IF(QLIQ.LT.0.)THEN
          DQICE=0.0-QLIQ
          QLIQ=0.0
        ENDIF
        QICE=QICE-RATIO2*DQ+DQICE
        IF(QICE.LT.0.)THEN
          DQLIQ=0.0-QICE
          QICE=0.0
        ENDIF
        QLIQ=QLIQ+DQLIQ
        QU=QS
        GOTO 96
      ELSE
        IF(RATIO2.LT.1.E-6)THEN
          RLL=XLV0-XLV1*T1
        ELSEIF(ABS(RATIO2-1.).LT.1.E-6)THEN
          RLL=XLS0-XLS1*T1
        ELSE
          RLL=RL
        ENDIF
        CP=1005.7*(1.+0.89*QU)
        IF(QTOT.LT.1.E-10)THEN
!
!                         If no liquid water or ice is available,
!                         the temperature is given by:
!
          T1=T1+RLL*(DQ/(1.+DQ))/CP
          T1=MAX(T1,150.)
          GOTO 96
!
!
!                         Second case:  not enough water and ice to
!                         saturate the updraft.
!
        ELSE
!
!
!                         If some liquid water/ice is available, but
!                         not enough to achieve saturation, the
!                         temperature is given by:
!
          T1=T1+RLL*((DQ-QTOT)/(1+DQ-QTOT))/CP
          T1=MAX(T1,150.)
          QU=QU+QTOT
          QTOT=0.
        ENDIF
        QLIQ=0
        QICE=0.
      ENDIF
!
!
!                         Partition the new condensate into liquid water
!                         and ice categories.
!
96    TU=T1
      QNEWLQ=(1.-RATIO2)*QNEW
      QNEWIC=RATIO2*QNEW
!
!
      IF(ITCNT.GT.10)PRINT *,'***** NUMBER OF ITERATIONS IN TPMIX =', &
                              ITCNT
!
!
      RETURN
      END
