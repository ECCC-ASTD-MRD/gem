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
!** S/P ENVIRTHT
      SUBROUTINE ENVIRTHT(P1,T1,Q1,THT1,R1,RL, &
         ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!
!
      implicit none
!!!#include <arch_specific.hf>
!
!
      REAL ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE
      REAL C1,C2,C3,C4,C5
      REAL EE
      REAL P00,P1
      REAL Q1
      REAL RL,R1
      REAL TDPT,TFPT,TLOG,TLOGIC,THT
      REAL THT1,TSAT,TSATIC,TSATLQ,T00,T1
!
!AUTHOR
!          Kain and Fritsch
!
!REVISION
!
!OBJECT
!          to calculate environmental equivalent potential
!          temperature...
!
!ARGUMENTS
!
!          - Input -
! P1       pressure
! T1       temperature
! Q1       specific humidity
!
!          - Output -
! THT1     environmental equivalent potential temperature
!
!          - Input -
! R1       RATIO2(degree of glaciation)
! RL       latent heat of vaporization
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
      DATA T00,P00,C1,C2,C3,C4,C5/273.16,1.E5,3374.6525,2.5403,3114.834, &
           0.278296,1.0723E-3/
!
         IF(R1.LT.1.E-6)THEN
!
!                                 Lower than the glaciation zone
!
           EE=Q1*P1/(0.622+Q1)
           TLOG=ALOG(EE/ALIQ)
           TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
           TSAT=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(T1-T00))* &
                 (T1-TDPT)
           THT=T1*(P00/P1)**(0.2854*(1.-0.28*Q1))
           THT1=THT*EXP((C1/TSAT-C2)*Q1*(1.+0.81*Q1))
         ELSEIF(ABS(R1-1.).LT.1.E-6)THEN
!
!                                 Above the glaciation zone
!
           EE=Q1*P1/(0.622+Q1)
           TLOG=ALOG(EE/AICE)
           TFPT=(CICE-DICE*TLOG)/(BICE-TLOG)
           THT=T1*(P00/P1)**(0.2854*(1.-0.28*Q1))
           TSAT=TFPT-(.182+1.13E-3*(TFPT-T00)-3.58E-4*(T1-T00))* &
                 (T1-TFPT)
           THT1=THT*EXP((C3/TSAT-C4)*Q1*(1.+0.81*Q1))
         ELSE
!
!                                In the glaciation zone
!
           EE=Q1*P1/(0.622+Q1)
           TLOG=ALOG(EE/ALIQ)
           TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
           TLOGIC=ALOG(EE/AICE)
           TFPT=(CICE-DICE*TLOGIC)/(BICE-TLOGIC)
           THT=T1*(P00/P1)**(0.2854*(1.-0.28*Q1))
           TSATLQ=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(T1-T00))* &
                  (T1-TDPT)
           TSATIC=TFPT-(.182+1.13E-3*(TFPT-T00)-3.58E-4*(T1-T00))* &
                  (T1-TFPT)
           TSAT=R1*TSATIC+(1.-R1)*TSATLQ
           THT1=THT*EXP(RL*Q1*C5/TSAT*(1.+0.81*Q1))
         ENDIF
!
!
!
      RETURN
      END
