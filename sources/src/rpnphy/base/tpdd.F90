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
!** FUNCTION TPDD
      FUNCTION TPDD(P,THTED,TGS,RS,RD,RH,XLV0,XLV1, &
         ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
!
      implicit none
#include <arch_specific.hf>
!
      INTEGER ITCNT
!
      REAL ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE,CP,DT
      REAL ES,F0,F1,RS,P,PI,T0,T1,THTGS,THTED,TGS
      REAL RL,XLV0,XLV1,RH,DQSDT,T1RH,RSRH,RD,TPDD
!
!Author
!          Jack Kain and JM Fritcsh (Oct 14,1990)
!
!Revision
! 001     Stephane Belair (1994)
! 002     Stephane Belair (2001) Impose minimum value for humidity
!
!Object
!          This function iteratively extracts temperature from
!          equivalent potential temperature. It is designed for use
!          with downdraft calculations. If relative humidity is specified
!          to be less than 100%, parcel temperature, specific humidity,
!          and liquid water content are iteratively calculated.
!
!Arguments
!
!          - Input -
! P        pressure
! THTED    environmental equivalent potential temperature
! TGS      temperature
!          - Input -
! RS       mixing ratio at specified relative humidity
! RD       actual mixing ratio
! RH       specified relative humidity
! XLV0     =3.147E+6 (constant for calculation of latent heating)
! XLV1     =2369. (constant for calculation of latent heating)
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
      ES=ALIQ*EXP((BLIQ*TGS-CLIQ)/(TGS-DLIQ))
      ES=MIN( ES , 0.5*P )
      RS=0.622*ES/(P-ES)
      RS=MIN( RS , 0.050 )
      RS=MAX( RS , 1.E-6 )
      PI=(1.E5/P)**(0.2854*(1.-0.28*RS))
      THTGS=TGS*PI*EXP((3374.6525/TGS-2.5403)*RS* &
           (1.+0.81*RS))
      F0=THTGS-THTED
      T1=TGS-0.5*F0
      T0=TGS
      CP=1005.7
!
!...ITERATE TO FIND WET-BULB TEMPERATURE...
!
      ITCNT=0
90    ES=ALIQ*EXP((BLIQ*T1-CLIQ)/(T1-DLIQ))
      ES=MIN( ES , 0.5*P )
      RS=0.622*ES/(P-ES)
      RS=MIN( RS , 0.050 )
      RS=MAX( RS , 1.E-6 )
      PI=(1.E5/P)**(0.2854*(1.-0.28*RS))
      THTGS=T1*PI*EXP((3374.6525/T1-2.5403)*RS* &
           (1.+0.81*RS))
      F1=THTGS-THTED
      IF(ABS(F1).LT.0.05)GOTO 50
      ITCNT=ITCNT+1
      IF(ITCNT.GT.10)GOTO 50
      DT=F1*(T1-T0)/(F1-F0)
      T0=T1
      F0=F1
      T1=T1-DT
      GOTO 90
50    RL=XLV0-XLV1*T1
!
!...IF RELATIVE HUMIDITY IS SPECIFIED TO BE LESS THAN 100%, ESTIMATE THE
!   TEMPERATURE AND MIXING RATIO WHICH WILL YIELD THE APPROPRIATE VALUE.
!
      IF(RH.EQ.1.)GOTO 110
        DQSDT=(CLIQ-BLIQ*DLIQ)/((T1-DLIQ)*(T1-DLIQ))
        DT=RL*RS*(1.-RH)/(CP+RL*RH*RS*DQSDT)
        T1RH=T1+DT
        ES=RH*ALIQ*EXP((BLIQ*T1RH-CLIQ)/(T1RH-DLIQ))
        ES=MIN( ES , 0.5*P )
        RSRH=0.622*ES/(P-ES)
        RSRH=MIN( RSRH , 0.050 )
        RSRH=MAX( RSRH , 1.E-6 )
!
!...CHECK TO SEE IF MIXING RATIO AT SPECIFIED RH IS LESS THAN ACTUAL
!...MIXING RATIO...IF SO, ADJUST TO GIVE ZERO EVAPORATION...
!
        IF(RSRH.LT.RD)THEN
          RSRH=RD
          T1RH=T1+(RS-RSRH)*RL/CP
        ENDIF
        T1=T1RH
        RS=RSRH
110   TPDD=T1
      IF(ITCNT.GT.10)PRINT *,'***** NUMBER OF ITERATIONS IN TPDD = ', &
                              ITCNT
      RETURN
      END
