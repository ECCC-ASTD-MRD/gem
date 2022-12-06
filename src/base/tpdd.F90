!-------------------------------------- LICENCE BEGIN -------------------------
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

module tpdd

   private
   public :: tpdd1

contains

   function tpdd1(P,THTED,TGS,RS,RD,RH,XLV0,XLV1, &
        ALIQ,BLIQ,CLIQ,DLIQ)
      implicit none
!!!#include <arch_specific.hf>

      integer ITCNT

      real ALIQ,BLIQ,CLIQ,DLIQ,CP,DT
      real ES,F0,F1,RS,P,PI,T0,T1,THTGS,THTED,TGS
      real RL,XLV0,XLV1,RH,DQSDT,T1RH,RSRH,RD,TPDD1

      !@Author Jack Kain and JM Fritcsh (Oct 14,1990)
      !@Revision
      ! 001     Stephane Belair (1994)
      ! 002     Stephane Belair (2001) Impose minimum value for humidity
      !@Object
      !          This function iteratively extracts temperature from
      !          equivalent potential temperature. It is designed for use
      !          with downdraft calculations. If relative humidity is specified
      !          to be less than 100%, parcel temperature, specific humidity,
      !          and liquid water content are iteratively calculated.
      !@Arguments
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

      ES=ALIQ*exp((BLIQ*TGS-CLIQ)/(TGS-DLIQ))
      ES=min( ES , 0.5*P )
      RS=0.622*ES/(P-ES)
      RS=min( RS , 0.050 )
      RS=max( RS , 1.E-6 )
      PI=(1.E5/P)**(0.2854*(1.-0.28*RS))
      THTGS=TGS*PI*exp((3374.6525/TGS-2.5403)*RS* &
           (1.+0.81*RS))
      F0=THTGS-THTED
      T1=TGS-0.5*F0
      T0=TGS
      CP=1005.7

      !...ITERATE TO FIND WET-BULB TEMPERATURE...

      ITCNT=0
90    ES=ALIQ*exp((BLIQ*T1-CLIQ)/(T1-DLIQ))
      ES=min( ES , 0.5*P )
      RS=0.622*ES/(P-ES)
      RS=min( RS , 0.050 )
      RS=max( RS , 1.E-6 )
      PI=(1.E5/P)**(0.2854*(1.-0.28*RS))
      THTGS=T1*PI*exp((3374.6525/T1-2.5403)*RS* &
           (1.+0.81*RS))
      F1=THTGS-THTED
      if(abs(F1).lt.0.05)goto 50
      ITCNT=ITCNT+1
      if(ITCNT.gt.10)goto 50
      DT=F1*(T1-T0)/(F1-F0)
      T0=T1
      F0=F1
      T1=T1-DT
      goto 90
50    RL=XLV0-XLV1*T1

      !...IF RELATIVE HUMIDITY IS SPECIFIED TO BE LESS THAN 100%, ESTIMATE THE
      !   TEMPERATURE AND MIXING RATIO WHICH WILL YIELD THE APPROPRIATE VALUE.

      if(RH.eq.1.)goto 110
      DQSDT=(CLIQ-BLIQ*DLIQ)/((T1-DLIQ)*(T1-DLIQ))
      DT=RL*RS*(1.-RH)/(CP+RL*RH*RS*DQSDT)
      T1RH=T1+DT
      ES=RH*ALIQ*exp((BLIQ*T1RH-CLIQ)/(T1RH-DLIQ))
      ES=min( ES , 0.5*P )
      RSRH=0.622*ES/(P-ES)
      RSRH=min( RSRH , 0.050 )
      RSRH=max( RSRH , 1.E-6 )

      !...CHECK TO SEE IF MIXING RATIO AT SPECIFIED RH IS LESS THAN ACTUAL
      !...MIXING RATIO...IF SO, ADJUST TO GIVE ZERO EVAPORATION...

      if(RSRH.lt.RD)then
         RSRH=RD
         T1RH=T1+(RS-RSRH)*RL/CP
      endif
      T1=T1RH
      RS=RSRH
110   TPDD1=T1
      if(ITCNT.gt.10)print *,'***** NUMBER OF ITERATIONS IN TPDD = ', &
           ITCNT
      return
   end function tpdd1

end module tpdd
