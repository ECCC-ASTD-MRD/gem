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

subroutine DTFRZNEW2(TU,P,THTEU,QVAP,QLIQ,QICE,RATIO2, &
     QNWFRZ,RL,FRC1,EFFQ,IFLAG,XLV0,XLV1,XLS0,XLS1, &
     ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE)
   implicit none
!!!#include <arch_specific.hf>

   integer IFLAG

   real A,B
   real C,CP,C5
   real DQVAP,DTFRZ2
   real EFFQ,ES,ESICE,ESLIQ,FRC1
   real P,PI,QICE,QLIQ,QLQFRZ,QNEW
   real QNWFRZ,QVAP,QVAP1,RATIO2,RL,RLC,RLF,RLS,RV
   real THTEU,TU,TU1
   real XLS0,XLS1,XLV0,XLV1
   real ALIQ,BLIQ,CLIQ,DLIQ,AICE,BICE,CICE,DICE

   !@AUTHOR Kain and Fritsch (Sept 25, 1986)
   !@REVISION
   ! 001      Stephane Belair (1994)
   !@OBJECT
   !          allow glaciation of the updraft to occur as an
   !          approximately linear function of temperature in the
   !          temperature range TTFRZ to TBFRZ...
   !@ARGUMENTS
   !          - INPUT/OUTPUT -
   ! TU       temperature of the updraft
   !          - INPUT -
   ! P        pressure of the level
   !          - OUTPUT -
   ! THTEU    potential temperature of the updraft
   !          - INPUT/OUTPUT -
   ! QVAP     water vapour content
   ! QLIQ     liquid water content
   ! QICE     cloud ice content
   !          - OUTPUT -
   ! RATIO2   degree of glaciation=(ESLIQ -ES)/(ESLIQ-ESICE)
   !          - INPUT/OUTPUT -
   ! QNWFRZ   fresh frozen condensate
   !          - OUPUT -
   ! RL       latent heat of vaporization
   !          - INPUT/OUTPUT -
   ! FRC1     fractional conversion to glaciation
   !          - INPUT -
   ! EFFQ     efficiency for the transformation of liquid water
   !          into cloud ice = (TTFRZ-TBFRZ)/(TTEMP-TBFRZ)
   !          = 0 ==> at the end of the glaciation zone
   !          = 1 ==> at the beginning of the glaciation zone
   !          - INPUT/OUTPUT -
   ! IFLAG    flag to determine if this layer is the first one
   !          of the updraft into the glaciation zone
   !          - Input -
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


   RV=461.51
   C5=1.0723E-3



   !                            Adjust the liquid water concentrations from
   !                            fresh condensate and that brought up from
   !                            lower levels to an amount that would be present
   !                            if no liquid water has frozen thus far.  This
   !                            is necessary because the expression for
   !                            temperature is multiplied by the fraction equal
   !                            to the parcel temperature decrease since the
   !                            last model level divided by the total glaciation
   !                            interval, so that effectively this approximatly
   !                            allows an amount of liquid water to freeze
   !                            which is equal to this same fraction of the
   !                            liquid water that was present before the
   !                            glaciation process was initiated.  Also, to
   !                            allow THETAU to convert approximately linearly
   !                            its value with respect to ice, we need to allow
   !                            a portion of the free condensate to contribute
   !                            to the glaciation process; the fractional
   !                            amount that applies to this portion is 1/2
   !                            of the fractional amount frozen of the "old"
   !                            condensate because this fresh condensate is
   !                            only produced gradually over the layer.
   !                            Note that in terms of the dynamics of the
   !                            precipitation process, i.e., precipitation
   !                            fallout, this fractional amount of fresh
   !                            condensate has already been included in the
   !                            ice category.

   QLQFRZ=QLIQ*EFFQ
   QNEW=QNWFRZ*EFFQ*0.5
   ESLIQ=ALIQ*exp((BLIQ*TU-CLIQ)/(TU-DLIQ))
   ESICE=AICE*exp((BICE*TU-CICE)/(TU-DICE))
   RLC=2.5E6-2369.276*(TU-273.16)
   RLS=2833922.-259.532*(TU-273.16)
   RLF=RLS-RLC
   CP=1005.7*(1.+0.89*QVAP)


   !                            A = D(es)/DT is that calculated from
   !                            Buck's (1981) empirical formulas
   !                            for saturation vapor pressure.

   A=(CICE-BICE*DICE)/((TU-DICE)*(TU-DICE))
   B=RLS*0.622/P
   C=A*B*ESICE/CP

   !                            Adjustment for water vapor and temperature.

   DQVAP=B*(ESLIQ-ESICE)/(RLS+RLS*C)-RLF*(QLQFRZ+QNEW)/(RLS+RLS/C)
   DTFRZ2=(RLF*(QLQFRZ+QNEW)+B*(ESLIQ-ESICE))/(CP+A*B*ESICE)

   !                            New temperature and water vapor.

   TU1=TU
   QVAP1=QVAP
   TU=TU+FRC1*DTFRZ2
   QVAP=QVAP-FRC1*DQVAP

   !                            Calculate RATIO2 (degree of glaciation).

   ES=QVAP*P/(0.622+QVAP)
   ESLIQ=ALIQ*exp((BLIQ*TU-CLIQ)/(TU-DLIQ))
   ESICE=AICE*exp((BICE*TU-CICE)/(TU-DICE))
   RATIO2=(ESLIQ-ES)/(ESLIQ-ESICE)


   !                            Typically, RATIO2 is very close to
   !                            (TTFRZ-TU)/(TTFRZ-TBFRZ).  Usually
   !                            within 1% (using TU before glaciation
   !                            effects are applied.  If the initial
   !                            temperature is below TBFRZ and RATIO2
   !                            is still less than 1, an adjustment to
   !                            FRC1 and RATIO2 is introduced so that
   !                            glaciation effects are not underestimated.
   !                            Conversely, if RATIO2 is greater than 1.
   !                            FRC1 is adjusted so that glaciation effects
   !                            are not overestimated.

   if(IFLAG.gt.0.and.RATIO2.lt.1)then
      FRC1=FRC1+(1.-RATIO2)
      TU=TU1+FRC1*DTFRZ2
      QVAP=QVAP1-FRC1*DQVAP
      RATIO2=1.
      IFLAG=1
      goto 20
   endif
   if(RATIO2.gt.1.)then
      FRC1=FRC1-(RATIO2-1.)
      FRC1=AMAX1(0.0,FRC1)
      TU=TU1+FRC1*DTFRZ2
      QVAP=QVAP1-FRC1*DQVAP
      RATIO2=1.
      IFLAG=1
   endif


   !                           Calculate a hybri value of THETAU, assuming
   !                           that the latent heat of vaporization/sublimation
   !                           can be estimated using the same weighting
   !                           function as that used to calculate saturation
   !                           vapor pressure.  Calculate new liquid water
   !                           and ice concentrations.

20 RLC=XLV0-XLV1*TU
   RLS=XLS0-XLS1*TU
   RL=RATIO2*RLS+(1.-RATIO2)*RLC
   PI=(1.E5/P)**(0.2854*(1.-0.28*QVAP))
   THTEU=TU*PI*exp(RL*QVAP*C5/TU*(1.+0.81*QVAP))

   if(IFLAG.eq.1)then

      !                           Above the glaciation zone, everything is
      !                           frozen, only have ice.

      QICE=QICE+FRC1*DQVAP+QLIQ
      QLIQ=0.
   else

      !                           In the transition zone.  Reduce QLIQ adn
      !                           increase QICE.  Remember also that water
      !                           vapor has to decrease due to the presence
      !                           of new ice.

      QICE=QICE+FRC1*(DQVAP+QLQFRZ)
      QLIQ=QLIQ-FRC1*QLQFRZ
   endif
   QNWFRZ=0.

   return
end subroutine DTFRZNEW2
