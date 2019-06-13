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
!** S/P PROF5
      SUBROUTINE PROF5(EQ,EE,UD)
!
!
      implicit none
!!!#include <arch_specific.hf>
!
!
      REAL A1,A2,A3,C1,C2,EE,EQ,EY,E45,FE,P,SIGMA,SQRT2P,T1,T2
      REAL UD,X,Y
!
!AUTHOR
!          Jack Kain and JM Fritsch (Oct 14,1990)
!
!REVISION
! 001      Stephane Belair (1994)
!
!OBJECT
!          integrates the area under the curve in the Gaussian
!          Distribution to determine the fractional entrainment
!          and detrainment rates given EQ, the critical fraction
!          of environmental air for neutral buoyancy.
!
!ARGUMENTS
!
!          - Input -
! EQ       is the fraction between 0 and 1 given to the equation
!
!          - Output -
! EE       fractional entrainment rate
! UD       fraction detrainment rate
!
!NOTES
!
!          The numerical approximation to the integral is taken
!          from "Handbook of Mathematical Functions with formulas,
!          graphs and mathematical tables" ED. by Abramowitz
!          and Stegun, Nat'l Bureau of Standards Applied Mathematics
!          Series. June, 1964., May, 1968.  Jack Kain
!
!*
!
      DATA SQRT2P,A1,A2,A3,P,SIGMA,FE/2.506628,0.4361836,-0.1201676, &
      0.9372980,0.33267,0.166666667,0.202765151/
!
      X=(EQ-0.5)/SIGMA
      Y=6.*EQ-3.
      EY=EXP(Y*Y/(-2))
      E45=EXP(-4.5)
      T2=1./(1.+P*ABS(Y))
      T1=0.500498
      C1=A1*T1+A2*T1*T1+A3*T1*T1*T1
      C2=A1*T2+A2*T2*T2+A3*T2*T2*T2
!
      IF(Y.GE.0.)THEN
        EE=SIGMA*(0.5*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY)) &
           -E45*EQ*EQ/2.
        UD=SIGMA*(0.5*(EY*C2-E45*C1)+SIGMA*(E45-EY)) &
           -E45*(0.5+EQ*EQ/2.-EQ)
      ELSE
        EE=SIGMA*(0.5*(EY*C2-E45*C1)+SIGMA*(E45-EY)) &
           -E45*EQ*EQ/2.
        UD=SIGMA*(0.5*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY)) &
           -E45*(0.5+EQ*EQ/2.-EQ)
      ENDIF
!
      EE=EE/FE
      UD=UD/FE
!
!
      RETURN
      END
