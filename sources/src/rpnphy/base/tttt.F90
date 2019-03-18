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
!**S/R TTTT
!
       SUBROUTINE TTTT(JJ,W,R1,R2,LMX,ILG,APAD,BPAD,D,AKI)
      implicit none
#include <arch_specific.hf>
      INTEGER I,II,ILG,JJ,KI,LMX
      REAL APAD(3,6),BPAD(3,6),D(3),AKI(2)
      REAL      W(LMX),R1(LMX),R2(LMX)
!
!Author
!          L.Garand (1989)
!
!Revision
! 001      G.Pellerin(Mar90)Standard documentation
!
!Object
!          to calculate transmission functions by Pade approximants
!          using Horner's algorithm (vectorized version)
!
!Arguments
!
!          - Input -
! JJ       =1, H2O for U between 1.E-05 and 100 g.cm-2
!          =2, CO2 for U between .3 and 10 g.cm-2
!          =3, O3 for U between .01 and 10 g.cm-2
! W        effective absorber amount
!
!          - Output -
! R1       result of transmission function 1
! R2       result of transmission function 2
!
!          - Input -
! LMX      maximum number of points in the function
! ILG      number of points in the function
! APAD     coefficients for Pade approximants
! BPAD     coefficients for Pade approximants
! D        coefficients for Pade approximants
! AKI      empirical grey absorption coefficients
!
!*
!
!
      DO 1 I=1,ILG
      R1(I)=APAD(JJ,6)
      R2(I)=BPAD(JJ,6)
 1    CONTINUE
      DO 3 II=1,5
      KI=6-II
      DO 2 I=1,ILG
      R1(I)=R1(I)*W(I)+APAD(JJ,KI)
      R2(I)=R2(I)*W(I)+BPAD(JJ,KI)
 2    CONTINUE
 3    CONTINUE
!
      DO 4 I=1,ILG
      R1(I)=(R1(I)/R2(I))*(1-D(JJ))+D(JJ)
 4    CONTINUE
!
      RETURN
      END
