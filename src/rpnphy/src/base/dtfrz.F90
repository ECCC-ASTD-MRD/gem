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
!** FUNCTION DTFRZ
      FUNCTION DTFRZ (T,P,Q,SD,D273)
      implicit none
!!!#include <arch_specific.hf>
      REAL DTFRZ,T,P,Q,SD,D273
!
!Author
!          Fritsch and Chappell  (Sept 25, 1986)
!
!Revision
!
!Object
!     to calculate the difference of T due to freezing
!
!Arguments
!
!          - Input -
! T        temperature
! P        pressure
! Q        specific humidity
! SD       difference of moisture
! D273     = 1/273.16
!
!Notes
!
!
!*
!
      REAL RLC,RLS,TOTRU,BT,DT,QS1,QS2,QS,X
!
!
!...COMPUTE FREEZING EFFECT.
!
      RLC = 2.5E6-2369.276*(T-273.16)
      RLS = 2833922.-259.532*(T-273.16)
      TOTRU = SD+Q
      BT = 1./T
      DT = D273-BT
      QS1 = 611.*EXP(RLC*DT/461.)
      QS2 = 611.*EXP(RLS*DT/461.)
      QS = .622*RLS*QS2/(P-QS2)
      X = 2.1672E-3*RLS*QS*BT*BT
      DTFRZ = (SD*(RLS-RLC)+.622*RLS*QS1/(P-QS1)-QS)/(1004.+1952.*TOTRU+ &
         X)
      RETURN
      END
