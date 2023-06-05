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
!**FUNCTION TP
      FUNCTION TP (PRESS,THETAE,TGS,D273,RL,QS,PI)
      implicit none
!!!#include <arch_specific.hf>
      REAL TP,PRESS,THETAE,TGS,D273,RL,QS,PI
!
!Author
!          Da-Lin Zhang     (Sept 25, 1986)
!
!Revision
! 001      Stephane Belair (Oct 2003) - Add clips to secure the code
!
!Object
!          to convert from equivalent potential temperature
!          to temperature
!
!Arguments
!
!          - Input -
! PRESS    pressure
! THETAE   equivalent potential temperature of environment
! TGS      input temperature
! D273     = 1/273.16
! RL       latent heat of vaporization
! QS       specific humidity in updraft
! PI       conversion factor from T to Theta
!
!Notes
!
!*
!
      INTEGER COUNT
      REAL RL461,RL1004,RP,ES,FO,T1,TGUESS,F1,DT
!
!
!...ITERATIVELY CONVERT TEMPERATURE FROM EQUIVALENT POTENTIAL
!...TEMPERATURE.
!
      COUNT = 0
      RL461 = RL/461.
      RL1004 = RL/1004.
      RP = THETAE/PI
      ES = 611.*EXP(RL461*(D273-1./TGS))
      ES = MIN(ES,0.5*PRESS)
      QS = .622*ES/(PRESS-ES+1.E-10)
      QS = MAX(QS,0.)
      QS = MIN(QS,0.050)
      FO = TGS*EXP(RL1004*QS/TGS)-RP
      T1 = TGS-.5*FO
      T1 = AMAX1( T1, 150.0 )
      TGUESS = TGS
   10 ES = 611.*EXP(RL461*(D273-1./T1))
      ES = MIN(ES,0.5*PRESS)
      QS = .622*ES/(PRESS-ES+1.E-10)
      QS = MAX(QS,0.)
      QS = MIN(QS,0.050)
      F1 = T1*EXP(RL1004*QS/T1)-RP
      IF (ABS(F1).LT..05) GOTO 20
      DT = F1*(T1-TGUESS)/(F1-FO+1.E-10)
      TGUESS = T1
      FO = F1
      T1 = T1-DT
      T1 = AMAX1( T1, 150.0 )
      COUNT = COUNT + 1
      IF (COUNT.GE.15) THEN
        TP = TGS
        RETURN
      ENDIF
      GOTO 10
!
   20 TP = T1
      RETURN
      END
