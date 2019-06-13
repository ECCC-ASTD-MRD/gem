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
!**s/r fonction shraes3  -  passage de hr a es
!
      Function shraes3(hr,tt,ps,swph)
      use tdpack, only: shrahu3, shuaes3
      implicit none
!!!#include <arch_specific.hf>
      Real shraes3, hr, tt, ps
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Object
!          to return the dew point depression(Celsius) calculated
!          from relative humidity, temperature and pressure
!
!Arguments
!
!          - Input -
! hr       relative humidity (fraction)
! tt       temperature or virtual temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
!*
!--------------------------------------------------------------------
      Real hus
!--------------------------------------------------------------------
!
      hus = shrahu3(hr,tt,ps,swph)
      shraes3 = shuaes3(hus,tt,ps,swph)
!
      End Function shraes3
