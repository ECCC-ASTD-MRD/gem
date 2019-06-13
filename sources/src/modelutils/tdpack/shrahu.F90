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
!**fonction shrahu3  -  passage de hr a hu
!
      Function shrahu3(hr,tt,ps,swph)
      use tdpack, only: foew, foewa, foqfe
      implicit none
!!!#include <arch_specific.hf>
      Real shrahu3, hr, tt, ps
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Revision 01 - N. Brunet (May 1994) - Improved iteration for tx = tv
!
!Object
!          to return specific humidity(kg/kg) calculated from
!          relative humidity, temperature and pressure
!
!Arguments
!
!          - Input -
! hr       relative humidity (fraction)
! tt       temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
!*
!--------------------------------------------------------------------
      Real e
!--------------------------------------------------------------------
!
      If(swph)Then
         e = dmin1(Dble(ps),hr * foew(tt))
      Else
         e = dmin1(Dble(ps),hr * foewa(tt))
      End If
      shrahu3 = foqfe(e,ps)
!
      End Function shrahu3
