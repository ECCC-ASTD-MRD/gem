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
!**fonction shuahr3  -  passage de hu a hr
!
      Function shuahr3(hu,tt,ps,swph)
      use tdpack, only: fohr, fohra
      implicit none
!!!#include <arch_specific.hf>
      Real shuahr3, hu, tt, ps
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Object
!          to return relative humidity(fraction) calculated from
!          specific humidity, temperature and pressure
!
!Arguments
!
!          - Input -
! hu       specific humidity in kg/kg
! tx       temperature or virtual temperature in K
! ps       pressure in Pa
! swtt     .true. to pass tt for argument
!          .false. to pass tv for argument
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
!*
!--------------------------------------------------------------------

      If(swph)Then
         shuahr3 = fohr(hu,tt,ps)
      Else
         shuahr3 = fohra(hu,tt,ps)
      End If
!
      End Function shuahr3
