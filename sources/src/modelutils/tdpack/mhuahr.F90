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
!**s/r mhuahr3  -  passage de hu a hr
!
      Subroutine mhuahr3(hr,hu,tt,ps,swph,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
!
      Integer ni, nk, n
      Real hr(ni,nk), hu(ni,nk), tt(ni,nk)
      Real ps(ni,*)
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      B. Bilodeau (January 2001) - Automatic arrays
!
!Object
!          to calculate relative humidity from specific humidity,
!          temperature and pressure
!
!Arguments
!
!          - Output -
! hr       relative humidity in fraction
!
!          - Input -
! hu       specific humidity in kg/kg
! tt       temperature or virtual temperature in Kelvins
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! ni       horizontal dimension
! nk       vertical dimension
! n        number of treated points
!*
!--------------------------------------------------------------------
      Integer i, k
!--------------------------------------------------------------------
!
      If(swph)Then
         Do k=1,nk
         Do i=1,n
               hr(i,k) = fohr(hu(i,k),tt(i,k),ps(i,k))
         Enddo
         Enddo
      Else
         Do k=1,nk
         Do i=1,n
               hr(i,k) = fohra(hu(i,k),tt(i,k),ps(i,k))
         Enddo
         Enddo
      Endif
!
      End Subroutine mhuahr3
