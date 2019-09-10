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
!**s/r mhrahu3  -  passage de hr a hu
!
      Subroutine mhrahu3(hu,hr,tt,ps,swph,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
!
      Integer ni, nk, n
      Real hu(ni,nk), hr(ni,nk), tt(ni,nk)
      Real ps(ni,*)
!
      Logical swph
!
!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      N. Brunet (May 1994) - Improved iteration for tx=tv
! 003      B. Bilodeau (January 2001) - Automatic arrays
!
!Object
!          to calculate the specific humidity from relative humidity,
!          temperature and pressure
!
!Arguments
!
!          - Output -
! hu       specific humidity in kg/kg
!
!          - input -
! hr       relative humidity in fraction
! tt       temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! ni       horizontal dimension
! nk       vertical dimension
! n        number of treated points
!*
!--------------------------------------------------------------------
      Real e
      Integer k, i
!--------------------------------------------------------------------
      If(swph)Then
         Do k=1,nk
         Do i=1,n
            e = dmin1(Dble(ps(i,k)),hr(i,k)*foew(tt(i,k)))
            hu(i,k) = foqfe(e,ps(i,k))
         Enddo
         Enddo
     Else
         Do k=1,nk
         Do i=1,n
            e = dmin1(Dble(ps(i,k)),hr(i,k)*foewa(tt(i,k)))
            hu(i,k) = foqfe(e,ps(i,k))
         Enddo
         Enddo
     End If
!
     End Subroutine mhrahu3
