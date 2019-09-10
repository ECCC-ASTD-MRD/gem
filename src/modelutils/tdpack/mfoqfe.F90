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
!**s/r mfoqfe3  -  calcule hum. sp. de tension de vap. et pres.
!
      Subroutine mfoqfe3(qq,ee,ps,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
!
      Integer ni, nk, n
      Real qq(ni,nk), ee(ni,nk)
      Real ps(ni,*)
!
!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      B. Bilodeau (January 2001) - Automatic arrays
!
!Object
!          to calculate specific humidity from vapour pressure and
!          pressure
!
!Arguments
!
!          - Output -
! qq       specific humidity in kg/kg
!
!          - Input -
! ee       vapour pressure in Pa
! ps       pressure in Pa
! ni       horizontal dimension
! nk       vertical dimension
! n        number of points to process
!*
!--------------------------------------------------------------------
      Integer i, k
!--------------------------------------------------------------------
!
!
      Do k=1,nk
      Do i=1,n
         qq(i,k) = foqfe(ee(i,k),ps(i,k))
      Enddo
      Enddo
!
!
      End Subroutine mfoqfe3
