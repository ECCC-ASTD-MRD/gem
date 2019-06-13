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
!**s/r mfottv  -  calcule temp. a partir de tv et hum. sp.
!
      Subroutine mfottv(tt,tv,qq,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
      Integer ni, nk, n
      Real tt(ni,nk), tv(ni,nk), qq(ni,nk)
!
!Author
!          N. Brunet  (Jan91)
!
!Object
!          to calculate temperature tt from virtual temperature tv
!          and specific humidity qq
!
!Arguments
!
!          - Output -
! tt       temperature in K
!
!          - Input -
! tv       virtual temperature in K
! qq       specific humidity in kg/kg
! ni       horizontal dimension
! nk       vertical dimension
! n        number of points to process
!*
!--------------------------------------------------------------------
      Integer k, i
!--------------------------------------------------------------------
!
      Do k=1,nk
      Do i=1,n
         tt(i,k) = fottv(tv(i,k),qq(i,k))
      Enddo
      Enddo
!
      End Subroutine mfottv
