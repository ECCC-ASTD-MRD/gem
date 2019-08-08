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
!**s/r mfotvht  -  calcule la temperature virtuelle a partir de la
!                  temperature, de l'humidite specifique et la masse
!                  specifique des hydrometeores.
!                  note : la temperature virtuelle est celle qui tient
!                         compte de la vapeur et des hydrometeores.
!
      Subroutine mfotvht(tv,tt,qq,qh,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
      Integer ni, nk, n
      Real tv(ni,nk), tt(ni,nk), qq(ni,nk), qh(ni,nk)
!
!Author
!          B. Bilodeau (Sep 2002), based on mfotvt from N. Brunet  (Jan91)
!
!Object
!          To calculate virtual temperature tv from temperature tt,
!          specific humidity qq and specific mass of hydrometeors qh.
!          Note: The virtual temperature here is the one that accounts
!                for the vapor and the hydrometeors.
!
!Arguments
!
!          - Output -
! tv       virtual temperature in K
!
!          - Input -
! tt       temperature in K
! qq       specific humidity in kg/kg
! qv       specific mass of hydrometeors in kg/kg
! ni       horizontal dimension
! nk       vertical dimension
! n        number of points to process
!*
!--------------------------------------------------------------------
      Integer i, k
!--------------------------------------------------------------------

      Do k= 1,nk
         Do i=1,n
            tv(i,k) = fotvht(tt(i,k),qq(i,k),qh(i,k))
         Enddo
      Enddo

      End Subroutine mfotvht
