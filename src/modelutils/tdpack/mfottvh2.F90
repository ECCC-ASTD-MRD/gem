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
!**s/p mfottvh2  -  Calcule la temperature a partir de la temperature
!                   virtuelle, l'humidite specifique et la masse specifique
!                   des hydrometeores si tt2vt=.false.
!                   L'inverse si tt2vt=.true.
!                   Note : La temperature virtuelle est celle qui tient
!                          compte de la vapeur et des hydrometeores.
!
      Subroutine mfottvh2(tt,tv,qq,qh,minx,maxx,miny,maxy,nk,i0,in,j0,jn,tt2vt)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
!
      Logical tt2vt
      Integer minx,maxx,miny,maxy,nk,i0,in,j0,jn
      Real tt(minx:maxx,miny:maxy,nk), tv(minx:maxx,miny:maxy,nk), &
           qq(minx:maxx,miny:maxy,nk), qh(minx:maxx,miny:maxy,nk)
!
!author
!          A. Plante (Apr 2002), based on mfottv from N. Brunet  (Jan 91)
!
!Arguments
!
!          - Output -
! tt       temperature in K
!
!          - Input -
! tv       virtual temperature in K
! qq       specific humidity in kg/kg
! qv       specific mass of hydrometeors in kg/kg
!*
!--------------------------------------------------------------------
      Integer i,j,k
!--------------------------------------------------------------------

      If (tt2vt) Then
         Do k= 1,nk
         Do j=j0,jn
         Do i=i0,in
            tv(i,j,k) = fotvht(tt(i,j,k),qq(i,j,k),qh(i,j,k))
         Enddo
         Enddo
         Enddo
      Else
         Do k= 1,nk
         Do j=j0,jn
         Do i=i0,in
            tt(i,j,k) = fottvh(tv(i,j,k),qq(i,j,k),qh(i,j,k))
         Enddo
         Enddo
         Enddo
      Endif

      End Subroutine mfottvh2
