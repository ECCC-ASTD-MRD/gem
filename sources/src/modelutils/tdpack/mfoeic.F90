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
!**s/r mfoeic  -  calcule tension de vapeur saturante ei.
!              -  (glace seulement)
!
      Subroutine mfoeic(ei,tt,ni,nk,n)
      use tdpack
      implicit none
#include <arch_specific.hf>
      Integer ni, nk, n
      Real ei(ni,nk), tt(ni,nk)
!
!Author
!          A. Plante (June 2003) - based on MFOEWA from N. Brunet  (Jan91)
!
!Revision
!
!Object
!          to calculate the saturation vapour pressure. (Ice phase
!          considered only for all temperatures)
!
!Arguments
!
!          - Output -
! ei       saturated vapour pressure in Pa
!
!          - Input -
! tt       temperature in K
! ni       horizontal dimension
! nk       vertical dimension
! n        number of points to process
!--------------------------------------------------------------------
      Integer i, k
      Real*8, Dimension(ni,nk) :: work
!***********************************************************************
      Do k=1,nk
         Do i=1,n
            work(i,k)=fesif(tt(i,k))
         Enddo
      Enddo
      Call vexp(work,work,n*nk)
      Do k=1,nk
         Do i=1,n
            ei(i,k)=aerk1i*work(i,k)
         Enddo
      Enddo
!
      End Subroutine mfoeic
