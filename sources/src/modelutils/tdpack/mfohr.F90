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
!**s/r mfohr  - compute relative humidity from specific humidity
!               temperature and pression 

      Subroutine mfohr4 (hr, qq, tt, ps, ni, nk, n, satuco)
      use tdpack_const
      implicit none
#include <arch_specific.hf>

      logical satuco
      Integer ni, nk, n
      Real hr(ni,*), qq(ni,*), tt(ni,*)
      Real ps(ni,*)

!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      B. Bilodeau  (January 2001) - Automatic arrays
! 003      G. Pellerin  (June 2003) - IBM conversion
!                  - calls to vexp routine (from massvp4 library)
!
!Object
!          to calculate relative humidity from specific humidity,
!          temperature and pressure(Water and ice phase
!          considered according to temperature). The definition
!          E/ESAT is used.
!
!Arguments
!
!          - Output -
! hr       relative humidity
!
!          - Input -
! qq       specific humidity
! tt       temperature in K
! ps       pressure in Pa
! ni       horizontal dimension
! nk       vertical dimension
! n        number of points to process

      Integer k, i
      real*8  temp1

#define __FORTRAN__
#include "tdpack_func.h"

!***********************************************************************
!
      if (satuco) then
         Do k=1,nk
            Do i=1,n
               temp1   = FOEWF(tt(i,k))
               temp1   = exp(temp1)
               temp1   = FOMULTS(temp1, tt(i,k))
               hr(i,k) = FOHRX(qq(i,k), ps(i,k), temp1)
            Enddo
         Enddo
      else
         Do k=1,nk
            Do i=1,n
               temp1   = FOEWAF(tt(i,k))
               temp1   = aerk1w*exp(temp1)
               hr(i,k) = FOHRX(qq(i,k), ps(i,k), temp1)
            Enddo
         Enddo
      endif

!
!***********************************************************************
!
      End Subroutine mfohr4
