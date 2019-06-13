!-------------------------------------- LICENCE BEGIN -------------------------
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
!-------------------------------------- LICENCE END ---------------------------

!/@
subroutine mfohr4 (hr, qq, tt, ps, ni, nk, n, satuco)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@objective compute relative humidity from specific humidity
   !           temperature and pression
   !@arguments
   !          - Output -
   ! hr       relative humidity
   !          - Input -
   ! qq       specific humidity
   ! tt       temperature in K
   ! ps       pressure in Pa
   ! ni       horizontal dimension
   ! nk       vertical dimension
   ! n        number of points to process
   logical, intent(in) :: satuco
   integer, intent(in) ::  ni, nk, n
   real, intent(in) ::  qq(ni,*), tt(ni,*), ps(ni,*)
   real, intent(out) ::  hr(ni,*)
   !@author N. Brunet  (Jan91)
   !@revision
   ! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
   ! 002      B. Bilodeau  (January 2001) - Automatic arrays
   ! 003      G. Pellerin  (June 2003) - IBM conversion
   !                  - calls to vexp routine (from massvp4 library)
   !@description
   !          to calculate relative humidity from specific humidity,
   !          temperature and pressure(Water and ice phase
   !          considered according to temperature). The definition
   !          E/ESAT is used.
   !@/
   integer :: k, i
   real(REAL64) :: temp1

!!$#define __FORTRAN__
!!$#include "tdpack_func.h"

   !---------------------------------------------------------------------
   if (satuco) then
      do k=1,nk
         do i=1,n
            temp1   = FOEWF(tt(i,k))
            temp1   = exp(temp1)
            temp1   = FOMULTS(temp1, tt(i,k))
            hr(i,k) = FOHRX(qq(i,k), ps(i,k), temp1)
         enddo
      enddo
   else
      do k=1,nk
         do i=1,n
            temp1   = FOEWAF(tt(i,k))
            temp1   = aerk1w*exp(temp1)
            hr(i,k) = FOHRX(qq(i,k), ps(i,k), temp1)
         enddo
      enddo
   endif
   !---------------------------------------------------------------------
   return
end subroutine mfohr4
