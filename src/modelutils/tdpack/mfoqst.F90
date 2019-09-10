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
subroutine mfoqst3(qs,tt,ps,ni,nk,n)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack
   implicit none
   !@object calcule humidite specifique saturante.
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"
   !@arguments
   !          - Output -
   ! qs       saturation specific humidity in kg/kg
   !
   !          - Input -
   ! tt       temperature in K
   ! ps       pression in Pa
   ! ni       horizontal dimension
   ! nk       vertical dimension
   ! n        number of points to process
   integer, intent(in) :: ni, nk, n
   real, intent(out) :: qs(ni,nk)
   real, intent(in) :: tt(ni,nk)
   real, intent(in) :: ps(ni,*)
   !@author N. Brunet  (Jan91)
   !@revision
   ! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
   ! 002      B. Bilodeau (January 2001) - Automatic arrays
   ! 003      L. Spacek   (May 2003)     - IBM conversion
   !                - calls to vexp routine (from massvp4 library)
   !
   !@description
   !          to calculate saturation specific humidity (water and ice
   !          phase considered according to temperature)
   !
   !@/
   integer k, i
   real(REAL64) :: dtemp
   real(REAL64), dimension(ni,nk) :: xt

!!$#define __FORTRAN__
!!$#include "tdpack_func.h"

   !--------------------------------------------------------------------
   do k=1,nk
      do i=1,ni
         xt(i,k) = FOEWF(tt(i,k))
      enddo
   enddo

   call vexp(xt,xt,ni*nk)

   do k=1,nk
      do i=1,n
         dtemp = FOMULTS(xt(i,k),tt(i,k))
         qs(i,k) = FOQSTX(ps(i,k),dtemp)
      enddo
   enddo
   !--------------------------------------------------------------------
   return
end subroutine mfoqst3
