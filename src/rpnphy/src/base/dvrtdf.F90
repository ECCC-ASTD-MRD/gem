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

subroutine DVRTDF(R , X , DS, N , MR , MX , NK)
   implicit none
!!!#include <arch_specific.hf>
   integer N, MR, MX, NK
   real R(MR,NK),X(MX,NK),DS(n,NK)

   !@Author R. Benoit RPN(Mar 1989)

   !@Revisions
   ! 001      R. Benoit (Aug 93) - DS(2D) for Local sigma

   !@Object calculate the vertical derivative by centred finite differences

   !@Arguments
   !          - Output -
   ! R        result
   !          - Input -
   ! X        variable to derive
   ! DS       distance between sigma levels 'U'
   ! N        horizontal dimensions
   ! MR       1st dimension of R
   ! MX       1st dimension of X
   ! NK       vertical dimension

   !@Notes
   !          R and X can share the same space, R(*,NK)=0

   integer J,K

   do K=1,NK-1
      do J=1,N
         R(J,K)=(X(J,K+1)-X(J,K))/DS(j,K)
      enddo
   enddo

   do J=1,N
      R(J,NK)=0
   enddo

   return
end subroutine DVRTDF
