!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

subroutine DIFUVD1(E, SC, A, B, C, U, D, N, NU, NK)
   implicit none
!!!#include <arch_specific.hf>
   integer N, NU, NK
   real E(N, NK), A(N, NK), B(N, NK), C(N, NK), U(NU, NK), D(N, NK)
   real SC

   !@Author R. Benoit(Mar 89)

   !@Object produce the tri-diagonal: E=SC*(A,B,C)*U + D

   !@Arguments
   !          - Output -
   ! E        result of (product +D)
   !          - Input -
   ! SC       factor of the operator
   ! A        lower-diagonal
   ! B        diagonal
   ! C        upper-diagonal
   ! U        vector to multiply on matrix(A,B,C)
   ! D        inhomogeneous term to add
   ! N        1st dimension (for all except for U)
   ! NU       1st dimension of U
   ! NK       second dimension

   integer K, I

   do K=1,NK
      !     FAIRE D'ABORD TERME DIAGONAL ET AJOUT DU TERME D DANS ACCUMULATEUR
      do I=1,N
         E(I,K)=D(I,K)+SC*B(I,K)*U(I,K)
      enddo
      if (K.gt.1) then
         do I=1,N
            E(I,K)=E(I,K)+SC*A(I,K)*U(I,K-1)
         enddo
      endif
      if (K.lt.NK) then
         do I=1,N
            E(I,K)=E(I,K)+SC*C(I,K)*U(I,K+1)
         enddo
      endif
   enddo

   return
end subroutine DIFUVD1
