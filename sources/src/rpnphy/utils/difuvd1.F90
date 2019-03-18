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
!**S/P DIFUVD1
!
      SUBROUTINE DIFUVD1 (E, SC, A, B, C, U, D, N, NU, NK)
!
      implicit none
#include <arch_specific.hf>
      INTEGER N, NU, NK
      REAL E(N, NK), A(N, NK), B(N, NK), C(N, NK), U(NU, NK), D(N, NK)
      REAL SC
!
!Author
!          R. Benoit(Mar 89)
!
!Object
!          to produce the tri-diagonal: E=SC*(A,B,C)*U + D
!
!Arguments
!
!          - Output -
! E        result of (product +D)
!
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
!
!
!LANGAGE CFT
!*
!
      INTEGER K, I
!
      DO 10 K=1,NK
!     FAIRE D'ABORD TERME DIAGONAL ET AJOUT DU TERME D DANS ACCUMULATEUR
         DO 1 I=1,N
1           E(I,K)=D(I,K)+SC*B(I,K)*U(I,K)
         IF (K.GT.1) THEN
            DO 2 I=1,N
2              E(I,K)=E(I,K)+SC*A(I,K)*U(I,K-1)
         ENDIF
         IF (K.LT.NK) THEN
            DO 3 I=1,N
3              E(I,K)=E(I,K)+SC*C(I,K)*U(I,K+1)
         ENDIF
10       CONTINUE
!
      RETURN
      END
