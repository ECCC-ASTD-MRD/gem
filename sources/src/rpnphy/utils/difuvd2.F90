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
!**S/P DIFUVD2
!
      SUBROUTINE DIFUVD2(P, A, B, C, D, DELTA, NP, N, NK)
!
      implicit none
#include <arch_specific.hf>
      INTEGER NP, N, NK
      REAL P(NP,NK), A(N,NK), B(N,NK), C(N,NK), D(N,NK), DELTA(N,NK)
!
!Author
!          R. Benoit(Mar 89)
!
!Object
!          this is the transverse vectorization of the subroutine
!          ROSSR1
!
!Arguments
!
!          - Output -
! P        result
!
!          - Input -
! A        lower-diagonal
! B        diagonal
! C        upper-diagonal
! D        right side
! DELTA    work space
! NP       1st dimension of P
! N        1st dimension of all other vectors
! NK       2nd dimension (order)
!
!Notes
!          This subroutine solves the tri-diagonal matrix problem
!          below using ROSSR:
!
!  *****                                      ****  ****     ***    ****
!  *****                                      ****  ****     ***    ****
!  ** B(1),C(1), 0  , 0  , 0  , - - - -  ,  0   **  **  P(1)  **    **D(
!  ** A(2),B(2),C(2), 0  , 0  , - - - -  ,  0   **  **  P(2)  **    **D(
!  **  0  ,A(3),B(3),C(3), 0  , - - - -  ,  0   **  **  P(3)  **    **
!  **  0  , 0  ,A(4),B(4),C(4), - - - -  ,  0   **  **  P(4)  **    **
!  **  0  , 0  , 0  ,A(5),B(5), - - - -  ,  0   **  **  P(5)  ** -- **
!  **  -                                    -   **  **    -   ** -- **
!  **  -                                    -   **  **    -   **    **
!  **  0  ,  , 0 ,A(NK-2),B(NK-2),C(NK-2),  0   **  ** P(NK-2)**    **
!  **  0  ,  , 0 ,  0   ,A(NK-1),B(NK-1),C(NK-1)**  ** P(NK-1)**    **
!  **  0  ,  , 0 ,  0   ,   0   , A(NK) , B(NK) **  **  P(NK) **    **D(
!  ****                                       ****  ****    ****    ****
!  ****                                       ****  ****    ****    ****
!
!          Delta is a working array of dimension NK.  If the array D
!          is not required subsequently then the call CALL ROSSR1(P.A,
!          B,C,D,D,NK) will reduce storage requirements.  If further
!          the array C is not required subsequently, then the call
!          will further reduce storage requirements.
!
!*
!
!-----------------------------------------------------------------------
!
      REAL X
      INTEGER I, K
!
      DO 10 I=1,N
      C(I,NK)=0.
      X=1./B(I,1)
      P(I,1)=-C(I,1)*X
      DELTA(I,1)=D(I,1)*X
   10 CONTINUE
!
      DO 20 K=2,NK
      DO 20 I=1,N
      X=1./(B(I,K)+A(I,K)*P(I,K-1))
      P(I,K)=-C(I,K)*X
      DELTA(I,K)=(D(I,K)-A(I,K)*DELTA(I,K-1))*X
   20 CONTINUE
!
      DO 30 I=1,N
      P(I,NK)=DELTA(I,NK)
   30 CONTINUE
      DO 40 K=NK-1,1,-1
      DO 40 I=1,N
      P(I,K)=P(I,K)*P(I,K+1)+DELTA(I,K)
   40 CONTINUE
!
      RETURN
      END
