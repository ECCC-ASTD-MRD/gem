
subroutine DIFUVD2(P, A, B, C, D, DELTA, NP, N, NK)

   implicit none
!!!#include <arch_specific.hf>
   integer NP, N, NK
   real P(NP,NK), A(N,NK), B(N,NK), C(N,NK), D(N,NK), DELTA(N,NK)

   !@Author R. Benoit(Mar 89)
   !@Object transverse vectorization of the subroutine ROSSR1
   !@Arguments
   !          - Output -
   ! P        result
   !          - Input -
   ! A        lower-diagonal
   ! B        diagonal
   ! C        upper-diagonal
   ! D        right side
   ! DELTA    work space
   ! NP       1st dimension of P
   ! N        1st dimension of all other vectors
   ! NK       2nd dimension (order)
   !@Notes
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
   !-----------------------------------------------------------------------

   real X
   integer I, K

   do I=1,N
      C(I,NK)=0.
      X=1./B(I,1)
      P(I,1)=-C(I,1)*X
      DELTA(I,1)=D(I,1)*X
   end do

   do K=2,NK
      do I=1,N
         X=1./(B(I,K)+A(I,K)*P(I,K-1))
         P(I,K)=-C(I,K)*X
         DELTA(I,K)=(D(I,K)-A(I,K)*DELTA(I,K-1))*X
      enddo
   enddo

   do I=1,N
      P(I,NK)=DELTA(I,NK)
   enddo

   do K=NK-1,1,-1
      do I=1,N
         P(I,K)=P(I,K)*P(I,K+1)+DELTA(I,K)
      enddo
   enddo

   return
end subroutine DIFUVD2
