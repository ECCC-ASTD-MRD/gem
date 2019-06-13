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
!**S/P DIFUVD4N
!
      SUBROUTINE DIFUVD4N (D, F, S, N, NK)
!
      implicit none
!!!#include <arch_specific.hf>
!
      INTEGER N,NK
      REAL D(N,NK), F(N,NK), S(N,NK)
!
!Author
!          S. Belair (February 1996)
!
!Object
!          Calculate a simple vertical derivative of a function
!          F.  The function F is defined on the sigma levels and
!          the resulting derivative is obtained on the
!          staggered levels (sigma').
!
!Revisions
! 001
!
!Arguments
!
!          - Output -
! D        vertical derivative
!
!          - Input -
! F        function to derive vertically
! S        sigma levels for U
! N        number of columns to process
! NK       vertical dimension
!
!*
!
      INTEGER I, K
      REAL HD
!
!
!
!
!                                  The derivatives
!
!
!                                  For K=1 to NK-1
!
!
!
      DO K=1,NK-1
        DO I=1,N
          HD = S(I,K+1)-S(I,K)
          D(I,K) = ( F(I,K+1)-F(I,K) ) / HD
        END DO
      END DO
!
!                                  For K=NK
!                                  NOTE:  This lowest level derivative
!                                         is not important for the case
!                                         of the atmospheric turbulent
!                                         fluxes since the boundary
!                                         condition is
!                                         flux = alpha + beta*var
!
      DO I=1,N
           D(I,NK) = D(I,NK-1)
      END DO
!
      RETURN
      END
