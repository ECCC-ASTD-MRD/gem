!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!**S/R MTXINV - INVERT A MATRIX 
! 
      SUBROUTINE  EZ_MTXINV8 (A,L,M,N) 
      implicit none
      INTEGER N 
      REAL*8 A(N,N) 
      INTEGER L(N), M(N)
! 
!Author
!          R. Sarrazin (1992)REAL*8 version for MTXINV
!
!Object 
!          to invert a matrix (in double precision)
!
!Revision
! 001      C. Thibeault (Apr 80) Documentation
! 002      M. Valin (Nov 86) More logical algorithm
! 003      R. Sarrazin (1992) Real*8 of MTXINV
!
!Arguments
!
!          - Input/Output -
! A        input matrix, destroyed in computation and replaced by 
!          resultant inverse.  A must be a general matrix.
!
!          - Input -
! L        work vector of length N
! M        work vector of length N
! N        order of matrix A
!
!Notes
!          This routine based closely upon IBM-SSP subroutine MTXINV.  
!          The standard Gauss-Jordan method is used.  If the matrix 
!          is singular, a message is printed and program execution is 
!          terminated.
!
!*
      INTEGER I, J, K 
      REAL*8 HOLD, BIGA 
! 
!        SEARCH FOR LARGEST ELEMENT 
! 
      DO 80 K=1,N 
        L(K)=K
        M(K)=K
        BIGA=A(K,K) 
        DO 20 J=K,N 
        DO 20 I=K,N 
          IF( ABS(BIGA).LT.ABS(A(I,J))) THEN
            BIGA=A(I,J) 
            L(K)=I
            M(K)=J
          ENDIF 
20      CONTINUE
! 
!     INTERCHANGE ROWS
! 
      J=L(K)
      IF(J.GT.K) THEN 
        DO 30 I=1,N 
          HOLD=-A(K,I)
          A(K,I)=A(J,I) 
30        A(J,I) =HOLD
      ENDIF 
! 
!     INTERCHANGE COLUMNS 
! 
      I=M(K)
      IF(I.GT.K) THEN 
        DO 40 J=1,N 
          HOLD=-A(J,K)
          A(J,K)=A(J,I) 
40        A(J,I) =HOLD
      ENDIF 
! 
!     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS 
!     CONTAINED IN BIGA)
! 
      IF(BIGA.EQ.0)THEN 
        WRITE (6,601) 
        STOP
      ENDIF 
      DO 51 I=1,K-1 
51    A(I,K)=A(I,K)/(-BIGA) 
      DO 52 I=K+1,N 
52    A(I,K)=A(I,K)/(-BIGA) 
! 
!        REDUCE MATRIX
! 
      DO 65 I=1,N 
        HOLD=A(I,K) 
        IF(I.NE.K)THEN
          DO 61 J=1,K-1 
61        A(I,J)=HOLD*A(K,J)+A(I,J) 
          DO 62 J=K+1,N 
62        A(I,J)=HOLD*A(K,J)+A(I,J) 
        ENDIF 
65    CONTINUE
! 
!     DIVIDE ROW BY PIVOT 
! 
      DO 71 J=1,K-1 
71    A(K,J)=A(K,J)/BIGA
      DO 72 J=K+1,N 
72    A(K,J)=A(K,J)/BIGA
! 
!     REPLACE PIVOT BY RECIPROCAL 
! 
      A(K,K)=1.0/BIGA 
   80 CONTINUE
! 
!     FINAL ROW AND COLUMN INTERCHANGE
! 
      DO 150 K=N-1,1,-1 
        I=L(K)
        IF(I.GT.K) THEN 
          DO 110 J=1,N
            HOLD=A(J,K) 
            A(J,K)=-A(J,I)
110         A(J,I) =HOLD
        ENDIF 
        J=M(K)
        IF(J.GT.K) THEN 
          DO 130 I=1,N
            HOLD=A(K,I) 
            A(K,I)=-A(J,I)
130         A(J,I) =HOLD
        ENDIF 
150   CONTINUE
      RETURN
! 
601   FORMAT ('0**ERROR** SUBROUTINE MTXINV WAS CALLED TO INVERT A',          ' SINGULAR MATRIX')
! 
      END 
