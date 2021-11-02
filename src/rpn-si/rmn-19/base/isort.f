*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
***S/R ISORT - RE=ARRANGES THE ELEMENTS OF A REAL ARRAY INTO
*             NON-DECREASING ORDER.
*
      SUBROUTINE ISORT (A,N)
      INTEGER A      (N)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*REVISION 002   M. LEPINE  -  VERSION SORT POUR UN TABLEAU D'ENTIERS
*Revision 003   M. Lepine - Avril 2007 - traitement du cas (ridicule) avec n <= 1
*
*LANGUAGE - fortran
*
*OBJECT(SORT)
*         - RE-ARRANGES THE ELEMENTS OF A INTEGER ARRAY INTO
*           NON-DECREASING ORDER.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=SORT
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL SORT(A,N)
*
*ARGUMENTS
* IN/OUT  - A - INTEGER ARRAY TO BE SORTED
*   IN    - N - NUMBER OF ELEMENTS IN A
*
*NOTES    - SORT RE-ORDERS A SUCH THAT A(1).LE.A(2).LE. ... .LE.A(N).
*         - AN EFFICIENT ALGORITHM, CALLED TREESORT3, IS USED BY SORT.
*           THE ALGORITHM PERFORMS THE SORTING IN A NUMBER OF OPERATIONS
*           PROPORTIONAL TO N LOG N WITHOUT USING ANY EXTRA STORAGE.
*           (CACM COLLECTED ALGORITHMS #245.). IT IS ALSO REASONABLY
*           COMPACT.
*         - ANOTHER SORTING ROUTINE (CALLED IPSORT) IS AVAILABLE.
*           IPSORT RETURNS IP SUCH THAT A(IP(1)).LE.A(IP(2)).LE. ...
*           .LE.A(IP(N)).
*
**


*-----------------------------------------------------------------------------
      INTEGER ASAVE,N,II,NN,IRETRN,I,J
*-----------------------------------------------------------------------------

      IF (n .le. 1) RETURN
*
* INITIALIZATION FOR FIRST LOOP (IF  N/2.GE.2 ONLY).
*
      IF (N.LE.3) GO TO 30
      II = N/2+1
      NN = N
      ASSIGN 20 TO IRETRN
*
   20 II = II-1
      IF (II.GE.2) GO TO 60
*
* INITIALIZATION OF SECOND LOOP.
*
   30 II = 1
      NN = N+1
      ASSIGN 50 TO IRETRN
*
   40 NN = NN-1
      GO TO 60
*
* EXCHANGE A(1),A(NN)
*
   50 ASAVE = A(1)
      A(1)  = A(NN)
      A(NN) = ASAVE
*
      IF (NN.GE.3) GO TO 40
*
      RETURN
*
* THE FOLLOWING IS THE CODE FOR SIFTUP(II,NN)
*
*
   60 ASAVE = A(II)
      I = II
   70 J = I+I
      IF (J-NN) 80,90,100
   80 IF (A(J+1).GT.A(J)) J=J+1
   90 IF (A(J).LE.ASAVE)  GO TO 100
      A(I) = A(J)
      I    =   J
      GO TO 70
  100 A(I) = ASAVE
      GO TO IRETRN,(20,50)
*
*----------------------------------------------------------------------------
*
      END
