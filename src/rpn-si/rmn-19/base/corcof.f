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
***S/R CORCOF - FINDS THE CORRELATION COEFFICIENT
*
      SUBROUTINE CORCOF(CC,FI,FF,FV,W,NI,NJ,IW1,IW2,NW1,NW2)
      implicit none
      integer :: NI,NJ,IW1,IW2,NW1,NW2
      REAL FI(NI,NJ),FF(NI,NJ),FV(NI,NJ),W(NI,NJ)
      REAL :: CC

*
*AUTHOR   - Y. BOURASSA  -  APR 75
*
*REVISION 001  C. THIBEAULT  -  SEP 79  DOCUMENTATION
*REVISION 002  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*REVISION 003  M. Valin      -  JUIN 2015  implicit none + utilisation de Real*8 pour les calculs
*
*OBJECT(CORCOF)
*         _ GIVEN AN INITIAL FIELD, A FORECAST AND A VERIFYING FIELD,
*           FINDS THE CORRELATION COEFFICIENT.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP      DECK=CORCOF
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL CORCOF (CC,FI,FF,FV,W,NI,NJ,IW1,IW2,NW1,NW2)
*
*ARGUMENTS
*   OUT   - CC  - CORRELATION COEFFICIENT
*   IN    - FI  - INITIAL FIELD
*         - FF  - FORECAST FIELD
*         - FV  - VERIFYING FIELD
*         - W   - WEIGHT FIELD
*         - NI  - X-DIMENSION
*         - NJ  - Y-DIMENSION
*         - IW1 - BOTTOM LEFT CORNER OF X-COORDINATE OF WINDOW
*         - IW2 - BOTTOM LEFT CORNER OF Y-COORDINATE OF WINDOW
*         - NW1 - UPPER RIGHT CORNER OF X-COORDINATE OF WINDOW
*         - NW2 - UPPER RIGHT CORNER OF Y-COORDINATE OF WINDOW
*
*NOTES    - IF THE COORDINATES OF THE WINDOW ARE OUTSIDE THE LIMITS OF
*           THE ARRAYS, CC IS SET TO 99999.
*         - CC IS THE CORRELATION COEFFICIENT BETWEEN THE ACTUAL
*           CHANGE IN THE FIELD (DA = FV-FI) AND THE PREDICTED
*           CHANGE IN THE FIELD (DF = FF-FI).  TO OBTAIN CC,
*           REMOVE THE (W-WEIGHTED) MEANS FROM DA AND DF AND
*           THEN COMPUTE
*            CC = MEAN(DA*DF) / SQRT(MEAN(DA)**2 * MEAN(DF)**2)
*           IF DA, DF OR W ARE PATHOLOGICAL, CC IS SET TO 99999.
*
*-----------------------------------------------------------------------------
*
*
*-------------------------------------------------------------------------------
*
      REAL *8 R(5), FFI, DA, DF, DW, FTW, X, Y, A, B, C
      integer :: I, J
      CC = 99999.
      DO 10 I=1,5
      R(I) = 0.0
   10 CONTINUE
*
*
*  TEST TO SEE IF THE WINDOW CO-ORDINATES ARE INSIDE THE
*  LIMITS OF THE ARRAYS DIMENSIONED (NI*NJ)
*
      IF (NI.LT.1.OR.NI.LT.IW1.OR.NI.LT.NW1) RETURN
      IF (NJ.LT.1.OR.NJ.LT.IW2.OR.NJ.LT.NW2) RETURN
      IF (NW1.LT.IW1.OR.NW2.LT.IW2) RETURN
*
      FTW = 0.
*
        DO 20 J=IW2,NW2
        DO 20 I=IW1,NW1
        FFI= FI(I,J)
        DA = FV(I,J)-FFI
        DF = FF(I,J)-FFI
        DW = W(I,J)
        FTW = FTW + DW
        X = DW*DA
        Y = DW*DF
        R(1) = R(1) + (X*DF)
        R(2) = R(2) + (Y*DF)
        R(3) = R(3) + (X*DA)
        R(4) = R(4) + (Y)
        R(5) = R(5) + (X)
   20   CONTINUE
*
      IF (FTW.EQ.0) RETURN
*
        DO 30 I=1,5
        R(I) = R(I)/FTW
   30   CONTINUE
*
      A = R(1) - R(4)*R(5)
      B = R(2) - R(4)*R(4)
      C = R(3) - R(5)*R(5)
      IF (B.EQ.0..OR.C.EQ.0.) RETURN
      CC = A/SQRT(B*C)
*
*-------------------------------------------------------------------------------
*
      RETURN
      END
