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
***S/R SQADSQ - COMPUTES THE SQUARES OF TWO ARRAYS, SCALES EACH RESULT
*               AND SUMS.
*                  R(I,J) = CONA*(A(I,J)**2) + CONB*(B(I,J)**2)
*
      SUBROUTINE SQADSQ (R,A,B,CONA,CONB,NI,NJ,NB)
      REAL R(NI,NJ), A(NI,NJ), B(NI,NJ)

*
*AUTHOR   - H. MITCHELL  -  D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - OCT 79  DOCUMENTATION AND CALLS TO VECLIB
*REVISION 002   C. THIBEAULT - MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE - fortran
*
*OBJECT(SQADSQ)
*         - COMPUTES THE SQUARES, ELEMENT WISE, OF TWO 2-DIMENSIONAL
*           ARRAYS, SCALES EACH RESULT AND SUMS THE SCALED ARRAYS.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=SQADSQ
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL SQADSQ(R,A,B,CONA,CONB,NI,NJ,NB)
*
*ARGUMENTS
*   OUT   - R    - ARRAY CONTAINING RESULT
*   IN    - A    - FIRST ARRAY TO BE OPERATED ON
*         - B    - SECOND ARRAY TO BE OPERATED ON
*         - CONA - SCALING FACTOR FOR ARRAY A
*         - CONB - SCALING FACTOR FOR ARRAY B
*         - NI   - X-DIMENSION
*         - NJ   - Y-DIMENSION
*         - NB  - 4-DIGIT DECIMAL NUMBER WHICH ALLOWS AT MOST NINE ROWS
*                 ON ANY SIDE OF THE GRID TO BE SKIPPED.
*               - NB = D1*10E3 + D2*10E2 + D3*10 + D4
*               - D1 = COLUMNS ARE SKIPPED FROM LEFT
*               - D2 = ROWS ARE SKIPPED FROM BOTTOM
*               - D3 = COLUMNS ARE SKIPPED FROM RIGHT
*               - D4 = ROWS ARE SKIPPED FROM TOP
*               - E.G.  NB=0000 INCLUDE ALL BORDERS (WRITTEN AS 0)
*                       NB=1111 OMIT ONE ROW ALL AROUND THE GRID
*                       NB=2241 INDICATES THAT THE 2 COLUMNS ON THE
*                       LEFT OF A, 2 ROWS AT THE BOTTOM, 4 COLUMNS
*                       AT THE RIGHT AND 1 ROW AT THE TOP ARE NOT
*                       OPERATED ON AND THE CORRESPONDING ROWS
*                       AND COLUMNS OF R ARE LEFT UNDEFINED BY
*                       THE ROUTINE.
*
*NOTES    - ANY OF R,A,B MAY CO-INCIDE.
*         - IF CONA = CONB = 1.0, SCALING IS NOT PERFORMED.
*         - NOTE THAT IF NB .EQ. 0, WE SET JB=JT=IL=1 AND IR=NI*NJ
*           AND WE HAVE THEREFORE A ONE TRIP LOOP WHICH IS MORE
*           EFFICIENT.
*         - WARNING *** VECLIB IS BANK SENSITIVE.
*
*-------------------------------------------------------------------------------
*
*
*-------------------------------------------------------------------------------
*
      IL=1+NB/1000
      JB=1+MOD(NB,1000)/100
      IR=NI-MOD(NB,100)/10
      JT=NJ-MOD(NB,10)
      IF (NB.EQ.0) JT = 1
      IF (NB.EQ.0) IR = NI*NJ
*
      IF (CONA.EQ.1.AND.CONB.EQ.1.) GO TO 20
*
      DO 10 J=JB,JT
      DO 10 I=IL,IR
   10 R(I,J)=CONA*A(I,J)**2+CONB*B(I,J)**2
      RETURN
*
   20 DO 30 J=JB,JT
      DO 30 I=IL,IR
      R(I,J) = A(I,J) * A(I,J) + B(I,J) * B(I,J)
   30 CONTINUE
      RETURN
*
*-------------------------------------------------------------------------------
*
      END
