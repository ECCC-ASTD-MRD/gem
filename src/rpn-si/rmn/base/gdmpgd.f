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
***S/R GDMPGD  - MULTIPLIES AN ARRAY BY ANOTHER AND SCALES THE
*                RESULTANT ARRAY BY A CONSTANT.
*                   R(I,J) = CON * A(I,J) * B(I,J)
*
      SUBROUTINE GDMPGD (R,A,B,CON,NI,NJ,NB)
      REAL R(NI,NJ), A(NI,NJ), B(NI,NJ)

*
* ENTRY GDDVGD - DIVIDES AN ARRAY BY ANOTHER AND SCALES THE
*                RESULTANT ARRAY BY A CONSTANT.
*                   R(I,J) = CON * A(I,J) / B(I,J)
*
*AUTHOR   - H. MITCHELL  -  D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - OCT 79  DOCUMENTATION AND CALLS TO VECLIB
*REVISION 002   C. THIBEAULT - MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE - fortran
*
*OBJECT(GDMDGD)
*         - MULTIPLIES, ELEMENT-WISE, ONE 2-DIMENSIONAL ARRAY BY ANOTHER
*           AND SCALES THE RESULTANT ARRAY BY A CONSTANT OMITTING NB
*           BORDER ROWS OR IF ENTRY (GDDVGD) IS CALLED IT DIVIDES AN
*           ARRAY BY ANOTHER AND SCALES THE RESULT BY A CONSTANT, ALSO
*           OMITTING NB BORDER ROWS.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=GDMPGD
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL GDMPGD(R,A,B,CON,NI,NJ,NB)
*         - CALL GDDVGD(R,A,B,CON,NI,NJ,NB)
*
*ARGUMENTS
*   OUT   - R   - ARRAY CONTAINING RESULT
*   IN    - A   - FIRST ARRAY TO BE OPERATED ON
*         - B   - SECOND ARRAY TO BE OPERATED ON
*         - CON - SCALING FACTOR
*         - NI  - X-DIMENSION
*         - NJ  - Y-DIMENSION
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
*         - MULTIPLICATION OR DIVISION BY 1.0 IS NOT PERFORMED.
*         - NOTE THAT IF NB .EQ. 0, WE SET JB=JT=IL=1 AND IR=NI*NJ
*           AND WE HAVE THEREFORE A ONE TRIP LOOP WHICH IS MORE
*           EFFICIENT.
*
*-------------------------------------------------------------------------------
*
*
*-------------------------------------------------------------------------------
*
*  SETS R = CON*A*B
*
      ASSIGN 20 TO KT
      GO TO 10
*
      ENTRY GDDVGD(R,A,B,CON,NI,NJ,NB)
*
*  SETS R = CON*A/B
*
      ASSIGN 60 TO KT
*
*  FIRST PICK OUT DIGITS IN NB AND SET LOOP INDEXES
*
   10 IL=1+NB/1000
      JB=1+MOD(NB,1000)/100
      IR=NI-MOD(NB,100)/10
      JT=NJ-MOD(NB,10)
      IF (NB.EQ.0) JT = 1
      IF (NB.EQ.0) IR = NI*NJ
*
      GO TO KT,(20,60)
*
*  GRID MULTIPLICATION LOOPS
*
   20 IF (CON.EQ.1.0) GO TO 40
      DO 30 J=JB,JT
      DO 30 I=IL,IR
      R(I,J) = CON * A(I,J) * B(I,J)
   30 CONTINUE
      RETURN
*
   40 DO 50 J=JB,JT
      DO 50 I=IL,IR
      R(I,J) = A(I,J) * B(I,J)
   50 CONTINUE
      RETURN
*
*  GRID DIVISION LOOPS
*
   60 IF (CON.EQ.1.0) GO TO 80
      DO 70 J=JB,JT
      DO 70 I=IL,IR
      R(I,J) = A(I,J) / B(I,J)
      R(I,J) = CON * R(I,J)
   70 CONTINUE
      RETURN
*
   80 DO 90 J=JB,JT
      DO 90 I=IL,IR
      R(I,J) = A(I,J) / B(I,J)
   90 CONTINUE
      RETURN
*
*-------------------------------------------------------------------------------
*
      END
