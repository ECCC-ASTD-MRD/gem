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
***S/R GDADCN - ADDS A CONSTANT TO AN ARRAY
*                   R(I,J) = CON + A(I,J)
*
      SUBROUTINE GDADCN (R,A,CON,NI,NJ,NB)
      REAL R(NI,NJ), A(NI,NJ)

*
*
* ENTRY GDMPCN - MULTIPLY AN ARRAY BY A CONSTANT
*                   R(I,J) = CON * A(I,J)
*
*AUTHOR   - H. MITCHELL  -  D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - OCT 79   DOCUMENTATION
*REVISION 002   C. THIBEAULT - MAR 83   CONVERSION AU CODE CRAY
*
*LANGUAGE - fortran
*
*OBJECT(GDADCN)
*         - ADDS A CONSTANT TO EVERY ELEMENT OF A 2-DIMENSIONAL ARRAY
*           OMITTING NB BORDER ROWS OR ENTRY (GDMPCN) MULTIPLY AN ARRAY
*           BY A CONSTANT.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=GDADCN
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL GDADCN(R,A,CON,NI,NJ,NB)
*         - CALL GDMPCN(R,A,CON,NI,NJ,NB)
*
*ARGUMENTS
*   OUT   - R   - ARRAY CONTAINING RESULT
*   IN    - A   - ARRAY TO BE OPERATED ON
*         - CON - CONSTANT
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
*NOTES    - R AND A MAY CO-INCIDE.  ADDITIONS OF 0.0 AND MULTIPLICATIONS
*           BY 1.0 ARE NOT PERFORMED.
*         - NOTE THAT IF NB .EQ. 0, WE SET JB=JT=IL=1 AND IR=NI*NJ
*           AND WE HAVE THEREFORE A ONE TRIP LOOP WHICH IS MORE
*           EFFICIENT.
*
*-------------------------------------------------------------------------------
*
*
*-------------------------------------------------------------------------------
*
      ASSIGN 20 TO KT
      GO TO 10
*
      ENTRY GDMPCN(R,A,CON,NI,NJ,NB)
*
      ASSIGN 40 TO KT
*
   10 IL=1+NB/1000
      JB=1+MOD(NB,1000)/100
      IR=NI-MOD(NB,100)/10
      JT=NJ-MOD(NB,10)
      IF (NB.EQ.0) JT = 1
      IF (NB.EQ.0) IR = NI*NJ
*
      GO TO KT,(20,40)
*
*  ADDITION LOOP
*
   20 IF (CON.EQ.0.0) GO TO 60
      DO 30 J=JB,JT
      DO 30 I=IL,IR
      R(I,J) = CON + A(I,J)
   30 CONTINUE
      RETURN
*
*  MULTIPLICATION LOOP
*
   40 IF (CON.EQ.1.0) GO TO 60
      DO 50 J=JB,JT
      DO 50 I=IL,IR
      R(I,J) = CON * A(I,J)
   50 CONTINUE
      RETURN
*
*  TRANSFER LOOP
*
   60 DO 70 J=JB,JT
      DO 70 I=IL,IR
      R(I,J) = A(I,J)
   70 CONTINUE
      RETURN
*
*-------------------------------------------------------------------------------
*
      END
