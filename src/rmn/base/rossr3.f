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
***S/R ROSSR3 - SOLVE THE LINEAR SYSTEM OF EQUATION.
*
      SUBROUTINE ROSSR3(P,A,DELTA,C,D,M)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001  C. THIBEAULT  -  MAR 80  DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(ROSSR3)
*         - SOLVE A TRI-DIAGONAL SYSTEM OF LINEAR EQUATIONS.
*
*LIBRARIES
*         - SOURCE LIBRARY  RMNSOURCELIB,ID=RMNP     DECK=ROSSR3
*         - OBJECT LIBRARY  RMNLIB,ID=RMNP
*
*USAGE    - CALL ROSSR3(P,A,DELTA,C,D,M)
*
*ARGUMENTS
*         - P     - SOLUTION VECTOR.
*         - A     - VECTOR OF SUB-DIAGONAL ELEMENTS.
*         - DELTA - WORK VECTOR.
*         - C     - VECTOR OF SUPER-DIAGONAL ELEMENTS.
*         - D     - RIGHT-HAND SIDE VECTOR.
*         - M     - NUMBER OF EQUATIONS TO BE SOLVED.
*
*NOTES    -   I) THE GENERAL TRI-DIAGONAL PROBLEM IS REDUCIBLE TO
*                THE ABOVE FORM BY SUITABLY SCALING I.E.
*                A(I)P(I-1)+B(I)P(I)+C(I)P(I+1) = D(I)
*                SHOULD BE MULTIPLIED BY 1/B(I) (ASSUMED NON-ZERO.)
*            II) ALGORITHM USED IS AN EFFICIENT SPECIALISED
*                APPLICATION OF GAUSSIAN ELIMINATION AND BACK-
*                SUBSTITUTION WITHOUT PIVOTING.
*           III) IT IS COMPUTATIONALLY STABLE WHEN THE MATRIX IS
*                DIAGONALLY DOMINANT.
*                I.E.  ABS(A(I) + C(I)) .LE. 1  WILL ENSURE
*                COMPUTATIONAL STABILITY.
*            IV) OPERATION COUNT - 5M MULTIPLICATIONS,
*                3M ADDITIONS, M DIVISIONS.
*             V) A(1) AND C(M) ARE NOT USED.
*            VI) IF ANY OF A, D OR C NEED NOT BE PRESERVED, THEN
*                THEY MAY BE USED FOR THE ACTUAL PARAMETER
*                SUPPLIED FOR WORK AND THUS REDUCE STORAGE
*                REQUIREMENTS.
*           VII) IF A SYSTEM OF LINEAR EQUATIONS HAS TO BE SOLVED
*                REPEATEDLY WITH DIFFERENT RIGHT-HAND-SIDES, THE
*                PAIR SETTRI - SOLTRI WILL OFFER A BETTER ALTERNATIVE.
*         - ALL PARAMETERS ARE IN SCM (LEVEL 1).
*
*-------------------------------------------------------------------------------
*
      REAL      P(M)
      REAL      A(M)
      REAL  DELTA(M)
      REAL      C(M)
      REAL      D(M)
*
*-------------------------------------------------------------------------------
*
*  DELTA IS A WORKING STORAGE ARRAY OF DIMENSION M
*  IF THE VECTOR D IS NOT REQUIRED SUBSEQUENTLY THEN THE CALL STATEMENT
*      CALL ROSSR3(P,A,D,C,D,M)
*  WILL USE THE ARRAY D AS WORKING STORAGE AND REDUCE THE OVERALL
*  STORAGE REQUIRED
*
*
*  IF THE ARRAY C IS NOT REQUIRED SUBSEQUENTLY THEN THE CALL STATEMENT
*      CALL ROSSR3(C,A,DELTA,C,D,M)
*  WILL REDUCE CORE STORAGE REQUIREMENTS
*  IF BOTH C AND D ARE NOT REQUIRED SUBSEQUENTLY THEN THE CALL STATEMENT
*      CALL ROSSR3(C,A,D,C,D,M)
*  WILL FURTHER REDUCE THE CORE REQUIREMENTS
*
*
      M1=M-1
      C(M)=0.
      P(1)=-C(1)
      DELTA(1)=D(1)
*
*
      DO 1 I=1,M1
      I1=I+1
      X = 1. / (1.+A(I1)*P(I))
      P(I1)=-C(I1)*X
    1 DELTA(I1)=(D(I1)-A(I1)*DELTA(I))*X
*
*
      P(M)=DELTA(M)
*
      DO 2 I=1,M1
      II=M-I
    2 P(II)=P(II)*P(II+1)+DELTA(II)
*
*-----------------------------------------------------------------------
*
      RETURN
      END
