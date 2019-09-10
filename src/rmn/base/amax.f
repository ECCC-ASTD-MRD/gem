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
***FUNCTION AMAX  - FIND MAXIMUM OF ARRAY A
*
      FUNCTION AMAX(A,NI,NJ,N)
      REAL  A(NI,NJ)

*
*AUTHOR   - HERSH MITCHELL
*
*REVISION 001:  C. THIBEAULT  -  JUL 79  DOCUMENTATION
*                                        CALL TO VECLIB
*REVISION 002:  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE  - fortran
*
*OBJECT(AMAX)
*         - FINDS MAXIMUM VALUE OF A(NI,NJ) OMITTING N BORDER ROWS
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP      DECK=AMAX
*           OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - FUNCTION AMAX(A,NI,NJ,N)
*
*ARGUMENTS
*   IN    - A  - ARRAY WHERE MAXIMUM VALUE IS TO BE FOUND
*         - NI - X-DIMENSION
*         - NJ - Y-DIMENSION
*         - N  - NUMBER OF BORDER ROWS TO BE OMITTED
*
*-----------------------------------------------------------------------
*
*
*-----------------------------------------------------------------------
*
*
*  WHEN OMITTING N BORDER ROWS, SET INDICES
*
      IL = 1+N
      JL = 1+N
      IH = NI-N
      JH = NJ-N
      IF (N.EQ.0) JH = 1
      IF (N.EQ.0) IH = NI*NJ
      AMAX = A(IL,JL)
*
      DO 30 J=JL,JH
      DO 30 I=IL,IH
      AMAX = AMAX1(AMAX,A(I,J))
   30 CONTINUE
*
*-----------------------------------------------------------------------
*
      RETURN
      END
