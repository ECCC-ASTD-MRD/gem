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
***FUNCTION AMEAN - CALCULATES THE MEAN OF ARRAY A
*
      REAL FUNCTION AMEAN (A,NI,NJ,N)
      implicit none
      integer :: NI,NJ,N
      REAL A(NI,NJ)

*
*AUTHOR   - HERSH MITCHELL
*
*REVISION 001:  C. THIBEAULT - JUL 79  DOCUMENTATION AND CALL TO VECLIB
*REVISION 002:  C. THIBEAULT - MAR 83  CONVERSION AU CODE CRAY
*REVISION 003:  M. Valin     - JUIN 2015  implicit none + utilisation de Real*8 pour les calculs
*
*LANGUAGE  - fortran
*
*OBJECT(AMEAN)
*         - CALCULATES MEAN OF A(NI,NJ) OMITTING N BORDER ROWS
*
*LIBRARIES
*         - SOURCE   RMNSOURCELIB,ID=RMNP     DECK=AMEAN
*         - OBJECT   RMNLIB,ID=RMNP
*
*USAGE    - FUNCTION AMEAN(A,NI,NJ,N)
*
*ARGUMENTS
*  IN     - A  - ARRAY WHERE THE MEAN IS CALCULATED FROM
*         - NI - X-DIMENSION
*         - NJ - Y-DIMENSION
*         - N  - NUMBER OF BORDER ROWS TO BE OMITTED
*
*
*
*-------------------------------------------------------------------------------
*
      REAL *8 :: SUM
      integer :: I, J, IH, IL, JH, JL, PTS
      SUM = 0.
*
*
*  WHEN OMITTING N BORDER ROWS, SET INDICES
*
      IL=1+N
      JL=1+N
      IH=NI-N
      JH=NJ-N
      IF (N.EQ.0) JH=1
      IF (N.EQ.0) IH=NI*NJ
      PTS=(IH-IL+1)*(JH-JL+1)
*
      DO 10 J=JL,JH
      DO 10 I=IL,IH
      SUM = SUM + A(I,J)
   10 CONTINUE
*
      AMEAN=SUM/PTS
*
*-------------------------------------------------------------------------------
*
      RETURN
      END
