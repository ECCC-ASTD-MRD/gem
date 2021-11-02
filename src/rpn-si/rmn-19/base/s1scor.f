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
***S/R S1SCOR - COMPUTES THE S1 (GRADIENT) SCORE
*
      SUBROUTINE S1SCOR(S1,F,V,NI,NJ,IW1,IW2,NW1,NW2,IGL)
      implicit none
      REAL F(NI,NJ), V(NI,NJ)
      integer :: NI,NJ,IW1,IW2,NW1,NW2,IGL
      REAL :: S1

*
*AUTHOR   - Y.R. BOURASSA  -  APR 75
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*REVISION 002   C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*REVISION 003   M. Valin      -  JUIN 2015  implicit none + utilisation de Real*8 pour les calculs
*
*LANGUAGE - fortran
*
*OBJECT(S1SCOR)
*         - GIVEN A FORECAST, A VERIFYING FIELD AND THE CO-ORDINATES
*           OF THE WINDOW, COMPUTES THE S1 (GRADIENT) SCORE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=S1SCOR
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL S1SCOR(S1,F,V,NI,NJ,IW1,IW2,NW1,NW2,IGL)
*
*ARGUMENTS
*   OUT   - S1  - S1 SCORE
*   IN    - F   - FORECAST FIELD
*         - V   - VERIFYING FIELD
*         - NI  - X-DIMENSION
*         - NJ  - Y-DIMENSION
*         - IW1 - BOTTOM LEFT CORNER OF X-COORDINATE OF WINDOW
*         - IW2 - BOTTOM LEFT CORNER OF Y-COORDINATE OF WINDOW
*         - NW1 - UPPER RIGHT CORNER OF X-COORDINATE OF WINDOW
*         - NW2 - UPPER RIGHT CORNER OF Y-COORDINATE OF WINDOW
*         - IGL - GRID POINT INTERVAL FOR VERIFICATION (SET TO 1 OR 2)
*
*NOTES    - IF THE CO-ORDINATES OF THE WINDOW ARE OUTSIDE THE LIMIT OF
*           THE ARRAYS OR IF IGL IS NOT EQUAL TO 1 OR 2, S1 IS SET
*           TO 99999.
*
*------------------------------------------------------------------------------
*
*
*------------------------------------------------------------------------
*
      REAL*8 A, B, VIJ, FIJ, FI1J, VI1J, FIJ1, VIJ1
      integer :: I, J, I2, J2
      S1 = 99999.
*
      IF (IGL.NE.1.AND.IGL.NE.2) RETURN
      IF (NW1-IW1.LT.IGL.OR.NW2-IW2.LT.IGL) RETURN
*
*  TEST TO SEE IF THE WINDOW CO-ORDINATES ARE INSIDE THE
*  LIMITS OF THE ARRAYS DIMENSIONED (NI*NJ)
*
      IF (NI.LT.1.OR.NI.LT.IW1.OR.NI.LT.NW1) RETURN
      IF (NJ.LT.1.OR.NJ.LT.IW2.OR.NJ.LT.NW2) RETURN
      IF (NW1.LT.IW1.OR.NW2.LT.IW2) RETURN
*
      A = 0.
      B = 0.
      I2 = NW1-IGL
      J2 = NW2-IGL
*
      DO 10 J = IW2,J2,IGL
      DO 10 I = IW1,I2,IGL
      VIJ  = V(I,J)
      FIJ  = F(I,J)
      FI1J = F(I+1,J)
      VI1J = V(I+1,J)
      FIJ1 = F(I,J+1)
      VIJ1 = V(I,J+1)
      A = A + ABS(FIJ-FI1J-VIJ+VI1J)
     1      + ABS(FIJ-FIJ1-VIJ+VIJ1)
      B = B + MAX( ABS(FIJ-FI1J), ABS(VIJ-VI1J))
     1      + MAX( ABS(FIJ-FIJ1), ABS(VIJ-VIJ1))
   10 CONTINUE
*
      IF (B.NE.0.) S1 = (A/B) * 100.
*
*--------------------------------------------------------------------------
*
      RETURN
      END
