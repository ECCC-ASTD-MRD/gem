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
***S/R MWERMS - COMPUTES THE MEAN WEIGHTED ERROR AND THE WEIGHTED
*               R.M.S. DEVIATION FROM THE MEAN
*
      SUBROUTINE MWERMS(BIAS,STD,F,V,W,NI,NJ,IW1,IW2,NW1,NW2)
      implicit none
      integer :: NI,NJ,IW1,IW2,NW1,NW2
      REAL F(NI,NJ), V(NI,NJ), W(NI,NJ)
      REAL :: BIAS, STD

*
*AUTHOR   - Y.R. BOURASSA  -  MAR 75
*
*REVISION 001   C. THIBEAULT  NOV 79  DOCUMENTATION
*REVISION 002   C. THIBEAULT  MAR 83  CONVERSION AU CODE CRAY
*REVISION 003   M. Valin      JUIN 2015  implicit none + utilisation de Real*8 pour les calculs
*
*LANGUAGE - fortran
*
*OBJECT(MWERMS)
*         - GIVEN A FORECAST, A VERIFICATION, A WEIGHT FIELD AND THE
*           CO-ORDINATES OF THE WINDOW, COMPUTES THE MEAN WEIGHTED ERROR
*           AND THE WEIGHTED R.M.S. DEVIATION FROM THE MEAN.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=MWERMS
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL MWERMS(BIAS,STD,F,V,W,NI,NJ,IW1,IW2,NW1,NW2)
*
*ARGUMENTS
*   OUT   - BIAS - MEAN WEIGHTED ERROR
*         - STD  - R.M.S. DEVIATION FROM THE MEAN
*   IN    - F    - FORECAST FIELD
*         - V    - VERIFYING FIELD
*         - W    - WEIGHT FIELD (SET TO 1'S IF WEIGHTING IS NOT DESIRED)
*         - NI   - X-DIMENSION
*         - NJ   - Y-DIMENSION
*         - IW1  - BOTTOM LEFT CORNER OF X-COORDINATE OF WINDOW
*         - IW2  - BOTTOM LEFT CORNER OF Y-COORDINATE OF WINDOW
*         - NW1  - UPPER RIGHT CORNER OF X-COORDINATE OF WINDOW
*         - NW2  - UPPER RIGHT CORNER OF Y-COORDINATE OF WINDOW
*
*NOTES    - IF THE CO-ORDINATES OF THE WINDOW ARE OUTSIDE THE LIMIT OF
*           THE ARRAYS OR THE TOTAL WEIGHT IS EQUAL TO ZERO, BIAS AND STD
*           ARE SET TO 99999.
*
*------------------------------------------------------------------------------
*
*
*--------------------------------------------------------------------------
*
      REAL *8 BIAS2, x, y, FTW, FMV, WIJ
      integer :: I, J

      BIAS = 99999.
      STD = 99999.
*
*  TEST TO SEE IF THE WINDOW CO-ORDINATES ARE INSIDE THE
*  LIMITS OF THE ARRAYS DIMENSIONED (NI*NJ)
*
      IF (NI.LT.1.OR.NI.LT.IW1.OR.NI.LT.NW1) RETURN
      IF (NJ.LT.1.OR.NJ.LT.IW2.OR.NJ.LT.NW2) RETURN
      IF (NW1.LT.IW1.OR.NW2.LT.IW2) RETURN
*
      FTW = 0.
      X = 0.
      Y = 0.
*
        DO 10 J=IW2,NW2
        DO 10 I=IW1,NW1
        WIJ = W(I,J)
        FMV = F(I,J)
        FMV = FMV-V(I,J)
        X = X + FMV*WIJ
        Y = Y + FMV*FMV*WIJ
        FTW = FTW + WIJ
   10   CONTINUE
*
      IF (FTW.EQ.0.) RETURN
      BIAS2 = X / FTW
      Y = Y / FTW
      STD = SQRT(max(0.0,Y-BIAS2*BIAS2))
      BIAS = BIAS2
*
*---------------------------------------------------------------------------
*
      RETURN
      END
