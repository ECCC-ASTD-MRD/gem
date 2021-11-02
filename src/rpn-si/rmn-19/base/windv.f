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
***S/R WINDV - COMPUTES WEIGHT R.M.S. ERROR
*
      SUBROUTINE WINDV(RMS,FU,FV,VU,VV,W,NI,NJ,IW1,IW2,NW1,NW2)
      implicit none
      integer :: NI,NJ,IW1,IW2,NW1,NW2
      REAL FU(NI,NJ),FV(NI,NJ),VU(NI,NJ),VV(NI,NJ),W(NI,NJ)
      REAL :: RMS

*
*AUTHOR   - Y.R. BOURASSA  -  APR 75
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*REVISION 002   C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*REVISION 003   M. Valin      -  JUIN 2015  implicit none + utilisation de Real*8 pour les calculs
*
*LANGUAGE - fortran
*
*OBJECT(WINDV)
*         - GIVEN A FORECAST, A VERIFICATION, A WEIGHT FIELD AND THE
*           CO-ORDINATES OF THE WINDOW, COMPUTES THE WEIGHTED
*           R.M.S. ERROR.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=WINDV
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL WINDV(RMN,FU,FV,VU,VV,W,NI,NJ,IW1,IW2,NW1,NW2)
*
*ARGUMENTS
*   OUT   - RMS - WEIGHTED R.M.S. ERROR OF VECTORS
*   IN    - FU  - X-COMPONENT OF FORECASTED WIND
*         - FV  - Y-COMPONENT OF FORECASTED WIND
*         - VU  - X-COMPONENT OF VERIFYING WIND
*         - VV  - Y-COMPONENT OF VERIFYING WIND
*         - W   - WEIGHT FIELD (SET TO 1'S IF WEIGHTING NOT DESIRED)
*         - NI  - X-DIMENSION
*         - NJ  - Y-DIMENSION
*         - IW1 - BOTTOM LEFT CORNER OF X-COORDINATE OF WINDOW
*         - IW2 - BOTTOM LEFT CORNER OF Y-COORDINATE OF WINDOW
*         - NW1 - UPPER RIGHT CORNER OF X-COORDINATE OF WINDOW
*         - NW2 - UPPER RIGHT CORNER OF Y-COORDINATE OF WINDOW
*
*NOTES    - IF THE CO-ORDINATES OF THE WINDOW ARE OUTSIDE THE LIMIT OF
*           THE ARRAYS OR IF THE TOTAL WEIGHT IS EQUAL TO ZERO, RMS IS SET
*           TO 99999.
*
*------------------------------------------------------------------------------
*
*
      real *8 :: FUMVU, FVMVV, WIJ, Y, FTW
      integer :: I, J
      RMS  = 99999.
*
*  TEST TO SEE IF THE WINDOW CO-ORDINATES ARE INSIDE THE
*  LIMITS OF THE ARRAYS DIMENSIONED (NI*NJ)
*
      IF (NI.LT.1.OR.NI.LT.IW1.OR.NI.LT.NW1) RETURN
      IF (NJ.LT.1.OR.NJ.LT.IW2.OR.NJ.LT.NW2) RETURN
      IF (NW1.LT.IW1.OR.NW2.LT.IW2) RETURN
*
      FTW = 0.
      Y = 0.
*
        DO 10 J=IW2,NW2
        DO 10 I=IW1,NW1
        FUMVU = FU(I,J) - VU(I,J)
        FVMVV = FV(I,J) - VV(I,J)
        WIJ   = W(I,J)
        Y = Y + FUMVU * FUMVU * WIJ + FVMVV * FVMVV * WIJ
        FTW = FTW + WIJ
   10   CONTINUE
*
      IF (FTW.EQ.0.) RETURN
      RMS = SQRT(Y/FTW)
*
*---------------------------------------------------------------------------
*
      RETURN
      END
