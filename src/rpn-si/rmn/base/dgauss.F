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
*
***S/P DGAUSS - CALCULATES THE ZEROES OF THE ORDINARY LEGENDRE
*               POLYNOMIAL OF ORDER N,  I.E. DEFINE GAUSSIAN GRID


      SUBROUTINE DGAUSS (N,ROOTS,KASE)
      REAL ROOTS(*)
*
*AUTEUR   - D. ROBERTSON
*
*LANGAGE  - FORTRAN
*
*OBJET(DGAUSS)
*         - CALCULATES THE ZEROES OF THE ORDINARY LEGENDRE
*           POLYNOMIAL OF ORDER N,  I.E. DEFINE GAUSSIAN GRID
*
*ALGORITHME
*         - THE POSITIVE ROOTS ARE APPROXIMATED BY THE BEST
*           ASYMPTOTIC FORMULA AVAILABLE TO THE AUTHOR, FOUND IN
*           ABRAMOWITZ AND STEGUN "HANDBOOK OF MATHEMATICAL FUNCTIONS".
*           CHAPTER 22 FORMULA 22.16.6.
*           NEWTON'S METHOD IS USED TO REFINE THE GUESS TO PRECISION
*           DEFINED BY THE CONSTANT TOL.  SINCE THE ROOTS ARE OF ORDER
*           OF MAGNITUDE UNITY, ABSOLUTE PRECISION IS ADEQUATE, RATHER
*           THAN A RELATIVE TEST.
*           A STANDARD IDENTITY IS USED TO DETERMINE THE DERIVATIVE OF
*           THE POLYNOMIAL IN TERMS OF THE VALUES OF P(N;X),P(N-1;X).
*           (X**2-1.0)*(DP/DX)=N*(X*P(N;,X)-P(N-1;X)).
*           SEE ABRAMOWITZ AND STEGUN FORMULA 22.8.5
*           NOTE THAT IN CONTRAST TO OTHER FORMULAS THIS REQUIRES ONLY
*           2 EVALUATIONS OF A LEGENDRE POLYNOMIAL PER ITERATION.
*           NOTE THAT THE COORDINATE USED IS CONVENTIONALLY REFERRED TO
*           AS MU=COS(THETA), RUNNING FROM +1 TO -1, FOR THETA FROM 0 TO
*           PIE. THE NEGATIVE ROOTS ARE  FILLED BY SYMMETRY.
*           FOR KASE=GLOBAL, ALL N ROOTS ARE FOUND, WHILE FOR
*           DASE=NORTH/SOUTH ONLY THE +VE/-VE ROOTS ARE FOUND,
*           (INCLUDING 0 IF N IS ODD)  I.E. N/2+MOD(N,2) ROOTS.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=DGAUSS
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL DGAUSS(N,ROOTS,KASE)
*
*ARGUMENTS
*   IN    - N     - ORDER OF THE POLYNOMIALS
*   OUT   - ROOTS - ARRAY CONTAINING THE ZEROES OF THE
*                   ORDINARY LEGENDRE POLYNOMIALS
*   IN    - KASE  - =0, GLOBAL
*                   =1, NORTH
*                   =2, SOUTH
*
*MODULES APPELES
*         - ORDLEG
*
*--------------------------------------------------------------------------
*
* THE ANSWERS ARE RETURNED IN ROOTS.
      REAL NORMN,NORMNM
#include "kases.cdk"
#if defined (ALL64)
      DATA TOL /1.0E-13/
#else
      DATA TOL /1.0E-6/
#endif
#include "pi.cdk"
*
*  ORDLEG RETURNS POLYNOMIALS NORMALIZED TO UNIT INTEGRAL.
*  NORMN,NORMNMN RESTORE THE CONVENTION NORMALIZATION, P(N;1.0)=1.0.
*
      NORMN =SQRT(2.0/(2.0*N+1.0))
      NORMNM=SQRT(2.0/(2.0*N-1.0))
      L = N/2
*
* CALCULATE ASYMPTOTIC APPROXIMATION
*
      DO 14 I=1,L
         IF(KASE.NE.SOUTH) J = I
         IF(KASE.EQ.SOUTH) J=I+L+MOD(N,2)
         T = (4*J-1)*PIE/FLOAT(4*N+2)
         IF(KASE.NE.SOUTH) IRT = I
         IF(KASE.EQ.SOUTH) IRT = I + MOD(N,2)
         ROOTS(IRT)=COS(T+1.0/(8.0*FLOAT(N**2)*TAN(T)))
   14    CONTINUE
*
      DO 15 I=1,L
*
*    REPEAT 1 NEWTON ITERATION
* **  BEGIN
*
    6       CALL ORDLEG(G,ROOTS(I),N)
            CALL ORDLEG(GM,ROOTS(I),N-1)
            PN = NORMN*G
            PNM= NORMNM*GM
*
* **         GUESS(K+1)=GUESS(K)-P/(DP/DX)
*
            RDPDX = (ROOTS(I)**2-1.0)/(N*(ROOTS(I)*PN-PNM))
            DELTA = -PN*RDPDX
            ROOTS(I) = ROOTS(I)+DELTA
            IF(ABS(DELTA).GT.TOL) GO TO 6
*
* **  END UNTIL ABS(DELTA).LE.TOL
*
               ROOTS(N+1-I) = -ROOTS(I)
   15    CONTINUE
*
      IF(MOD(N,2).EQ.0) RETURN
      IF(KASE.NE.SOUTH) IRT = L+1
      IF(KASE.EQ.SOUTH) IRT = 1
      ROOTS(IRT) = 0.0
      RETURN
      END


