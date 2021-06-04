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
***FONCTION SLL - INTERPOLATION A PARTIR D'UNE GRILLE LAT-LON


      FUNCTION SLL(LAMBDA,PHI,F,LI,LJ,XLA0,XLO0,DLAT,DLON)
      REAL F(LI,LJ), T(4)
*
*AUTEUR   - M. VALIN, DRPN
*
*REVISION 001  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGAGE  - FORTRAN
*
*OBJET(SLL)
*         - EXTRAIRE LA VALEUR DE LA FONCTION F AU POINT DE COORDONNEES
*           LAMBDA(LONGITUDE) ET PHI(COLATITUDE) DANS UNE GRILLE
*           LAT-LON.
*
*APPEL    - R=SLL(LAMBDA,PHI,F,LI,LJ,XLA0,XLO0,DLAT,DLON)
*
*LIBRAIRIES
*         - RMNSOURCELIB,ID=RMNP     DECK=SLL
*         - RMNLIB,ID=RMNP
*
*ARGUMENTS
*   OUT   - SLL    - VALEUR DE LA FONCTION AU POINT DEMANDE
*   IN    - LAMBDA - LONGITUDE DU POINT  (REEL)
*   IN    - PHI    - COLATITUDE(LATITUDE-90.0) DU POINT (REEL)
*   IN    - F      - CHAMP SOURCE  (LI,LJ)
*   IN    - LI     - PREMIERE DIMENSION DE LA GRILLE FF (PS)
*   IN    - LJ     - DEUXIEME DIMENSION DE LA GRILLE FF (PS)
*   IN    - XLA0   - LATITUDE DU COIN INFERIEUR GAUCHE (DEGRES).
*   IN    - XLO0    - LONGITUDE DU COIN INFERIEUR GAUCHE (DEGRES).
*   IN    - DLAT   - ESPACEMENT EN LATITUDE (DEGRES).
*   IN    - DLON   - ESPACEMENT EN LONGITUDE (DEGRES).
*
*--------------------------------------------------------------------------
*
#include "ordr.cdk"
      REAL LAMBDA
*      DATA ORDINT /3/
#include "lagrange.cdk"
#include "cubic.cdk"


*
      ORDINT = 3
      X = 1 + (LAMBDA-XLO0) / DLON
      Y = 1 + (90.-PHI-XLA0) / DLAT
      IF(ORDINT.EQ.0) THEN
         I = NINT(X)
         J = NINT(Y)
         SLL = F(MAX(1,MIN(I,LI)),MAX(1,MIN(J,LJ)))
         RETURN
      ELSE IF(ORDINT.EQ.1) THEN
         I = X
         I = MAX(1,MIN(I,LI-1))
         J = Y
         J = MAX(1,MIN(J,LJ-1))
         SLL = ZLIN( ZLIN(F(I,J),F(I+1,J),X-I),
     %               ZLIN(F(I,J+1),F(I+1,J+1),X-I),
     %               Y-J)
         RETURN
      ENDIF
      I = X
      J = Y
      I = MIN0(MAX0(I,2),LI-2)
      J = MIN0(MAX0(J,2),LJ-2)
      DX = X - I
      DY = Y - J


      IF( X.LT.2.0 .OR. X.GT.LI-1) GOTO 10
*
* A L'INTERIEUR DE LA GRILLE DANS LA DIRECTION X
*
        T(1)=CUBIC(F(I-1,J-1),F(I  ,J-1),F(I+1,J-1),F(I+2,J-1),DX)
        T(2)=CUBIC(F(I-1,J  ),F(I  ,J  ),F(I+1,J  ),F(I+2,J  ),DX)
        T(3)=CUBIC(F(I-1,J+1),F(I  ,J+1),F(I+1,J+1),F(I+2,J+1),DX)
        T(4)=CUBIC(F(I-1,J+2),F(I  ,J+2),F(I+1,J+2),F(I+2,J+2),DX)
        IF( Y.LT.2.0 .OR. Y.GT.LJ-1 ) GOTO 5
*
* A L'INTERIEUR DE LA GRILLE DANS LA DIRECTION Y.
*
          SLL=CUBIC(T(1),T(2),T(3),T(4),DY)
          RETURN
*
* A L'EXTERIEUR DE LA GRILLE DANS LA DIRECTION Y.
*
5         SLL=POLY(T(1),T(2),T(3),T(4),Y,FLOAT(J-1),FLOAT(J),
     X          FLOAT(J+1),FLOAT(J+2))
          RETURN


10    IF( Y.LT.2.0 .OR. Y.GT.LJ-1 ) GOTO 15
*
* A L'INTERIEUR DE LA GRILLE DANS LA DIRECTION Y
*
        T(1)=CUBIC(F(I-1,J-1),F(I-1,J  ),F(I-1,J+1),F(I-1,J+2),DY)
        T(2)=CUBIC(F(I  ,J-1),F(I  ,J  ),F(I  ,J+1),F(I  ,J+2),DY)
        T(3)=CUBIC(F(I+1,J-1),F(I+1,J  ),F(I+1,J+1),F(I+1,J+2),DY)
        T(4)=CUBIC(F(I+2,J-1),F(I+2,J  ),F(I+2,J+1),F(I+2,J+2),DY)
        GOTO 20
*
* PAS A L'INTERIEUR DE LA GRILLE DANS LA DIRECTION Y
*
15      T(1)=POLY(F(I-1,J-1),F(I-1,J  ),F(I-1,J+1),F(I-1,J+2),Y,
     X            FLOAT(J-1),FLOAT(J),FLOAT(J+1),FLOAT(J+2))
        T(2)=POLY(F(I  ,J-1),F(I  ,J  ),F(I  ,J+1),F(I  ,J+2),Y,
     X            FLOAT(J-1),FLOAT(J),FLOAT(J+1),FLOAT(J+2))
        T(3)=POLY(F(I+1,J-1),F(I+1,J  ),F(I+1,J+1),F(I+1,J+2),Y,
     X            FLOAT(J-1),FLOAT(J),FLOAT(J+1),FLOAT(J+2))
        T(4)=POLY(F(I+2,J-1),F(I+2,J  ),F(I+2,J+1),F(I+2,J+2),Y,
     X            FLOAT(J-1),FLOAT(J),FLOAT(J+1),FLOAT(J+2))
*
* PAS A L'INTERIEUR DE LA GRILLE DANS LA DIRECTION X
*
20      SLL=POLY(T(1),T(2),T(3),T(4),X,
     X           FLOAT(I-1),FLOAT(I),FLOAT(I+1),FLOAT(I+2))


      RETURN
      END
