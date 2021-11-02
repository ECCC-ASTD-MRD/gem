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


***FONCTION SPS - INTERPOLATION A PARTIR D'UNE GRILLE POLAIRE STEREOGRAPHIQUE
*
      FUNCTION SPS(LAMBDA,PHI,F,LI,LJ,PI,PJ,D60,DGRW,HEM)
      REAL F(LI,LJ), T(4)
*
*AUTEUR   - M. VALIN, DRPN
*
*REVISION 001  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGAGE  - FORTRAN
*
*OBJET(SPS)
*         - EXTRAIRE LA VALEUR DE LA FONCTION F AU POINT DE COORDONNEES
*           LAMBDA(LONGITUDE) ET PHI(COLATITUDE) DANS UNE GRILLE
*           POLAIRE STEREOGRAPHIQUE.
*
*APPEL    - R=SPS(LAMBDA,PHI,F,LI,LJ,PI,PJ,D60,DGRW,HEM)
*
*LIBRAIRIES
*         - RMNSOURCELIB,ID=RMNP     DECK=SPS
*         - RMNLIB,ID=RMNP
*
*ARGUMENTS
*   OUT   - SPS    - VALEUR DE LA FONCTION AU POINT DEMANDE
*   IN    - LAMBDA - LONGITUDE DU POINT  (REEL)
*   IN    - PHI    - COLATITUDE(LATITUDE-90.0) DU POINT (REEL)
*   IN    - F      - CHAMP SOURCE  (LI,LJ)
*   IN    - LI     - PREMIERE DIMENSION DE LA GRILLE F (PS)
*   IN    - LJ     - DEUXIEME DIMENSION DE LA GRILLE F (PS)
*   IN    - PI     - COORDONNEE X DU POLE  (REEL)
*   IN    - PJ     - COORDONNEE Y DU POLE  (REEL)
*   IN    - D60    - DISTANCE EN METRES ENTRE LES POINTS DE LA GRILLE
*                    A 60 DEGRES DE LATITUDE
*   IN    - DGRW   - ANGLE ENTRE L'AXE X ET LE MERIDIEN DE GRENWICH
*   IN    - HEM    - CODE D'HEMISPHERE  (1=NORD , 2=SUD)
*
*MODULES  - XYFLL
*
*--------------------------------------------------------------------------
*
#include "ordr.cdk"
      REAL LAMBDA
      INTEGER HEM
*      DATA ORDINT /3/
#include "lagrange.cdk"
#include "cubic.cdk"


*
      ORDINT = 3
      CALL XYFLL(X,Y,90.0-PHI,LAMBDA,D60,DGRW,HEM)
1     X = X + PI
      Y = Y + PJ
      IF(ORDINT.EQ.0) THEN
         I = NINT(X)
         J = NINT(Y)
         SPS = F(MAX(1,MIN(I,LI)),MAX(1,MIN(J,LJ)))
         RETURN
      ELSE IF(ORDINT.EQ.1) THEN
         I = X
         I = MAX(1,MIN(I,LI-1))
         J = Y
         J = MAX(1,MIN(J,LJ-1))
         SPS = ZLIN( ZLIN(F(I,J),F(I+1,J),X-I),
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
          SPS=CUBIC(T(1),T(2),T(3),T(4),DY)
          RETURN
*
* A L'EXTERIEUR DE LA GRILLE DANS LA DIRECTION Y.
*
5         SPS=POLY(T(1),T(2),T(3),T(4),Y,FLOAT(J-1),FLOAT(J),
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
20      SPS=POLY(T(1),T(2),T(3),T(4),X,
     X           FLOAT(I-1),FLOAT(I),FLOAT(I+1),FLOAT(I+2))


      RETURN
      ENTRY SLC(LAMBDA,PHI,F,LI,LJ,PI,PJ,D60,DGRW,HEM,PHI1,PHI2)
      CALL XYFLCC(X,Y,90.-PHI,LAMBDA,D60,DGRW,1.,PHI1,PHI2)
      GOTO 1
      END
