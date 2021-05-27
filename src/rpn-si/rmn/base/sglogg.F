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
***FONCTION SGLOGG - INTERPOLATION A PARTIR D'UNE GRILLE GAUSSIENNE


      FUNCTION SGLOGG(G,LI,LJ,ROOTS,LAMBDA,PHI)
      REAL G(LI,LJ)
*
*AUTEUR   - D.ROBERTSON RPN 1981
*
*REVISION 001   M. VALIN  -  1981  SYMETRIES ET COMMENTAIRES
*REVISION 002   C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGAGE  - FORTRAN
*
*OBJET(SGLOGG)
*         - TROUVER LA VALEUR D'UNE FONCTION EN UN POINT QUELCONQUE D'UNE
*           GRILLE GAUSSIENNE
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=SGLOGG
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - R=SGLOGG(G,LI,LJ,ROOTS,LAMBDA,PHI)
*
*ARGUMENTS
*   OUT   - SGLOGG - VALEUR INTERPOLEE AU POINT DEMANDE
*   IN    - G      - CHAMP SOURCE
*   IN    - LI     - PREMIERE DIMENSION DE G
*   IN    - LJ     - DEUXIEME DIMENSION DE G
*   IN    - ROOTS  - TABLE SERVANT A CALCULER LES RACINES DES POLYNOMES DE
*                    LEGENDRE SEVANT A DEFINIR LA GRILLE GAUSSIENNE. SI
*                    ROOTS CONTIENT DEJA CES RACINES (I.E. SI ROOTS(1) EST
*                    POSITIF ET INFERIEUR A 1.0) LES VALEURS DE ROOTS SERONT
*                    UTILISEES DIRECTEMENT.
*   IN    - LAMBDA - LONGITUDE DU POINT DESIRE
*   IN    - PHI    - COLATITUDE (90.0-LATITUDE) DU POINT DESIRE
*
*MODULES  - DGAUSS,IDXES
*
*----------------------------------------------------------------
*
#include "ordr.cdk"
      REAL LAMBDA,MU
      REAL ROOTS(LJ) , T(4) , P(4)
      INTEGER IDX(4,4), JDX(4)
#include "kases.cdk"
#include "pi.cdk"
*      DATA ORDINT /3/
#include "lagrange.cdk"
*
#include "cubic.cdk"
*
      ORDINT = 3
      KASE = GLOBAL
      SIGN = 1.0
      NPOLY = LJ
      GO TO 99
*
*     SCALAR SYMMETRIC NORTHERN HEMISPHERE
*
      ENTRY SSNHGG(G,LI,LJ,ROOTS,LAMBDA,PHI)
      KASE = NORTH
      SIGN = 1.0
      NPOLY = 2*LJ
      GO TO 99
*
*     SCALAR ANTISYMMETRIC NORTHERN HEMISPHERE
*
      ENTRY SANHGG(G,LI,LJ,ROOTS,LAMBDA,PHI)
      KASE = NORTH
      NPOLY = 2*LJ
      SIGN = -1.0
*
* START OF ACTUAL INTERPOLATION
* TAKE ABS(MU) IF HEMISPHERIC
*
   99 CONTINUE
*
      MU = COS(PHI*DGTORD)
      SIGN2 = 1.0
      IF(MU.LT.0) SIGN2 = -SIGN2
      IF (KASE.EQ.NORTH) MU = ABS(MU)
*
* DEFINE GAUSSIAN GRID IF NECESSARY
*
      IF (LEGVAR(ROOTS(1)).NE.0)
     1  CALL DGAUSS(NPOLY,ROOTS,KASE)
      IF(ROOTS(1).GT.1.0)
     1  CALL DGAUSS(NPOLY,ROOTS,KASE)
*
* THERE ARE 5 CASES TO CONSIDER DEPENDING ON MU
* 2 NEAR THE N POLE, 1 INTERIOUR, AND 2 NEAR THE S POLE.
*
      IF(MU.GT.ROOTS(1)) ICASE = 1
      IF(MU.LE.ROOTS(1) .AND. MU.GT.ROOTS(2)) ICASE=2
      IF(MU.LE.ROOTS(2).AND.MU.GT.ROOTS(NPOLY-1)) ICASE = 3
      IF(MU.LE.ROOTS(NPOLY-1).AND.MU.GT.ROOTS(NPOLY)) ICASE = 4
      IF(MU.LE.ROOTS(NPOLY)) ICASE = 5
      GO TO (1,2,3,4,5),ICASE
*
*     BETWEEN N POLE AND FIRST ROW
*
    1 JJ = 0
      A = 2.0 -ROOTS(2)
      B = 2.0 -ROOTS(1)
      C = ROOTS(1)
      D = ROOTS(2)
      GO TO 30
*
*     BETWEEN TOP TWO ROWS
*
    2 JJ = 1
      A = 2.0 -ROOTS(1)
      B = ROOTS(1)
      C = ROOTS(2)
      D = ROOTS(3)
      GO TO 30
*
*     FAR ENOUGH FROM POLES
*
    3 FRAC = 0.5*(1.0-MU)
      JJ = MIN0(NPOLY-1,MAX0(2,INT(FRAC*NPOLY)))
*
* ONE MUST LOCATE MU, I.E. ROOTS(JJ).GE.MU.GT.ROOTS(JJ+1)
* THE ABOVE IS A SUITABLE 1ST GUESS.
* WHILE ROOTS(JJ).LT.MU DO JJ=JJ-1
*
   10 IF(ROOTS(JJ).GE.MU) GO TO 15
         JJ=JJ-1
         GO TO 10
*
* WHILE ROOTS(JJ+1).GE.MU DO JJ=JJ+1
*
   15 IF(ROOTS(JJ+1).LT.MU) GO TO 17
         JJ=JJ+1
         GO TO 15
   17 CONTINUE
*
* ROOTS(JJ).GE.MU.GT.ROOTS(JJ+1)
*
      A = ROOTS(JJ-1)
      B = ROOTS(JJ)
      C = ROOTS(JJ+1)
      D = ROOTS(JJ+2)
      GO TO 30
*
*     BETWEEN LAST TWO ROWS
*
    4 JJ = NPOLY-1
      A = ROOTS(NPOLY-2)
      B = ROOTS(NPOLY-1)
      C = ROOTS(NPOLY)
      IF(KASE.EQ.GLOBAL) D = -2.0-ROOTS(NPOLY)
      IF(KASE.EQ.NORTH)  D = -ROOTS(NPOLY)
      GO TO 30
*
*     BETWEEN LAST ROW AND S POLE
*
    5 JJ = NPOLY
      A = ROOTS(NPOLY-1)
      B = ROOTS(NPOLY)
      IF(KASE.EQ.GLOBAL) C = -2.0 - ROOTS(NPOLY)
      IF(KASE.EQ.NORTH)  C = -ROOTS(NPOLY)
      IF(KASE.EQ.GLOBAL) D = -2.0 - ROOTS(NPOLY-1)
      IF(KASE.EQ.NORTH)  D = -ROOTS(NPOLY-1)
      GO TO 30
   30 CONTINUE
*
*     WE KNOW WHERE WE ARE, NOW DO THE WORK
*
      RDLAM = FLOAT(LI)/360.0
      XX = 1.0 + LAMBDA*RDLAM
      II = XX
      DX = XX-II
      CALL IDXES(IDX,JDX,II,JJ,LI,LJ,KASE)
      IF(ORDINT.EQ.0) THEN
         II = NINT(DX+2.0)
         JJ = NINT((MU-B)/(C-B)+2.0)
         SGLOGG = G(IDX(II,JJ),ABS(JDX(JJ)))
         IF(JDX(JJ).LT.0) SGLOGG = SGLOGG*SIGN
      ELSE IF(ORDINT.EQ.1) THEN
         DO 90 JJ = 2 , 3
            DO 80 II = 2 , 3
80             P(II) = G(IDX(II,JJ),ABS(JDX(JJ)))
            IF(JDX(JJ).LT.0) THEN
               T(JJ) = SIGN * ZLIN(P(2),P(3),DX)
            ELSE
               T(JJ) = ZLIN(P(2),P(3),DX)
            ENDIF
90       CONTINUE
         SGLOGG = ZLIN(T(2),T(3),(MU-B)/(C-B))
      ELSE
*
*     FOUR INTERPOLATIONS WITH CONSTANT GRID LENGTH IN THE E-W
*     DIRECTION, FOLLOWED BY ONE IN THE N-S DIRECTION
*
      DO 47 J=1,4
         DO 45 K = 1 , 4
            P(K) = G(IDX(K,J),IABS(JDX(J)))
            IF(JDX(J).LT.0) P(K) = SIGN*P(K)
   45      CONTINUE
         T(J)=CUBIC(P(1),P(2),P(3),P(4),DX)
   47    CONTINUE
*
      SGLOGG = POLY(T(1),T(2),T(3),T(4),MU,A,B,C,D)
      ENDIF
*
*     APPLY SYMMETRIES IF NECESSARY
*
      IF(KASE.EQ.NORTH .AND. SIGN.LT.0) SGLOGG = SIGN2*SGLOGG
      RETURN
      END
