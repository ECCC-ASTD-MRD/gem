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
***S/P IDXES - CALCULATES IDX, JDX NECESSARY TO INTERPOLATE
*              TO A POINT I.DI, J.DJ IN GRID CONTAINING
*              LJ LATITUDES CIRCLES.


      SUBROUTINE IDXES (IDX,JDX,I,J,LI,LJ,KASE)
*
*AUTEUR   - D. ROBERTSON
*
*LANGAGE  - FORTRAN
*
*OBJET(IDXES)
*         - CALCULATES IDX, JDX NECESSARY TO INTERPOLATE
*           TO A POINT I.DI, J.DI IN GRID CONTAINING
*           LJ LATITUDES CIRCLES, EACH CIRCLE HOLDING
*           LI POINTS.
*
*ALGORITHME
*         - ON EXIT THE 4 VALUES IN JDX LOCATE THE 4 LATITUDE
*           CIRCLES TO BE USED, FOR A PARTICULAR J VALUE, JPT SAY,
*           THE 4 POINTS WITH INDICES (IDX(II,JPT),JDX(JJ)),
*           FOR II= 1 TO 4, ARE THE ONES TO INDEXED AS
*           (IDX(IPT,JJ),JDX(JJ)) ARE ON A LONGITUDE LINE, AND
*           CAN BE USED TO INTERPOLATE TO THE FINAL LATITUDE.
*           FOR HEMISPHERIC GRIDS, KASE=NORTH OR SOUTH, POINTS
*           BELOW OR ABOVE THE EQUATOR THE J VALUE IS REFLECTED
*           BACK ACROSS THE EQUATOR, AND GIVEN A -VE SIGN TO INDICATE
*           THIS.  I.E. IF JDX(JJ) IS < 0, THE SYMMETRICALLY LOCATED
*           ROW HAS J = -JDX(JJ).
*
*LIBRAIRIES
*         - SOURCE   RMNSOURCELIB,ID=RMNP     DECK=IDXES
*         - OBJET    RMNLIB,ID=RMNP
*
*APPEL    - CALL IDXES(IDX,JDX,I,J,LI,LJ,KASE)
*
*ARGUMENTS
*   OUT   - IDX  - ARRAY CONTAINING INDICES OF POINTS ON A LONGITUDE
*                  LINE TO INTERPOLATE TO THE FINAL LATUTUDE.
*   OUT   - JDX  - VALUES THAT LOCATE THE LATITUDE CIRCLE
*                  TO BE USED.
*   IN    - I    - INTEGER PART OF FRACTIONAL INDEX COORDINATE OF TARGET POINT
*         - J    - INTEGER PART OF FRACTIONAL INDEX COORDIANTE OF TARGET POINT
*         - LI   - NUMBER OF POINTS ON EACH LATITUDE CIRCLE.
*         - LJ   - NUMBER OF LATITUDES CIRCLES.
*         - KASE - = 0, GLOBAL
*                  = 1, NORTH
*                  = 2, SOUTH
*
*----------------------------------------------------------------------------
*
      INTEGER IDX (4,4) , JDX(4)
       LOGICAL RANGE
#include "kases.cdk"
*
      LIO2 = LI/2
*
*-----------------------------------------------------------------------
*              TEST IF (I,J) IN RANGE
      LJJ = LJ
      IF(KASE.EQ.NORTH) LJJ = 2*LJ
       IF(I.GT.0.AND.I.LE.LI.AND.J.GE.0.AND.J.LE.LJJ) GO TO 80
                    WRITE(6,630)  I,  J
*-----------------------------------------------------------------------
*
   80 RANGE = (I.GT.1.AND.I.LT.LI-1.AND.J.GT.1.AND.J.LT.LJ-1)
                         ISGN = +1
      IF (KASE.EQ.SOUTH) ISGN = -1
*
                         JSGN = 1
      IF (KASE.EQ.NORTH) JSGN = -1
*
      DO 95 JJ = 1 , 4
        JDX(JJ) = J-2+JJ
        DO 85 II = 1 , 4
           IDX(II,JJ) = I-2+II
   85      CONTINUE
   95   CONTINUE
*-----------------------------------------------------------------------
*
*              INTERPOLATION DO NOT GO ACROSS THE POLE
*                   OR GREENWICH -RETURNM-
*
                    IF(RANGE) RETURN
*-----------------------------------------------------------------------
*
* CORRECT VALUES IF ONE HAS GONE ACROSS POLE, THEN REDUCE MODULO LI
*
      DO  500 JJ = 1 , 4
      IF  ((JDX(JJ).GT.0 .AND. JDX(JJ).LE.LJ))
     $    GO TO 400
      IF (KASE .EQ. NORTH .AND. JDX(JJ).GT.0)
     $    GO TO 400
      IF (KASE.EQ.SOUTH .AND. JDX(JJ).LE.LJ)
     $    GO TO 400
*
* MOVE I VALUES HALFWAY AROUND
*
        DO 300 II = 1 , 4
           IDX(II,JJ) = IDX(II,JJ)+LIO2
  300      CONTINUE
  400   DO 450 II = 1 , 4
           IF (IDX(II,JJ) .LT.1)
     $         IDX(II,JJ) = IDX(II,JJ)+LI
           IF (IDX(II,JJ) .GT.LI)
     $         IDX(II,JJ) = IDX(II,JJ)-LI
  450     CONTINUE
*
        IF (JDX(JJ).LT.1)
     $      JDX(JJ) = ISIGN(1-JDX(JJ),ISGN)
        IF (JDX(JJ).GT.LJ)
     $            JDX(JJ) = ISIGN(2*LJ+1-JDX(JJ),JSGN)
*
  500 CONTINUE
      RETURN
*
  630 FORMAT(/" WRONG VALUE OF I OR J IN SBRT IDXES",2I5)
      END
