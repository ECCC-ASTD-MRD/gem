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
***S/P VGLOGG - RETURN 2  SCALAR QUANTITIES FROM U AND V FIELD
*             AT POINT CLAMDA (LONGITUDE) AND PHI(LATITUDE)
*
      SUBROUTINE VGLOGG(SPEED , PSI ,CLAMDA, PHI ,U,V,LI,LJ ,ROOTS)
      REAL U(LI,LJ), V(LI,LJ)
*
*AUTEUR   - M. VALIN  -  MODIFICATION DE VGLO
*
*REVISION 001  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE - FORTRAN
*
*OBJECT(VGLOGG)
*         - RETURN 2 SCALARS QUANTITIES USING U-V VECTOR FIELD
*         - ON LAT LON COMPUTE SPEED AND PSI WITH A CUBIC
*         - INTERPOLATION ON EAST TO WEST  AND  NORTH TO SOUTH
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=VGLOGG
*         - OBJECT  RMNLIB,ID=RMNP
*
*APPEL    - CALL VGLOGG(SPEED,PSI,CLAMDA,PHI,U,V,LI,LJ,ROOTS)
*
*ARGUMENTS
*   OUT   - SPEED  - WIND SPEED, SAME UNITS AS U,V
*   OUT   - PSI    - ANGLE BETWEEN THE WIND AND GREENWICH ON P.S. GRID
*   IN    - CLAMDA - LONGITUDE OF PT / TO PI OR -PI TO +PI EAST = +
*   IN    - PHI    - COLATITUDE OF U OR V (0 = N.P. 180 = S.P.)
*   IN    - U,V    - ARE UNIFORM LONG LAT GRIDS
*   IN    - LI     - POINTS PER CIRCLE, MUST BE EVEN
*                    FIRST CIRCLE DISPLACE FROM POLE BY OFFSET (90.0/LJ DEG)
*   IN    - LJ     - LATITUDE CIRCLES
*   IN    - ROOTS  - TABLE SERVANT A CALCULER LES RACINES DES POLYNOMES
*                    DE LEGENDRE SERVANT A DEFINIR LA LATITUDE DES
*                    NOEUDS DE LA GRILLE GAUSSIENNE. SI ROOTS CONTIENT
*                    DEJA CES RACINES (I.E. SI ROOTS(1) EST POSITIF ET
*                    INFERIEUR A 1.0) LES VALEURS DE ROOTS SERONT
*                    UTILISES DIRECTEMENT.
*
*MODULES  - NPSPIN,IDXES,ABORT,DGAUSS
*
*MESSAGES - "LAMBDA PHI LI LJ IN SBRT VGLOGG (VWNHGG) WRONG
*            VALUE, ..-..-..-.."
*
*NOTES    - IF SPEED = 0.0, PSI = ARBPSI
*         - SPEED AND PSI ARE GIVEN RELATIVE TO A P.S.
*           GRID WITH DGRW=0.0
*
*-------------------------------------------------------------------------------
*


      REAL MU


**            -VALUE IN JDX(4) LOCATE THE 4 LAT CIRCLES ON WHICH THE
**            -POINT NEEDED LIE
**            -IDX LOCATE  THE POINT ON LATITUDE CIRCLE
      DIMENSION  IDX(4,4) , JDX(4)
*
* SUPS,SVPS HOLD 4 VALUES OF U,V ALONG A MERIDION OF LONGITUDE,
* WHICH ARE USED TO INTERPOLATE A VALUE AT THE FINAL LATITUDE.
* THE 4 VALUES AT THE FINAL LATITUDE ARE HELD IN FUPS,FVPS.
*
      DIMENSION SUPS(4) , SVPS(4) , FUPS(4) , FVPS(4) , ROOTS(LJ)
* POLAR-STEREOGRAPHIC COMPONENTS ARE USED FOR THE INTERPOLATION.


*             -EPSPHI AN ARBITRARY LIMIT


#include "kases.cdk"


      DATA EPSPHI / 1.0E-8 /
#include "pi.cdk"
      DATA ARBPSI /99999.0/


#include "cubic.cdk"
#include "lagrange.cdk"
       KASE = GLOBAL
      NPOLY = LJ
      SIGN0 = 1.0
      GOTO 11
*
* HEMISPHERIC, WALL AT THE EQUATOR
*
      ENTRY VWNHGG(SPEED,PSI,CLAMDA,PHI,U,V,LI,LJ,ROOTS)
      KASE = NORTH
      NPOLY = 2*LJ
      SIGN0 = -1.0
11     LILJ = LI * LJ
*
*                   CHECK LIMIT OF PHI AND CLAMDA
*                   MINIMUM AND MAXIMUM VALUE OF FIELD
*
               IF(PHI.GE.0.0 .AND.PHI.LE.180.0.AND.
     A   CLAMDA.GE. -180.0 .AND. CLAMDA .LT. 360.0
     B.AND.LI.GE.4  .AND. LJ.GE.4 .AND. LILJ.LT.70000) GO TO 12
               WRITE(6,600) CLAMDA, PHI ,LI,LJ
      CALL ABORT
                    STOP


   12          IF( PHI.GT.EPSPHI.AND.PHI.LT. (180. -EPSPHI))  GO TO 13


*                   SPECIAL INTER NEAR POLE
* ------------------------------------------------------------------------------
      CALL NPSPIN (SPEED,PSI,U,V,LI,LJ,KASE,PHI,ARBPSI,DELLAM)
       RETURN


* ------------------------------------------------------------------------------
*                   DEFINE PARAMETER AND COMPUTE COORDINATE
* ------------------------------------------------------------------------------


13    CONTINUE
      DELLAM = 360./LI
      MU = COS(PHI*DGTORD)
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
* BETWEEN POLE AND FIRST ROW
*
    1 JJ = 0
      A = 2.0 -ROOTS(2)
      B = 2.0 -ROOTS(1)
      C = ROOTS(1)
      D = ROOTS(2)
      GO TO 30
*
* BETWEEN FIRST AND SECOND ROW
*
    2 JJ = 1
      A = 2.0 -ROOTS(1)
      B = ROOTS(1)
      C = ROOTS(2)
      D = ROOTS(3)
      GO TO 30
*
* FAR ENOUGH FROM POLES
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
* BETWEEN LAST TWO ROWS
*
    4 JJ = NPOLY-1
      A = ROOTS(NPOLY-2)
      B = ROOTS(NPOLY-1)
      C = ROOTS(NPOLY)
      D = -2.0-ROOTS(NPOLY)
      GO TO 30
*
* BETWEEN LAST ROW AND SOUTH POLE
*
    5 JJ = NPOLY
      A = ROOTS(NPOLY-1)
      B = ROOTS(NPOLY)
      C = -2.0 - ROOTS(NPOLY)
      D = -2.0 - ROOTS(NPOLY-1)
      GO TO 30
*
* POSITION IS NOW KNOWN, LET'S DO THE WORK
*
   30 CONTINUE
      RDLAM = FLOAT(LI)/360.0
      XX = 1.0 + CLAMDA*RDLAM
      II = XX
      DX = XX-II
      CALL IDXES(IDX,JDX,II,JJ,LI,LJ,KASE)


*
* INTERPOLATION DONE HERE, FIRST 4 IN E-W DIRECTION, THEN 1 IN N-S
* DIRECTION. THE WIND COMPONENTS ARE CONVERTED TO VALUES ON A P-S
* PROJECTION ALIGNED WITH GREENWICH.
*
      DO 50 JJ = 1 , 4
        DO 40 II = 1 , 4
          III = IDX(II,JJ)
          JJJ = IABS(JDX(JJ))
          UU  = U(III,JJJ)
          VV  = V(III,JJJ)
          IF(JDX(JJ).LT.0) VV = VV * SIGN0
          PLAMD = (III-1)*DELLAM*DGTORD
          SS = - SIN(PLAMD)
          CC =   COS(PLAMD)
          SUPS(II) = SS*UU - CC*VV
          SVPS(II) = CC*UU + SS*VV
   40     CONTINUE
*
        FUPS(JJ) = CUBIC (SUPS(1),SUPS(2),SUPS(3),SUPS(4),DX)
       FVPS(JJ) = CUBIC (SVPS(1),SVPS(2),SVPS(3),SVPS(4),DX)
*
   50   CONTINUE
*
      UFIN = POLY(FUPS(1),FUPS(2),FUPS(3),FUPS(4),MU,A,B,C,D)
      VFIN = POLY(FVPS(1),FVPS(2),FVPS(3),FVPS(4),MU,A,B,C,D)
*
      SPEED = SQRT (UFIN**2 + VFIN**2)
                        PSI = ARBPSI
      IF (SPEED.NE.0.0) PSI = RDTODG*ATAN2(VFIN,UFIN)
      IF (PSI.LT.0.0) PSI = PSI + 360.0
      RETURN
  600 FORMAT(//41H LAMBDA PHI LI LJ SBRT VGLOGG WRONG VALUE,2E12.4,2I7)
      END
