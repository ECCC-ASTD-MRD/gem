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
***S/P VGLO - RETURN 2  SCALAR QUANTITIES FROM U AND V FIELD
*             AT POINT CLAMDA (LONGITUDE) AND PHI(LATITUDE)
*
      SUBROUTINE VGLO(SPEED ,  PSI , CLAMDA,  PHI , U, V, LI, LJ )
      REAL U(LI,LJ), V(LI,LJ)
*
*AUTEUR   - P.SARRAZIN DEC 78 DRPN DORVAL QUEBEC CANADA
*
*REVISION 001  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE - FORTRAN
*
*OBJECT(VGLO)
*         - RETURN 2 SCALARS QUANTITIES USING U-V VECTOR FIELD
*         - ON LAT LON COMPUTE SPEED AND PSI WITH A CUBIC
*         - INTERPOLATION ON EAST TO WEST  AND  NORTH TO SOUTH
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=VGLO
*         - OBJECT  RMNLIB,ID=RMNP
*
*APPEL    - CALL VGLO(SPEED,PSI,CLAMDA,PHI,U,V,LI,LJ)
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
*
*MODULES  - NPSPIN,IDXES,ABORT
*
*MESSAGES - "LAMBDA PHI LI LJ IN SBRT VGLO(VWNH) WRONG VALUE, ..-..-..-.."
*
*NOTES    - IF SPEED = 0.0, PSI = ARBPSI
*         - SPEED AND PSI ARE GIVEN RELATIVE TO A P.S.
*           GRID WITH DGRW=0.0
*
*-------------------------------------------------------------------------------
*
*
* - VALUE IN JDX(4) LOCATE THE 4 LAT CIRCLES ON WHICH THE
* - POINT NEEDED LIE
* - IDX LOCATE  THE POINT ON LATITUDE CIRCLE
      DIMENSION  IDX(4,4) , JDX(4)
*
* SUPS,SVPS HOLD 4 VALUES OF U,V ALONG A MERIDION OF LONGITUDE,
* WHICH ARE USED TO INTERPOLATE A VALUE AT THE FINAL LATITUDE.
* THE 4 VALUES AT THE FINAL LATITUDE ARE HELD IN FUPS,FVPS.
*
      DIMENSION SUPS(4) , SVPS(4) , FUPS(4) , FVPS(4)
* POLAR-STEREOGRAPHIC COMPONENTS ARE USED FOR THE INTERPOLATION.


*             -EPSPHI AN ARBITRARY LIMIT


#include "kases.cdk"


      DATA EPSPHI / 1.0E-8 /
#include "pi.cdk"
      DATA ARBPSI /99999.0/


#include "cubic.cdk"


       KASE = GLOBAL
      SIGN0 = 1.0
      NLATS = LJ
      GOTO 1
*
* HEMISPHERIC, WALL AT THE EQUATOR
*
      ENTRY VWNH(SPEED,PSI,CLAMDA,PHI,U,V,LI,LJ)
      KASE = NORTH
      SIGN0 = -1.0
      NLATS = LJ*2
1      LILJ = LI * LJ


*                   CHECK LIMIT OF PHI AND CLAMDA
*                   MINIMUM AND MAXIMUM VALUE OF FIELD
               IF(PHI.GE.0.0 .AND.PHI.LE.180.0.AND.
     A   CLAMDA.GE. -180.0 .AND. CLAMDA .LT. 360.0
     B.AND.LI.GE.4  .AND. LJ.GE.4 .AND. LILJ.LT.70000) GO TO 10
               WRITE(6,600) CLAMDA, PHI ,LI,LJ
      CALL ABORT
                    STOP


   10          IF( PHI.GT.EPSPHI.AND.PHI.LT. (180. -EPSPHI))  GO TO 30


*                   SPECIAL INTER NEAR POLE
* ------------------------------------------------------------------------------
      CALL NPSPIN (SPEED,PSI,U,V,LI,LJ,KASE,PHI,ARBPSI,DELLAM)
      RETURN


* ------------------------------------------------------------------------------
*                   DEFINE PARAMETER AND COMPUTE COORDINATE
* ------------------------------------------------------------------------------


   30  OFFSET  =  90.0 /NLATS
       DELPHI  =  2.0 * OFFSET
       DELLAM  =  360.0 / LI


*                   CLONG 0 TO 2 PI-DELLAM EAST +


      CLONG = CLAMDA
      IF( CLAMDA .LT.0.0 )  CLONG = 360.0 +CLAMDA


      XI  =  1.0 + CLONG / DELLAM
      YJ  =  1.0 + (PHI - OFFSET) / DELPHI


      II = XI
      JJ = YJ


      DX = XI - II
      DY = YJ - JJ


*                   FIND INDICES FOR 16 PTS


      CALL IDXES (IDX, JDX, II, JJ, LI,LJ, KASE)
*
* INTERPOLATION DONE HERE, FIRST 4 IN N-S DIRECTION, THEN 1 IN E-W
* DIRECTION. THE WIND COMPONENTS ARE CONVERTED TO VALUES ON A P-S
* PROJECTION ALIGNED WITH GREENWICH.
*
      DO 50 II = 1 , 4
        DO 40 JJ = 1 , 4
          III = IDX(II,JJ)
          JJJ = IABS(JDX(JJ))
          UU  = U(III,JJJ)
          VV  = V(III,JJJ)
          IF(JDX(JJ).LT.0) VV = VV * SIGN0
          PLAMD = (III-1)*DELLAM*DGTORD
          S = - SIN(PLAMD)
          C =   COS(PLAMD)
          SUPS(JJ) = S*UU - C*VV
          SVPS(JJ) = C*UU + S*VV
   40     CONTINUE
*
        FUPS(II) = CUBIC (SUPS(1),SUPS(2),SUPS(3),SUPS(4),DY)
       FVPS(II) = CUBIC (SVPS(1),SVPS(2),SVPS(3),SVPS(4),DY)
*
   50   CONTINUE
*
      UFIN = CUBIC (FUPS(1),FUPS(2),FUPS(3),FUPS(4),DX)
      VFIN = CUBIC (FVPS(1),FVPS(2),FVPS(3),FVPS(4),DX)
*
      SPEED = SQRT (UFIN**2 + VFIN**2)
                        PSI = ARBPSI
      IF (SPEED.NE.0.0) PSI = RDTODG*ATAN2(VFIN,UFIN)
      IF (PSI.LT.0.0) PSI = PSI + 360.0
      RETURN
  600 FORMAT( //39H LAMBDA PHI LI LJ SBRT VGLO WRONG VALUE, 2E12.4,2I7)
      END
