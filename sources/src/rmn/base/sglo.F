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
***FONCTION SGLO - RETURNS THE INTERPOLATED VALUE FROM THE
*                  FIELD F, AT THE POINT WITH CO-LATITUDE
*                  PHI, LONGITUDE LAMBDA.


      FUNCTION SGLO (LAMBDA,PHI,F,LI,LJ)
      REAL F(LI,LJ)
*
*AUTEUR   - D. ROBERTSON
*
*REVISION 001  C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGAGE  - FORTRAN
*
*OBJET(SGLO)
*         - RETURNS THE INTERPOLATED VALUE FROM THE FIELD F,
*           AT THE POINT WITH CO-LATITUDE PHI, LONGITUDE LAMBDA.
*           RANGES PHI 0 TO PI, 0 AT NORTH POLE, LAMBDA, EITHER
*           0 TO 2 PI OR -PI TO +PI, 0 IS GREENWICH, +VE TOWARDS
*           THE EAST.
*
*ALGORITHME
*         - A BI-CUBIC INTERPOLATION IS USED, WHICH TAKES INTO
*           ACCOUNT THE PERIODICITY IN  BOTH COORDINATES.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=SGLO
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - R=SGLO(LAMBDA,PHI,F,LI,LJ)
*
*ARGUMENTS
*   OUT   - SGLO   - INTERPOLATED VALUE
*   IN    - LAMBDA - LONGITUDE
*         - PHI    - CO-LATITUDE (90. - LAT)
*         - F      - A UNIFORM LON-LAT. GRID (STAGGERED)
*         - LI     - POINTS PER CIRCLE. MUST BE EVEN. THE FIRST
*                    CIRCLE IS DISPLACED FROM THE NORTH POLE BY
*                    OFFSET = 90.0 / LJ DEG.
*         - LJ     - NUMBER OF LATITUDE CIRCLES
*
*MODULES  - ABORT,IDXES
*
*MESSAGES - "LAMBDA PHI, LI LJ IN FUNCTION SGLO (SSNH,SANH)
*            WRONG VALUE, ..-..-..-.."
*
*----------------------------------------------------------------------------
*
      REAL LAMBDA , LONG
*
* THE VALUES IN JDX(4) LOCATE THE 4 LATITUDE CIRCLES ON WHICH THE
* POINTS NEEDED LIE. THE CORRESPONDING COLUMNOF IDX, LOCATES THE
* POINTS ON THE LATITUDE CIRCLE. ALTERNATELY A ROW OF IDX, TOGETHER WITH
* THE CORRESPONDING VALUES IN JDX, GIVES 4 POINTS ON LINE OF CONSTANT
* LONGITUDE.
*
#include "ordr.cdk"
      DIMENSION IDX(4,4) , JDX(4)
*
* USED TO HOLD THE 4 POINTS ALONG A LINE OF CONSTANT LONGITUDE, USED
* TO INTERPOLATE TO THE FINAL LATITUDE CIRCLE.
*
      DIMENSION  P(4)
*
* USED TO HOLD THE 4 POINTS ON THE FINAL LATITUDE CIRCLE, USED TO
* INTERPOLATE TO THE FINAL VALUE.
*
      DIMENSION Z(4)
*
* SPECIAL TREATMENT FOR POINTS AT OR VERY NEAR POLES
* EPSPHI-AN ARBITRARY LIMIT. POINTS WITHIN A CONE OF HALF-ANGLE
* EPSPHI AROUND THE POLE ARE CONSIDERED TO BE AT THE POLE. THIS
* SIMPLE SCHEME NEEDS A BIT OF EXAMINATION.
#include "kases.cdk"
*
      DATA  EPSPHI /1.0E-8/
*      DATA ORDINT /3/
#include "cubic.cdk"
*
      ORDINT = 3
      KASE = GLOBAL
      SIGN =+1.0
      GO TO 5
*
*     SCALAR SYMMETRIC NORTHERN HEMISPHERE
*
      ENTRY SSNH(LAMBDA,PHI,F,LI,LJ)
      KASE = NORTH
      SIGN =+1.0
      GO TO 5
*
*     SCALAR ANTISYMMETRIC NORTHERN HEMISPHERE
*
      ENTRY SANH(LAMBDA,PHI,F,LI,LJ)
      KASE = NORTH
      SIGN =-1.0
    5 CONTINUE
*    -------------------------------------------------------------------
*
*         CHECK LIMIT OF PHI LAMDA,LI,LJ
       LILJ = LI * LJ
                    IF(PHI.GE.0.0.AND.PHI.LE.180.0.AND.
     A                 LAMBDA.GE.-180.0.AND.LAMBDA.LT. 360.0
     B   .AND.LI.GE.4  .AND. LJ.GE.4  .AND.LILJ.LT.70000) GO TO 60
                    WRITE(6,620) LAMBDA,PHI,LI,LJ
       CALL ABORT
                    STOP
*    -------------------------------------------------------------------
*
   60  IF(PHI.GT.EPSPHI.AND. PHI.LT.(180.0- EPSPHI)) GO TO 70
*----------------------------------------------------------------------
*         CASE  FOR POLE
*
*
      J1 = 1
      J2 = 2
      IF(PHI.GT.90.0 .AND. KASE.EQ.GLOBAL)
     $         J1 = LJ
      IF(PHI.GT.90.0 .AND. KASE.EQ.GLOBAL)
     $        J2 = LJ - 1
      SUM1 = 0.
      SUM2 = 0.
      DO 65 I = 1,LI
      SUM1 = SUM1 + F(I,J1)
65    SUM2 = SUM2 + F(I,J2)
      SGLO = (SUM1*1.125 - SUM2*.125)/LI
       IF(KASE.EQ.NORTH.AND.PHI.GT.90.0)  SGLO = SIGN * SGLO
      RETURN
*
*----------------------------------------------------------------------
*
   70  OFFSET  =  90.0 / LJ
      IF (KASE.NE.GLOBAL) OFFSET = 0.5*OFFSET
      DELPHI = 2.0*OFFSET
      DELLAM = 360.0/LI
*
* LONG RUNS FROM 0 TO 2PI EAST.
*
      LONG = LAMBDA
      IF (LAMBDA .LT.0.0)
     $     LONG = 360.0 + LAMBDA
*
      XI = 1.0 + LONG/DELLAM
      YJ = 1.0+(PHI-OFFSET)/DELPHI
*
* POINT HAS COORDINATES II.DX, JJ.DY
*
      II = XI
      JJ = YJ
*
      DX = XI-II
      DY = YJ-JJ
*
*
* FIND THE INDICES OF THE 16 POINTS TO BE USED.
*
      CALL IDXES (IDX,JDX,II,JJ,LI,LJ,KASE)
      IF(ORDINT.EQ.0) THEN
         II = NINT(DX+2.0)
         JJ = NINT(DY+2.0)
         SGLO = F(IDX(II,JJ),ABS(JDX(JJ)))
         IF(JDX(JJ).LT.0) SGLO = SGLO*SIGN
      ELSE IF(ORDINT.EQ.1) THEN
         DO 90 JJ = 2 , 3
            DO 80 II = 2 , 3
80             P(II) = F(IDX(II,JJ),ABS(JDX(JJ)))
            IF(JDX(JJ).LT.0) THEN
               Z(JJ) = SIGN * ZLIN(P(2),P(3),DX)
            ELSE
               Z(JJ) = ZLIN(P(2),P(3),DX)
            ENDIF
90       CONTINUE
         SGLO = ZLIN(Z(2),Z(3),DY)
      ELSE
*
* INTERPOLATE 4 POINTS ON THE FINAL LATITUDE CIRCLE, STORE IN Z.
* APPLY SYMETRY IF OVER THE EQUATOR
*
      DO 200 II = 1 , 4
        DO 100 JJ = 1 , 4
          P(JJ) = F(IDX(II,JJ),IABS(JDX(JJ)))
          IF (JDX(JJ).LT.0)
     $       P(JJ) = SIGN*P(JJ)
  100     CONTINUE
       Z(II)= CUBIC(P(1),P(2),P(3),P(4),DY)
  200    CONTINUE
*
* FINAL INTERPOLATION.
*
       SGLO=CUBIC(Z(1),Z(2),Z(3),Z(4),DX)
      ENDIF
      RETURN
*    -------------------------------------------------------------------
*
  620 FORMAT(/66H LAMDA PHI, LI LJ IN FUNCTION SGLO(SSNH,SANH)
     1WRONG VALUE, 2E12.4,2I7)
*    -------------------------------------------------------------------
      END
