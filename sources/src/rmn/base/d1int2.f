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
***S/R D1INT2 - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE D1INT2 (FI,IFI,JFI,F,IF,JF,FX,FY,FXY,XI,YI,X,Y,HX,HY,
     1                   ZA,ZB,ZC,ZD)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(D1INT2)
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=D1INT2
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL D1INT2(FI,IFI,JFI,F,IF,JF,FX,FY,FXY,XI,YI,X,Y,
*                       HX,HY,ZA,ZB,ZC,ZD)
*
*ARGUMENTS
*   OUT   - FI   - INTERPOLATED FIELD.
*   IN    - IFI  - X-DIMENSION OF FIELD FI.
*         - JFI  - Y-DIMENSION OF FIELD FI.
*         - F    - ORIGINAL FIELD.
*         - IF   - X-DIMENSION OF FIELD F.
*         - JF   - Y-DIMENSION OF FILED F.
*         - FX   - ARRAY OF SIZE (IF,JF) THAT CONTAINS COMPUTED
*                  PARTIAL DERIVATIVE OF FIELD F WITH RESPECT TO X.
*         - FY   - PARTIAL DERIVATIVE WITH RESPECT TO Y.
*         - FXY  - PARTIAL SECOND DERIVATIVE WITH RESPECT TO X AND Y.
*         - XI   - X LOCATION OF THE VALUES OF THE INTERPOLATED FIELD FI.
*         - YI   - Y LOCATION OF THE VALUES OF THE INTERPOLATED FIELD FI.
*         - X    - X LOCATION OF THE VALUES OF THE ORIGINAL FIELD F.
*         - Y    - Y LOCATION OF THE VALUES OF THE ORIGINAL FIELD F.
*         - HX   - GRID-LENGTHS DEFINED BY HX(I)=X(I+1)-X(I),I=1,IF-1
*         - HY   - GRID-LENGTHS DEFINED BY HY(I)=Y(I+1)-Y(I),I=1,JF-1
*         - ZA   - WORKING VECTOR OF LENGTH (JFI).
*         - ZB   - WORKING VECTOR OF LENGTH (JFI).
*         - ZC   - WORKING VECTOR OF LENGTH (JFI).
*         - ZD   - WORKING VECTOR OF LENGTH (JFI).
*
*NOTES    - A SET OF PROGRAMS IS AVAILABLE TO DO 1 AND 2 DIMENSION
*           SPLINE INTERPOLATION. (A 1-DIMENSIONAL SPLINE IS A CURVE
*           CONSTRAINED TO PASS THROUGH CERTAIN GIVEN POINTS, VARYING
*           CUBICALLY BETWEEN THE GIVEN POINTS. THE DIFFERENT CUBIC
*           SEGMENTS JOIN UP SMOOTHLY, WITH CONTINUOUS 1ST AND 2ND DERIVATIVES.
*           THE 2-DIMENSIONAL SPLINE VARIES BICUBICALLY IN GRID SQUARES,
*           PASSES THROUGH ALL THE ORIGINAL POINTS, AND HAS CONTINUOUS
*           FIRST AND SECOND DERIVATIVES, (AS FAR AS F(XXYY)).)
*         - THE SPACING OF THE POINTS IS ARBITRARY. (THE 2-DIMENSIONAL
*           GRID IS ASSUMED TO HAVE ORTHOGONAL COORDINATE AXES).
*
*---------------------------------------------------------------------------
*
      REAL FI     (IFI,JFI)
      REAL F      (IF ,JF )
      REAL FX     (IF ,JF )
      REAL FY     (IF ,JF )
      REAL FXY    (IF ,JF )
      REAL XI     (IFI)
      REAL YI     (JFI)
      REAL X      (IF )
      REAL Y      (JF )
      REAL HX     (IF )
      REAL HY     (JF )
      REAL ZA     (JFI)
      REAL ZB     (JFI)
      REAL ZC     (JFI)
      REAL ZD     (JFI)
*
*----------------------------------------------------------------------
*
      LL=2
      DO 15 J=1,JFI
      DO 22 L=LL,JF
      IF(YI(J).LE.Y(L)) GO TO 23
   22 CONTINUE
*
   23 L1=L-1
      LL=L
      WN=YI(J)-Y(L1)
      WE=WN/HY(L1)
      WE1=1.-WE
      WE2=WE1*WE1
      WW=2.*WE
      ZA(J)=WE2*WN
      ZB(J)=WE1*WN*WE
      ZC(J)=WE2*(1.+WW)
      ZD(J)=WE*WE*(3.-WW)
   15 CONTINUE
*
      KK=2
      DO 11 I=1,IFI
      DO 12 K=KK,IF
      IF(XI(I).LE.X(K)) GO TO 13
   12 CONTINUE
*
   13 K1=K-1
      KK=K
      WM=XI(I)-X(K1)
      WZ=WM/HX(K1)
      WZ1=1.-WZ
      WZ2=WZ1*WZ1
      ZZ=2.*WZ
      ZE=WZ2*WM
      ZF=WZ1*WM*WZ
      ZG=WZ2*(1.+ZZ)
      ZL=WZ*WZ*(3.-ZZ)
      LL=2
      DO 115 J=1,JFI
      DO 122 L=LL,JF
      IF(YI(J).LE.Y(L)) GO TO 123
  122 CONTINUE
*
  123 L1=L-1
      LL=L
      Z=ZA(J)*(ZE*FXY(K1,L1)-ZF*FXY(K,L1)+ZG*FY(K1,L1)+ZL*FY(K,L1))
      Z=Z-ZB(J)*(ZE*FXY(K1,L)-ZF*FXY(K,L)+ZG*FY(K1,L)+ZL*FY(K,L))
      Z=Z+ZC(J)*(ZE*FX(K1,L1)-ZF*FX(K,L1)+ZG*F(K1,L1)+ZL*F(K,L1))
      Z=Z+ZD(J)*(ZE*FX(K1,L)-ZF*FX(K,L)+ZG*F(K1,L)+ZL*F(K,L))
      FI(I,J)=Z
  115 CONTINUE
*
   11 CONTINUE
*
*--------------------------------------------------------------------
*
      RETURN
      END
