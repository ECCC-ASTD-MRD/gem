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
***S/R D1INTR - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE D1INTR(FI,IFI,JFI,F,IF,JF,FX,FY,FXY,XI,YI,
     $                  X,Y,HX,HY)
*
*AUTHOR   - D. ROBERTSON
*
*
*LANGUAGE - fortran
*
*OBJECT(D1INTR)
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*           THIS VERSION EVALUATES THE SPLINE AT A LIST OF (X,Y) PAIRS
*           HELD IN THE ARRAYS (XI(I,J),YI(I,J),I=1,IFI,J=1,JFI).
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=D1INTR
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL D1INTR(FI,IFI,JFI,F,IF,JF,FX,FY,FXY,XI,YI,X,Y,
*                       HX,HY)
*           D1INTR IS NOT NORMALLY USER CALLED, BUT IS INVOKED FROM INTRPR.
*
*ARGUMENTS
*   OUT   - FI   - INTERPOLATED FIELD.
*   IN    - IFI  - X-DIMENSION OF FIELD FI.
*         - JFI  - Y-DIMENSION OF FIELD FI.
*         - F    - ORIGINAL FIELD.
*         - IF   - X-DIMENSION OF FIELD F.
*         - JF   - Y-DIMENSION OF FIELD F.
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
      REAL XI     (IFI,JFI)
      REAL YI     (IFI,JFI)
      REAL X      (IF )
      REAL Y      (JF )
      REAL HX     (IF )
      REAL HY     (JF )
*
*----------------------------------------------------------------------
      DO 750 J = 1 , JFI
      DO 725 I = 1 , IFI
*
* FIND K,L S.T. X(K-1).LT.XPT.LE.X(K), Y(L-1).LT.YPT.LE.Y(L)
*
        XPT = XI(I,J)
        YPT = YI(I,J)
        DO 100 K=2,IF
          IF (XPT.LE.X(K))
     $        GO TO 123
  100     CONTINUE
        K = IF
  123 CONTINUE
        DO 200 L=2,JF
           IF(YPT.LE.Y(L))
     $         GO TO 225
  200      CONTINUE
        L = JF
  225   CONTINUE
*
* EVALUATE SPLINE
*
        WN = YPT-Y(L-1)
        WE = WN/HY(L-1)
        WE1 = 1.0-WE
        WE2 = WE1**2
        WW = 2.0*WE
        ZA  = WE2*WN
        ZB  = WE1*WN*WE
        ZC  = WE2*(1.0+WW)
        ZD  = (3.0-WW)*WE**2
        WM  = XPT-X(K-1)
        WZ  = WM/HX(K-1)
        WZ1 = 1.0-WZ
        WZ2 = WZ1**2
        ZZ  = 2.0*WZ
        ZE  = WZ2*WM
        ZF  = WZ1*WM*WZ
        ZG  = WZ2*(1.0+ZZ)
        ZL  = (3.0-ZZ)*WZ**2
        Z   = ZA*(ZE*FXY(K-1,L-1)-ZF*FXY(K,L-1)+ZG*FY(K-1,L-1)+ZL*
     $          FY(K,L-1))
        Z=Z-ZB*(ZE*FXY(K-1,L)-ZF*FXY(K,L)+ZG*FY(K-1,L)+ZL*FY(K,L))
        Z=Z+ZC*(ZE*FX(K-1,L-1)-ZF*FX(K,L-1)+ZG*F(K-1,L-1)+ZL*F(K,L-1))
        FI(I,J)=Z+ZD*(ZE*FX(K-1,L)-ZF*FX(K,L)+ZG*F(K-1,L)+ZL*F(K,L))
C
  725   CONTINUE
  750   CONTINUE
      RETURN
      END
