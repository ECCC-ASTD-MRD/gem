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
***S/R D1INT1 - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE D1INT1(FI,F,XI,X,FX,H,M,MI)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(D1INT1)
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           HOW TO USE THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=D1INT1
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL D1INT1(FI,F,XI,X,FX,H,M,MI)
*
*ARGUMENTS
*   OUT   - FI - VECTOR OF INTERPOLATED VALUES.
*   IN    - F  - VECTOR OF VALUES OF THE FUNCTION.
*         - XI - LOCATION OF THE VALUES OF THE INTERPOLATED VECTOR.
*         - X  - LOCATION OF THE VALUES OF INPUT VECTOR F.
*         - FX - VECTOR CONTAINING THE SPLINE DERIVATIVES OF THE
*                FUNCTION AT THE ORIGINAL POINT (X).
*         - H  - VECTOR WHICH MUST BE INITIALIZED WITH THE
*                CO-ORDINATE DIFFERENCES, I.E. H(I)=X(I+1)-X(I) FOR
*                 I=1 TO (M-1). (THE M TH ELEMENT OF H IS NOT USED.)
*         - M  - LENGTH OF VECTORS F, X AND H.
*         - MI - LENGTH OF VECTORS FI AND XI.
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
      REAL     FI(MI)
      REAL      F(M)
      REAL     XI(MI)
      REAL      X(M)
      REAL     FX(M)
      REAL      H(M)
*
*------------------------------------------------------------------
*
      JJ=2
      DO 11 I=1,MI
      DO 12 J=JJ,M
      IF(XI(I).LE.X(J)) GO TO 13
   12 CONTINUE
*
   13 J1=J-1
      JJ=J
      V1=X(J)-XI(I)
      W1=XI(I)-X(J1)
      V2=V1*V1
      W2=W1*W1
      H1=H(J1)
      H2=1./(H1*H1)
      H3=H2/H1
      Z=(FX(J1)*V2*W1-FX(J)*W2*V1)*H2
      FI(I)=Z+(F(J1)*V2*(2.*W1+H1)+F(J)*W2*(2.*V1+H1))*H3
   11 CONTINUE
*
*--------------------------------------------------------------------
*
      RETURN
      END
