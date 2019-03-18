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
***S/R D1 - ROUTINE USED IN SPLINE INTERPOLATION.
*
      FUNCTION D1 (Y,H,N)
*
* ENTRY D1
* ENTRY DN
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(D1)
*         - D1 RETURNS A ONE-SIDED ESTIMATE OF THE DERIVATIVE AT THE
*           FIRST POINT, 4-TH ORDER ACCURATE.
*         - DN RETURNS A ONE-SIDED ESTIMATE OF THE DERIVATIVE AT THE
*           LAST POINT, 4-TH ORDER ACCURATE.
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=D1
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - R = D1(Y,H,N)
*           R = DN(Y,H,N)
*
*ARGUMENTS
*         - Y - CONTAINS N POINTS EQUALLY SPACED.
*         - H - SEPARATION.
*         - N - LENGTH OF VECTOR Y.
*
*NOTES    - A SET OF PROGRAMS IS AVAILABLE TO DO 1 AND 2 DIMENSION
*           SPLINE INTERPOLATION. (A 1-DIMENSIONAL SPLINE IS A CURVE
*           CONSTRAINED TO PASS THROUGH CERTAIN GIVEN POINTS, VARYING
*           CUBICALLY BETWEEN THE GIVEN POINTS. THE DIFFERENT CUBIC
*           SEGMENTS JOIN UP SMOOTHLY, WITH CONTINUOUS 1ST AND 2ND DERIVATIVES.
*           THE 2-DIMENSIONAL SPLINE VARIES BICUBICALLY IN GRID SQUARES,
*           PASSES THROUGH ALL THE ORIGINAL POINTS, AND HAS CONTINUOUS
*           FIRST AND SECOND DERIVATIVES, (AS FAR AS F(XXYY)).)
*         - THE SPACING OF THE POINTS IS ARBITRARY FOR THIS SET OF PROGRAMS.
*           (THE 2-DIMENSIONAL GRID IS ASSUMED TO HAVE ORTHOGONAL
*           COORDINATE AXES).
*         - N.B.  D1 AND DN ASSUME CONSTANT SPACING.
*
*---------------------------------------------------------------------------
*
      REAL Y(N)
*
* CONSTANTS
*
      DATA  C1 /-21.0/  ,
     1      C2 /+13.0/  ,
     2      C3 /+17.0/  ,
     3      C4 /- 9.0/  ,
     4     R20 /0.05 /
*
*---------------------------------------------------------------------
*
* D1
*
      I = 1
      IP = 1
      SIGN = 1.0
      GO TO 10
*
* DN
*
      ENTRY DN(Y,H,N)
*
      I = -1
      IP = N
      SIGN = -1.0
*
   10 D1 = SIGN * ( C1*Y(IP)+C2*Y(IP+I)+C3*Y(IP+2*I)+C4*Y(IP+3*I))
     1                  * R20/H
*
*-------------------------------------------------------------------
*
      RETURN
      END
