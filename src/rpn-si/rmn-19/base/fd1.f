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
***S/R FD1 - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE FD1(GD1,F,H)
      REAL F(*),H(*)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(FD1)
*         - THIS ROUTINE RETURNS ONE-SIDED ESTIMATE OF THE DERIVATIVE
*           AT THE FIRST POINT, 4TH ORDER ACCURATE.
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=FD1
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL FD1(GD1,F,H)
*
*ARGUMENTS
*         - GD1 - RESULT.
*         - F   - VECTOR OF VALUES OF THE FUNCTION.
*         - H   - VECTOR OF COORDINATE DIFFERENCES.
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
*
*------------------------------------------------------------------
*
      H1=H(1)
      H2=H(2)
      H3=H(3)
      X1=H1+H2
      X2=H2+H3
      X3=X1+H3
      GD1=-(1./H1+1./X1+1./X3)*F(1)+X1*X3*F(2)/(H1*H2*X2)
      GD1=GD1-H1*X3*F(3)/(X1*H2*H3)+H1*X1*F(4)/(X3*X2*H3)
*
*---------------------------------------------------------------------
*
      RETURN
      END
