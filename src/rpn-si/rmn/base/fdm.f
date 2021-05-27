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
***S/R FDM - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE FDM(GDM,F,H,M)
      REAL F(*),H(*)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(FDM)
*         - THIS ROUTINE RETURNS AN ESTIMATE OF THE DERIVATIVE
*           AT THE LAST POINT, 4-TH ORDER ACCURATE.
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=FDM
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL FDM(GDM,F,H,M)
*
*ARGUMENTS
*         - GDM - RESULT.
*         - F   - VECTOR OF VALUES OF THE FUNCTION.
*         - H   - VECTOR OF COORDINATE DIFFERENCES.
*         - M   - DIMENSION OF F AND H.
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
*----------------------------------------------------------------
*
      H1 = H(M-1)
      H2 = H(M-2)
      H3 = H(M-3)
      X1 = H1+H2
      X2 = H2+H3
      X3 = X1+H3
      GDM=-X1*H1*F(M-3)/(H3*X2*X3)+X3*H1*F(M-2)/(H3*H2*X1)
      GDM= GDM-X3*X1*F(M-1)/(X2*H2*H1)+(1.0/X3+1.0/X1+1.0/H1)*F(M)
*
*------------------------------------------------------------------
*
      RETURN
      END
