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
***S/R DEFVEC - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE DEFVEC (A,L,DELTAA,A1)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(DEFVEC)
*         - THIS ROUTINE IS USED TO DEFINE A VECTOR OF COORDINATES
*           WITH CONSTANT SPACING.
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=DEFVEC
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL DEFVEC(A,L,DELTAA,A1)
*
*ARGUMENTS
*         - A      - VECTOR OF COORDINATES.
*         - L      - LENGTH OF VECTOR A.
*         - DELTAA - SPACING CONSTANT.
*         - A1     - 1ST VALUE OF VECTOR OF COORDINATES.
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
      REAL A(L)
*
*-----------------------------------------------------------------------------
*
      AX = A1
*
      DO 50 I=1,L
      A(I) = AX
      AX = AX + DELTAA
   50 CONTINUE
*
*-------------------------------------------------------------------------------
*
      RETURN
      END
