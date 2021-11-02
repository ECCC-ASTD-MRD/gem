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
***S/R SPD - ROUTINE USED IN SPLINE INTERPOLATION.
*
      SUBROUTINE SPD (P,S,M,H,CMU1,C1,CLMDAM,CM,A,C,D)
*
*AUTHOR   - D. ROBERTSON
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(SPD)
*         - THIS ROUTINE IS USED TO COMPUTE SPLINE DERIVATIVES.
*         - THIS ROUTINE IS PART OF THE SPLINE INTERPOLATION PACKAGE.
*           SEE THE WRITE-UPS ON INTRP AND INT1D1 FOR INFORMATION ABOUT
*           USING THE PACKAGE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=SPD
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL SPD(P,S,M,H,CMU1,C1,CLMDAM,CM,A,C,D)
*
*ARGUMENTS
*         - P      - SOLUTION VECTOR.
*         - S      - VECTOR OF VALUES OF THE FUNCTION.
*         - M      - LENGTH OF VECTORS P, S, A, C, D AND H.
*         - H      - VECTOR WHICH MUST BE INITIALIZED WITH THE
*                    COORDINATE DIFFERENCES.
*         - CMU1   - CONSTANT USED TO SPECIFY THE END CONDITIONS
*                    FOR THE SPLINE.
*         - C1     - CONSTANT USED TO SPECIFY THE END CONDITIONS
*                    FOR THE SPLINE.
*         - CLMDAM - CONSTANT USED TO SPECIFY THE END CONDITIONS
*                    FOR THE SPLINE.
*         - CM     - CONSTANT USED TO SPECIFY THE END CONDITIONS
*                    FOR THE SPLINE.
*         - A      - WORKING VECTOR OF LENGTH M.
*         - C      - WORKING VECTOR OF LENGTH M.
*         - D      - WORKING VECTOR OF LENGTH M.
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
*MODULES  - ROSSR3
*
*---------------------------------------------------------------------------
*
      REAL      P(M)
      REAL      S(M)
      REAL      A(M)
      REAL      C(M)
      REAL      D(M)
      REAL      H(M)
*
*---------------------------------------------------------------------
*
      M1=M-1
      DO 1 I=2,M1
      I1=I-1
      X1=0.5/(H(I)+H(I1))
      X2=H(I)/H(I1)
      A(I)=H(I)*X1
      C(I)=H(I1)*X1
    1 D(I)=3.*((S(I)-S(I1))*X2+(S(I+1)-S(I))/X2)*X1
*
      C(1)=CMU1
      A(M)=CLMDAM
      D(1)=C1
      D(M)=CM
      CALL ROSSR3(P,A,D,C,D,M)
*
*-------------------------------------------------------------------
*
      RETURN
      END
