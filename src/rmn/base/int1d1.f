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
***S/R INT1D1 - ONE-DIMENSIONAL SPLINE INTERPOLATION.
*
      SUBROUTINE INT1D1(FI,F,XI,X,FX,H,M,MI,CMU1,C1,CLMDAM,CM,A,C,D)
*
*AUTHOR   - D. ROBERTSON - A. STANIFORTH
*
*REVISION 001   C. THIBEAULT - FEB 80 - DOCUMENTATION
*
*LANGUAGE - fortran
*
*OBJECT(INT1D1)
*         - INTERPOLATE BETWEEN A NUMBER OF ARBITRARILY SPACED DATA
*           POINTS USING CUBIC SPLINES.
*
*ALGORITHM
*         - A 1-DIMENSIONAL INTERPOLATING CUBIC SPLINE IS DEFINED
*           AS FOLLOWS...
*           THE INTERPOLATING CURVE PASSES THROUGH ALL THE ORIGINAL
*           DATA POINTS, AND VARIES CUBICALLY BETWEEN POINTS. THE
*           CUBICS BETWEEN THE POINTS JOIN UP SMOOTHLY, IN THE SENSE
*           THAT THE FIRST AND SECOND DERIVATIVES ARE CONTINUOUS.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=INT1D1
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL INT1D1 (FI,F,XI,X,FX,H,M,MI,CMU1,C1,CLMDAM,CM,A,C,D)
*
*ARGUMENTS
*   OUT   - FI     - VECTOR OF INTERPOLATED VALUES.
*   IN    - F      - VECTOR OF VALUES OF THE FUNCTION.
*         - XI     - LOCATION OF THE VALUES OF THE INTERPOLATED VECTOR.
*         - X      - LOCATION OF THE VALUES OF INPUT VECTOR F.
*   OUT   - FX     - VECTOR CONTAINING THE SPLINE DERIVATIVES OF THE
*                    FUNCTION AT THE ORIGINAL POINTS (X).
*   IN    - H      - VECTOR WHICH MUST BE INITIALIZED WITH THE
*                    CO-ORDINATE DIFFERENCES, I.E. H(I)=X(I+1)-X(I) FOR
*                     I=1 TO (M-1). (THE M TH ELEMENT OF H IS NOT USED.)
*         - M      - LENGTH OF VECTORS F, FX, X, H, A, D AND C.
*         - MI     - LENGTH OF VECTORS FI AND XI.
*         - CMU1   - CMU1,C1,CLMDAM,CM ARE CONSTANTS WHICH ARE USED TO
*         - C1       SPECIFY THE END CONDITIONS FOR THE SPLINE. ONE
*         - CLMDAM   SUITABLE CHOICE IS TO SPECIFY THE END POINT SLOPES,
*         - CM       WHICH IS DONE BY SETTING CMU1=CLMDAM=0.0, AND
*                    PUTTING THE DESIRED SLOPES IN C1 AND CM. SUBROUTINES
*                    FD1 AND FDM, AND FUNCTIONS D1 AND DN ARE AVAILABLE
*                    TO CALCULATE END POINT SLOPES. FD1 AND FDM USE THE
*                    ANALYTICALLY DIFFERENTIATED FORM OF THE CUBIC
*                    LAGRANGE INTERPOLATING POLYNOMIAL, AND MAY BE USED
*                    FOR GENERALLY SPACED POINTS. THEY ARE USED AS...
*                          CALL FD1 (GD1,F,H)
*                          CALL FDM (GDM,F,H,M)
*                   THE DERIVATIVES ARE RETURNED IN GD1 OR GDM. (F,H,M
*                    AS BEFORE). THE FUNCTIONS D1 AND DN, (DN IS AN ENTRY
*                    POINT IN D1), ARE SPECIALIZED TO UNIFORM SPACING,
*                    AND ARE USED AS ...
*                          GD1 = D1(F,H,M)
*                          GDM = DN(F,H,M)
*                    WHERE H IS NOW THE CONSTANT VALUE OF X(I+1) - X(I).
*                    (F,M AS BEFORE.)
*         - A      - WORKING VECTOR OF LENGTH M.
*         - C      - WORKING VECTOR OF LENGTH M.
*         - D      - WORKING VECTOR OF LENGTH M.
*
*NOTES    -    I) ALL VARIABLES ARE IN SCM (LEVEL 1).
*         -   II) THE LOCATION AND SPACING OF THE POINTS X(I), XI(I) ARE
*                 COMPLETELY ARBITRARY. THE X AND XI VECTORS MUST BE
*                 GIVEN IN STRICTLY INCREASING ORDER X(I)<X(I+1), (NO XS
*                 OR XIS EQUAL). SUBROUTINE SORT IS AVAILABLE FOR THIS.
*         -  III) ONE CAN USE THE SUBROUTINE GRAPH, TO CHECK THE OUTPUT
*                 FROM THE INTERPOLATION.
*         -   IV) A TWO-DIMENSIONAL CUBIC SPLINE PACKAGE IS ALSO
*                 AVAILABLE - SEE THE WRITE-UP ON INTRP FOR INFORMATION.
*         - *** WARNING *** EXTRAPOLATION IS NOT ALLOWED. THAT IS, ONE
*                           MUST HAVE X(1).LE.XI(J).LE.X(M) FOR J=1,MI.
*         - FOR UNIFORMLY SPACED POINTS THE SUBROUTINE DEFVEC IS
*           AVAILABLE TO INITIATE THE VECTORS X,XI APPROPRIATELY.
*           CALL DEFVEC(X,M,DELTAX,X1) INITIALIZES X(M) SUCH THAT
*           X(1) = X1 AND X(I+1) = X(I)+DELTAX, I=1,...M-1.
*
*MODULES  - SPD    TO GET SPLINE DERIVATIVES.
*         - ROSSR3 SOLVES A TRIDIAGONAL SYSTEM OF LINEAR EQUATIONS.
*                  WITH AN ALGORITHM WHICH IS EFFICIENT AND
*                  COMPUTATIONALLY STABLE WHEN THE MATRIX IS DIAGONALLY
*                  DOMINANT, (AS IT IS IN THIS CASE).
*         - D1INT1 TO CALCULATE THE INTERPOLATED VALUES.
*         - D1, DN, FD1 AND FDM USED AS OUTLINED ABOVE.
*
*--------------------------------------------------------------------------
*
      REAL     FI(MI)
      REAL      F(M)
      REAL     XI(MI)
      REAL      X(M)
      REAL     FX(M)
      REAL      H(M)
      REAL      A(M)
      REAL      D(M)
      REAL      C(M)
*
*--------------------------------------------------------------------------
*
      CALL SPD(FX,F,M,H,CMU1,C1,CLMDAM,CM,A,C,D)
      CALL D1INT1(FI,F,XI,X,FX,H,M,MI)
*
*----------------------------------------------------------------------------
*
      RETURN
      END
