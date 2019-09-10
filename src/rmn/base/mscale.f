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
***S/R MSCALE - COMPUTES THE MAP-SCALE FACTOR
*
      SUBROUTINE MSCALE (R,D60,PI,PJ,NI,NJ)
      REAL R(NI,NJ)

*
*AUTHOR   - HERCHELL MITCHELL  -  FEB 75
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*REVISION 002   C. THIBEAULT  -  MAR 83  CONVERSION AU CODE CRAY
*
*LANGUAGE - fortran
*
*OBJECT(MSCALE)
*         - COMPUTES THE MAP-SCALE FACTOR AT ALL POINTS OF A GRID OF
*           UNIFORM MESH-LENGTH. THE GRID IS ASSUMED TO BE SUPERIMPOSED
*           ON A POLAR STEREOGRAPHIC PROJECTION OF THE EARTH TRUE AT
*           LATITUDE 60 DEGREES.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=MSCALE
*         - OBJECT RMNLIB,ID=RMNP
*
*USAGE    - CALL MSCALE(R,D60,PI,PJ,NI,NJ)
*
*ARGUMENTS
*         - R   - ARRAY WHERE THE RESULTS ARE STORED
*         - D60 - GRID-LENGTH AT LATITUDE 60 DEGREES (IN METRES)
*         - PI  - X-CO-ORDINATE OF POLE RELATIVE TO BOTTOM LEFT
*                 CORNER OF GRID
*         - PJ  - Y-CO-ORDINATE OF POLE RELATIVE TO BOTTOM LEFT
*                 CORNER OF GRID
*         - NI  - X-DIMENSION
*         - NJ  - Y-DIMENSION
*
*NOTES    - THE POLE NEED NOT BE LOCATED ON A GRID-POINT.
*
*-------------------------------------------------------------------------------
*
*
*------------------------------------------------------------------------
*
*     *    R = (1+SIN(60))/(1+SIN(LAT))
*     *    SIN(LAT) = (RE2-R2)/(RE2+R2)
*     *    RE=(EARTH RADIUS)*(1+SIN(60))/D60
*     *    MEAN EARTH RADIUS = 6371 KM.
*     *    CM = (1+SIN(60))/(2.*RE2)
*
      RE = 1.866025*6.371E6/D60
      RE2= RE*RE
      CM = 1.866025/(2.*RE2)
*
      DO 10 J=1,NJ
      Y = FLOAT(J) - PJ
      Y2= Y*Y
      DO 10 I=1,NI
      X = FLOAT(I) - PI
      R2= X*X + Y2
      R(I,J)=CM*(R2+RE2)
   10 CONTINUE
*
*------------------------------------------------------------------------
*
      RETURN
      END
