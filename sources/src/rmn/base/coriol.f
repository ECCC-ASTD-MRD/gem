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
***S/R CORIOL - COMPUTES THE CORIOLIS PARAMETER
*
      SUBROUTINE CORIOL (FF,D60,PI,PJ,NI,NJ)
      REAL FF(NI,NJ)

*
*AUTHOR   - HERSCHEL MITCHELL - FEB 1975
*
*REVISION 001  C. THIBEAULT  -  OCT 1979   DOCUMENTATION
*REVISION 002  C. THIBEAULT  -  MAR 1983   CONVERSION AU CODE CRAY
*
*OBJECT(CORIOL)
*         - COMPUTES THE CORIOLIS PARAMETER AT ALL POINTS OF A GRID OF
*           INIFORM MESH-LENGTH. THE GRID IS ASSUMED TO BE SUPERIMPOSED
*           ON A POLAR STEREOGRAPHIC PROJECTION OF THE EARTH TRUE AT
*           LATITUDE 60 DEGREES.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=CORIOL
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL CORIOL (FF,D60,PI,PJ,NI,NJ)
*
*ARGUMENTS
*   OUT   - FF  - ARRAY CONTAINING THE CORIOLIS PARAMETER
*   IN    - D60 - GRID-LENGTH AT 60 DEG. LAT. IN METRES
*         - PI  - X-COORDINATE OF POLE RELATIVE TO BOTTOM
*                 LEFT-HAND CORNER OF GRID
*         - PJ  - Y-COORDINATE OF POLE RELATIVE TO BOTTOM
*                 LEFT-HAND CORNER OF GRID
*         - NI  - X-DIMENSION
*         - NJ  - Y-DIMENSION
*
*NOTES    - SIN(LAT) = (RE2-R2) / (RE2+R2)
*         - RE = (EARTH RADIUS) * (1+SIN(60))/D60
*         - MEAN EARTH RADIUS = 6371 KM.
*         - CF = 2. * OMEGA
*         - NOTE THAT THE POLE NEED NOT BE LOCATED ON A GRID POINT
*
*-------------------------------------------------------------------------------
*
      DATA CF  /1.458E-4/
*
*-------------------------------------------------------------------------------
*
      RE = 1.866025*6.371E6/D60
      RE2= RE*RE
*
      DO 10 J=1,NJ
      Y = FLOAT(J) - PJ
      Y2= Y*Y
*
      DO 10 I=1,NI
      X = FLOAT(I) - PI
      R2= X*X + Y2
      FF(I,J) = CF * (RE2-R2) / (RE2+R2)
   10 CONTINUE
*
*-------------------------------------------------------------------------------
*
      RETURN
      END
