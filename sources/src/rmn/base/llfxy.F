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
***S/R LLFXY - COMPUTES LATITUDE AND LONGITUDE OF A POINT
*
      SUBROUTINE LLFXY(DLAT,DLON,X,Y,D60,DGRW,NHEM)
*
*AUTHOR   - J.D. HENDERSON  -  FEB 75
*
*REVISION 001   C. THIBEAULT  -  NOV 79  DOCUMENTATION
*
*LANGUAGE - FORTRAN
*
*OBJECT(LLFXY)
*         - COMPUTES LATITUDE AND LONGITUDE OF A POINT IN A POLAR
*           STEREOGRAPHIC GRID FROM CO-ORDINATES IN THE GRID MEASURED
*           FROM THE POLE.
*
*LIBRARIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=LLFXY
*         - OBJECT  RMNLIB,ID=RMNP
*
*USAGE    - CALL LLFXY(DLAT,DLON,X,Y,D60,DGRW,NHEM)
*
*ARGUMENTS
*   OUT   - DLAT - LATITUDE IN DEGREES (-90 TO +90, POSITIVE N).
*         - DLON - LONGITUDE IN DEGREES (-180 TO +180, POSITIVE E).
*   IN    - X    - X-CO-ORDINATE OF THE POINT AS MEASURED WITH POLE
*                  AS ORIGIN
*         - Y    - Y-CO-ORDINATE OF THE POINT AS MEASURED WITH POLE
*                  AS ORIGIN
*         - D60  - GRID LENGTH (IN METRES) OF THE POLAR STEREOGRAPHIC
*                  GRID AT 60 DEGREES
*         - DGRW - ORIENTATION OF GREENWICH MERIDIAN WITH RESPECT TO
*                  THE GRID (IN DEGREES)
*         - NHEM - 1 FOR NORTHERN HEMISPHERE
*                  2 FOR SOUTHERN HEMISPHERE
*
*NOTES    - THE COMPANION ROUTINE XYFLL, WHICH COMPUTES THE GRID
*           CO-ORDINATES GIVEN THE LATITUDE AND LONGITUDE, IS ALSO
*           AVAILABLE.
*
*----------------------------------------------------------------------------
*
*     * 1.866025=(1+SIN60),   6.371E+6=EARTH RADIUS IN METERS.
#include "pi.cdk"
*
      RE=1.866025*6.371E+6/D60
      RE2=RE**2
*
*     * IF POINT IS AT POLE SET COORD TO (0.,90.).
*
      DLAT=90.
      DLON=0.
      IF(X.EQ.0. .AND. Y.EQ.0.) GO TO 39
*
*     * CALCULATE LONGITUDE IN MAP COORDINATES.
*
      IF(X.EQ.0.) DLON=SIGN(90.,Y)
      IF(X.NE.0.) DLON=ATAN(Y/X)*RDTODG
      IF(X.LT.0.) DLON=DLON+SIGN(180.,Y)
*
*     * ADJUST LONGITUDE FOR GRID ORIENTATION.
*
      DLON=DLON-DGRW
      IF(DLON.GT.+180.) DLON=DLON-360.
      IF(DLON.LT.-180.) DLON=DLON+360.
*
*     * CALCULATE LATITUDE.
*
      R2=X**2+Y**2
      DLAT=(RE2-R2)/(RE2+R2)
      DLAT= ASIN(DLAT)*RDTODG
*
*     * CHANGE SIGNS IF IN SOUTHERN HEMISPHERE.
*
   39 IF(NHEM.EQ.2) DLAT=-DLAT
      IF(NHEM.EQ.2) DLON=-DLON
*
*-------------------------------------------------------------------------
*
      RETURN
      END
