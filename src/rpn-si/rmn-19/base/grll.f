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
***S/P GRLL - CALCULE LA LATITUDE ET LA LONGITUDE D'UNE
*             GRILLE LAT. LON.


      SUBROUTINE GRLL(XLAT,XLON,NI,NJ,XLA0,XLO0,DLA0,DLO0)
      REAL XLAT(NI,NJ), XLON(NI,NJ)
*
*AUTEUR   - M. VALIN  -  JAN 82
*
*REVISION 001  C. THIBEAULT  -  MARS 83   CONVERSION AU CODE CRAY
*REVISION 002  M. LEPINE     -  NOV  94   TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(GRLL)
*         - CALCULE LA LATITUDE ET LA LONGITUDE DE CHACUN
*           DES POINT D'UNE GRILLE LAT.-LON.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=GRLL
*         - RMNLIB,ID=RMNP
*
*APPEL    - CALL GRLL(XLAT,XLON,NI,NJ,XLA0,XLO0,DLA0,DLO0)
*
*ARGUMENTS
*   OUT   - XLAT - CHAMP DE LATITUDES (NI,NJ). (NIVEAU 1 OU 2)
*   OUT   - XLON - CHAMP DE LONGITUDES (NI,NJ). (NIVEAU 1 OU 2)
*   IN    - NI   - NOMBRE DE POINTS PAR CERCLE DE LATITUDE.
*   IN    - NJ   - NOMBRE DE CERCLES DE LATITUDE.
*   IN    - XLA0 - LATITUDE DU COIN INFERIEUR GAUCHE (DEGRES).
*   IN    - XLO0 - LONGITUDE DU COIN INFERIEUR GAUCHE (DEGRES).
*   IN    - DLA0 - ESPACEMENT EN LATITUDE (DEGRES).
*   IN    - DLO0 - ESPACEMENT EN LONGITUDE (DEGRES).
*
*----------------------------------------------------------------------
*
*
      DO 20 J=1,NJ
        XLA = XLA0 +(J-1) * DLA0

        DO 10 I=1,NI
          XLAT(I,J) = XLA
          XLON(I,J) = AMOD(XLO0 +(I-1) * DLO0,360.0)
 10     CONTINUE
 20   CONTINUE

      RETURN
      END
