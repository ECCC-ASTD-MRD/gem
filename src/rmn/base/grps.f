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
***S/P GRPS - CALCULE LA LATITUDE ET LA LONGITUDE D'UNE
*             GRILLE POLAIRE STEREOGRAPHIQUE.


      SUBROUTINE GRPS(XLAT,XLON,NI,NJ,PI,PJ,D60,DGRW,HEM)
      REAL XLAT(NI,NJ), XLON(NI,NJ)
*
*AUTEUR   - M. VALIN  -  JAN 82
*
*REVISION 001  C. THIBEAULT  -  MARS 83  CONVERSION AU CODE CRAY
*REVISION 002  M. LEPINE     -  NOV  94   TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(GRPS)
*         - CALCULE LA LATITUDE ET LA LONGITUDE DE CHACUN
*           DES POINTS D'UNE GRILLE POLAIRE STEREOGRAPHIQUE.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=GRPS
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL GRPS(XLAT,XLON,NI,NJ,PI,PJ,D60,DGRW,HEM)
*
*ARGUMENTS
*   OUT   - XLAT - CHAMP DE LATITUDES (NI,NJ). (NIVEAU 1 OU 2)
*   OUT   - XLON - CHAMP DE LONGITUDES (NI,NJ). (NIVEAU 1 OU 2)
*   IN    - NI   - PREMIERE DIMENSION DES CHAMPS XLAT,XLON.
*   IN    - NJ   - DEUXIEME DIMENSION DES CHAMPS XLAT,XLON.
*   IN    - PI   - COORDONNEE X DU POLE (REEL).
*   IN    - PJ   - COORDONNEE Y DU POLE (REEL).
*   IN    - D60  - DISTANCE EN METRES ENTRE LES POINTS DE LA
*                  GRILLE A 60 DEGRES DE LATITUDE.
*   IN    - DGRW - ANGLE ENTRE L'AXE X ET LE MERIDIEN DE
*                  GREENWICH.
*   IN    - HEM  - CODE D'HEMISPHERE (1=NORD, 2=SUD).
*
*MODULES  - LLFXY
*
*-----------------------------------------------------------------
*
      INTEGER HEM

      DO 20 J=1,NJ
        Y = J-PJ

        DO 10 I=1,NI
          CALL LLFXY(XLA,XLO,I-PI,Y,D60,DGRW,HEM)
          XLAT(I,J) = XLA
          IF (XLO.LT.0) XLO = XLO + 360.
          XLON(I,J) = XLO
 10     CONTINUE

 20   CONTINUE
      RETURN
      END
