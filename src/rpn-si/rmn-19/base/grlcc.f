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
***S/P GRLCC - CALCULE LA LATITUDE ET LA LONGITUDE D'UNE
*             GRILLE LAMBERT CONFORME CENTREE NORD


      SUBROUTINE GRLCC(XLAT,XLON,NI,NJ,PHI1,PHI2,LAMBDA0,PHI0,DS)
      INTEGER NI, NJ, I, J
      REAL XLAT(NI,NJ),XLON(NI,NJ),PHI1,PHI2,LAMBDA0,PHI0,DS
*
*AUTEUR   - M. VALIN  -  NOV 86
*REVISION 001 M. LEPINE - NOV 94 - TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(GRLCC)
*         - CALCULE LA LATITUDE ET LA LONGITUDE DE CHACUN
*           DES POINTS D'UNE GRILLE LAMBERT CONFORME CENTREE
*           AUTOUR DU POINT LAMBDA0,PHI0.
*
*ARGUMENTS
*   OUT   - XLAT - CHAMP DE LATITUDES (NI,NJ).
*   OUT   - XLON - CHAMP DE LONGITUDES (NI,NJ).
*   IN    - NI   - PREMIERE DIMENSION DES CHAMPS XLAT,XLON.
*   IN    - NJ   - DEUXIEME DIMENSION DES CHAMPS XLAT,XLON.
*   IN    - PHI1 - LATITUDE DU PARALLELE STANDARD LE PLUS AU NORD
*   IN    - PHI2 - LATITUDE DU PARALLELE STANDARD LE PLUS AU SUD
*   IN    - LAMDA0 - MERIDIEN CENTRAL DE LA PROJECTION ET CENTRE DE LA GRILLE
*   IN    - PHI0 - LATITUDE DU POINT CENTRAL DE LA GRILLE
*   IN    - DS   - DISTANCE EN METRES ENTRE LES POINTS DE LA
*                  GRILLE A PHI1 DEGRES DE LATITUDE (OU PHI2).
*
*MODULES
      EXTERNAL LCCFXY, XYFLCC
*
**
      REAL X0,Y0
      CALL XYFLCC(X0,Y0,PHI0,LAMBDA0,PHI1,PHI2,LAMBDA0,DS,1.)
      X0 = (NI+1)*.5 - X0
      Y0 = (NJ+1)*.5 - Y0

      DO 10 I=1,NI
      DO 10 J=1,NI
        CALL LCCFXY(X0+I,Y0+J,XLAT(I,J),XLON(I,J),
     %              PHI1,PHI2,LAMBDA0,DS,1.)
 10   CONTINUE
      RETURN
      END
