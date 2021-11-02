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
*
***S/P VLL - INTERPOLATION DE VECTEURS A PARTIR D'UNE GRILLE LL
*
      SUBROUTINE VLL(S,PSI,LAMBDA,PHI,U,V,LI,LJ,XLA0,XLO0,DLAT,DLON)
      REAL U(LI,LJ), V(LI,LJ)
*
*AUTEUR   - M.VALIN DRPN 1981
*
*OBJECT(VLL)
*           INTERPOLATION DE VECTEURS A PARTIR D'UNE GRILLE LL
*
*LANGAGE  - FORTRAN
*
*LIBRAIRIES
*         SOURCE  RMNSOURCELIB,ID=RMNP     DECK=VLL
*         OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL VLL(S,PSI,LAMBDA,PHI,U,V,LI,LJ,XLA0,XLO0,DLAT,DLON)
*
*ARGUMENTS
*   OUT   - S      - VITESSE DU VENT AU POINT DEMANDE
*   OUT   - PSI    - DIRECTION DU VENT
*   IN    - LAMBDA - LONGITUDE DU POINT  (REEL)
*   IN    - PHI    - COLATITUDE(LATITUDE-90.0) DU POINT (REEL)
*   IN    - U,V    - COMPOSANTES U ET V DU VENT (LI,LJ)
*   IN    - LI     - PREMIERE DIMENSION DES CHAMPS U ET V (PS)
*   IN    - LJ     - DEUXIEME DIMENSION DES CHAMPS U ET V (PS)
*   IN    - XLA0   - LATITUDE DU COIN INFERIEUR GAUCHE (DEGRES).
*   IN    - XLO0   - LONGITUDE DU COIN INFERIEUR GAUCHE (DEGRES).
*   IN    - DLAT   - ESPACEMENT EN LATITUDE (DEGRES).
*   IN    - DLON   - ESPACEMENT EN LONGITUDE (DEGRES).
*
*---------------------------------------------------------------------------
*
      REAL LAMBDA
#include "pi.cdk"
      ARBPSI = 99999.0


      UFIN = SLL(LAMBDA,PHI,U,LI,LJ,XLA0,XLO0,DLAT,DLON)
      VFIN = SLL(LAMBDA,PHI,V,LI,LJ,XLA0,XLO0,DLAT,DLON)


      S = SQRT(UFIN**2 + VFIN**2)


      PSI = ARBPSI
      SS = -SIN (LAMBDA*DGTORD)
      CC = COS (LAMBDA*DGTORD)
      UF = SS*UFIN - CC*VFIN
      VF = CC*UFIN + SS*VFIN
      IF (S.NE.0) PSI = RDTODG * ATAN2(VF,UF)
      IF(PSI.LT.0.0) PSI = PSI + 360.0


      RETURN
      END
