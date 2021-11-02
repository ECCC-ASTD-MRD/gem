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
***S/P VPS - INTERPOLATION DE VECTEURS A PARTIR D'UNE GRILLE PS


      SUBROUTINE VPS(S,PSI,LAMBDA,PHI,U,V,LI,LJ,PI,PJ,D60,DGRW,HEM)
      REAL U(LI,LJ), V(LI,LJ)
*
*AUTEUR   - M.VALIN DRPN 1981
*
*OBJECT(VPS)
*            INTERPOLATION DE VECTEURS A PARTIR D'UNE GRILLE PS
*
*LANGAGE  - FORTRAN
*
*LIBRAIRIES
*         SOURCE  RMNSOURCELIB,ID=RMNP     DECK=VPS
*         OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL VPS(S,PSI,LAMBDA,PHI,U,V,LI,LJ,PI,PJ,D60,DGRW,HEM)
*
*ARGUMENTS
*   OUT   - S      - VITESSE DU VENT AU POINT DEMANDE
*   OUT   - PSI    - DIRECTION DU VENT
*   IN    - LAMBDA - LONGITUDE DU POINT  (REEL)
*   IN    - PHI    - COLATITUDE(LATITUDE-90.0) DU POINT (REEL)
*   IN    - U,V    - COMPOSANTES U ET V DU VENT (LI,LJ)
*   IN    - LI     - PREMIERE DIMENSION DES CHAMPS U ET V (PS)
*   IN    - LJ     - DEUXIEME DIMENSION DES CHAMPS U ET V (PS)
*   IN    - PI     - COORDONNEE X DU POLE  (REEL)
*   IN    - PJ     - COORDONNEE Y DU POLE  (REEL)
*   IN    - D60    - DISTANCE EN METRES ENTRE LES POINTS DE LA GRILLE
*                    A 60 DEGRES DE LATITUDE
*   IN    - DGRW   - ANGLE ENTRE L'AXE X ET LE MERIDIEN DE GRENWICH
*   IN    - HEM    - CODE D'HEMISPHERE  (1=NORD , 2=SUD)
*
*---------------------------------------------------------------------------
*
      REAL LAMBDA
      integer HEM

#include "pi.cdk"
      ARBPSI = 99999.0


      UFIN = SPS(LAMBDA,PHI,U,LI,LJ,PI,PJ,D60,DGRW,HEM)
      VFIN = SPS(LAMBDA,PHI,V,LI,LJ,PI,PJ,D60,DGRW,HEM)


      S = SQRT(UFIN**2 + VFIN**2)


      PSI = ARBPSI
      IF(S.NE.0) PSI=RDTODG*ATAN2(VFIN,UFIN) - DGRW
      IF(PSI.LT.0.0) PSI = PSI + 360.0


      RETURN
      END
