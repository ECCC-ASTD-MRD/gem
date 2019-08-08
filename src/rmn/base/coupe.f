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
***S/R COUPE - COUPE TRANSVERSALE D'UN CHAMP BIDIMENSIONNEL
*
      SUBROUTINE COUPE(R,NP,Z,NI,NJ,X1,Y1,X2,Y2,OK)
      REAL R(NP), Z(NI,NJ)

*
*AUTEUR   - M. VALIN - JUIN 81
*
*LANGAGE  - fortran
*
*OBJET(COUPE)
*         - CALCULE UNE COUPE DE NP POINTS DANS LE CHAMP Z
*         - A PARTIR DES INDICES X1, Y1 AUX POINTS X2, Y2.
*           X1, Y1, X2 ET Y2 SONT DES INDICES FRACTIONNAIRES
*
*LIBRAIRIES
*         - SOURCE   RMNSOURCELIB,ID=RMNP
*         - OBJET    RMNLIB,ID=RMNP
*
*APPEL    - CALL COUPE(R,NP,Z,NI,NJ,X1,Y1,X2,Y2,OK)
*
*ARGUMENTS
*         - R  - CHAMP QUI CONTIENT LES RESULTATS DE LA COUPE (NP)
*         - NP - NOMBRE DE POINTS DANS LA COUPE
*         - Z  - CHAMP A PARTIR DUQUEL ON CALCULE LA COUPE (NI,NJ)
*         - NI - X-DIMENSION DU CHAMP Z
*         - NJ - Y-DIMENSION DU CHAMP Z
*         - X1 - COORDONNEE X DU COIN GAUCHE DE LA COUPE
*         - Y1 - COORDONNEE Y DU COIN GAUCHE DE LA COUPE
*         - X2 - COORDONNEE X DU COIN DROIT DE LA COUPE
*         - Y2 - COORDONNEE Y DU COIN DROIT DE LA COUPE
*         - OK - CLE LOGIQUE, TRUE  = LES COORDONNEES SONT DANS LES
*                                     LIMITES DU CHAMP Z
*         -                   FALSE = LES COORDONNEES NE SONT PAS DANS
*                                     LES LIMITES DU CHAMP Z
*
*-----------------------------------------------------------------------
*
      LOGICAL OK


      OK = .TRUE.
      DELTAX = (X2-X1) / (NP-1)
      DELTAY = (Y2-Y1) / (NP-1)


      DO 50 I=1,NP
      X = X1 + (I-1) * DELTAX
      Y = Y1 + (I-1) * DELTAY
      R(I) = BILIN(Z,NI,NJ,X,Y,OK)
   50 CONTINUE


      RETURN
      END
