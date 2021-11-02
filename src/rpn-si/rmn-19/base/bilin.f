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
***FUNCTION BILIN - INTERPOLATEUR BI-LINEAIRE
*
      FUNCTION BILIN(Z,NI,NJ,X,Y,OK)
      REAL Z(NI,NJ)

*
*AUTEUR   - M. VALIN - JUIN 81
*
*LANGAGE  - fortran
*
*OBJET(BILIN)
*         - CALCULE LA VALEUR AUX POINTS X, Y (INDICES FRACTIONNAIRES)
*           DANS UN CHAMP Z DE DIMENSION NI*NJ. L'INTERPOLATION
*           EST FAITE DE FACON BI-LINEAIRE.
*
*LIBRAIRIES
*         - SOURCE   RMNSOURCELIB,ID=RMNP
*         - OBJET    RMNLIB,ID=RMNP
*
*APPEL    - R = BILIN(Z,NI,NJ,X,Y,OK)
*
*ARGUMENTS
*         - Z  - CHAMP A PARTIR DUQUEL ON CALCULE LA COUPE (NI,NJ)
*         - NI - X-DIMENSION DU CHAMP Z
*         - NJ - Y-DIMENSION DU CHAMP Z
*         - X  - COORDONNEE X DU POINT A CALCULER
*         - Y  - COORDONNEE Y DU POINT A CALCULER
*         - OK - CLE LOGIQUE, TRUE  = LES COORDONNEES SONT DANS LES
*                                     LIMITES DU CHAMP Z
*                             FALSE = LES COORDONNEES NE SONT PAS
*                                     DANS LES LIMITES DU CHAMP Z
*
*----------------------------------------------------------------------
*
      LOGICAL OK


      OK = OK.AND.X.GE.1.AND.X.LE.NI.AND.Y.GE.1.AND.Y.LE.NJ
      IF (.NOT.OK) BILIN = 0.
      IF (.NOT.OK) RETURN


      I = AINT(X)
      J = AINT(Y)


      IF (I.EQ.NI) I=NI-1
      IF (J.EQ.NJ) J=NJ-1
      DX = X-I
      DY = Y-J


      Z1 = Z(I,J) + DY * (Z(I,J+1) - Z(I,J))
      Z2 = Z(I+1,J) + DY * (Z(I+1,J+1) -  Z(I+1,J))
      BILIN = Z1 + DX * (Z2-Z1)


      RETURN
      END
