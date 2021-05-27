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
***FONCTION POLAGR - INTERPOLATION DE LAGRANGE
*
      FUNCTION POLAGR(ZZ,F,Z,NPLUS1)
      REAL Z(NPLUS1), F(NPLUS1)
*
*AUTEUR   - M. VALIN, RPN, MAR 81
*REVISION 001  M. LEPINE     -  NOV  94   TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(POLAGR)
*         - TROUVE LA VALEUR AU POINT ZZ DE L'INTERPOLATION
*           DE LAGRANGE DE DEGRE N, (NPLUS1 = N+1), QUI
*           PREND LES VALEURS DE F AUX POINTS DE Z.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=POLAGR
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - YY=POLAGR(ZZ,F,Z,NPLUS1)
*
*ARGUMENTS
*         - ZZ     - POINT POUR LEQUEL ON CHERCHE F(ZZ)
*         - F      - VALEURS DE LA FONCTION
*         - Z      - LISTE DE COORDONNEES
*         - NPLUS1 - DEGRE MAXIMUM DES POLYNOMES DE LAGRANGE + 1
*


      POLAGR = 0.0

      DO 20 I=1,NPLUS1
      POLY = 1.0


      DO 10 J=1,NPLUS1
        IF (J .NE. I) THEN
          POLY = POLY*(ZZ-Z(J))/(Z(I)-Z(J))
        ENDIF
 10   CONTINUE

      POLAGR = POLAGR + F(I)*POLY

 20   CONTINUE
      RETURN
      END
