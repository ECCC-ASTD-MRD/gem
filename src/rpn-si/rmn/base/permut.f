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
***S/P PERMUT - ROTATION D'UNE MATRICE AUTOUR DE LA LIGNE
*               DU MILIEU
      SUBROUTINE PERMUT(Z,NI,NJ)
      REAL Z(NI,NJ)
*
*AUTEUR   - M. VALIN  -  FEV 82
*
*REVISION 001  C. THIBEAULT  -  MARS 83  CONVERSION AU CODE CRAY
*REVISION 002  M. LEPINE     -  NOV  94   TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(PERMUT)
*         - ROTATION D'UNE MATRICE AUTOUR DE LA LIGNE DU MILIEU.
*           CETTE ROUTINE EST UTILISEE POUR RE=ARRANGER LES
*           DONNEES DANS UN CHAMP. EX: POUR CONTOURER, ON A
*           DES DONNEES DANS L'ORDRE SUIVANT, DU SUD AU NORD
*           ALORS QUE POUR L'INTERPOLATION ELLES DOIVENT ETRE
*           ETALEES DU NORD AU SUD.
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=PERMUT
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL PERMUT(Z,NI,NJ)
*
*ARGUMENTS
* IN/OUT  - Z  - CHAMP QUI SUBIT LA ROTATION
*   IN    - NI - PREMIERE DIMENSION DE Z
*   IN    - NJ - DEUXIEME DIMENSION DE Z
*
*-----------------------------------------------------------------------
*


      NCC = NJ/2


      DO 20 J=1,NCC 

        DO 10 I=1,NI 
          T = Z(I,NJ+1-J)
          Z(I,NJ+1-J) = Z(I,J)
          Z(I,J) = T
 10     CONTINUE

 20   CONTINUE

      RETURN
      END
