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
***S/P VPOLAGR - INTERPOLATION DE LAGRANGE
*
      SUBROUTINE VPOLAGR(FF,ZZ,F,Z,NA,NPLUS1)
      REAL FF(NA), F(NA,NPLUS1), Z(NPLUS1)
*AUTEUR   - M. VALIN,  RPN, MAR 81
*
*REVISION 001  C. THIBEAULT  -  MARS 83  CONVERSION AU CODE CRAY
*REVISION 002  M. LEPINE     -  NOV  94   TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(VPOLAGR)
*         - ETANT DONNE F(*,J) = FONCTION DE Z(J), ON VEUT
*           EVALUER FF(*) QUI EST EGALE A LA VALEUR DE LA
*           FONCTION A Z = ZZ
*
*LIBRAIRIES
*         - SOURCE  RMNSOURCELIB,ID=RMNP     DECK=VPOLAGR
*         - OBJET   RMNLIB,ID=RMNP
*
*APPEL    - CALL VPOLAGR(FF,ZZ,F,Z,NA,NPLUS1)
*
*ARGUMENTS
*   OUT   - FF     - CHAMP DE SORTIE (NA)
*         - ZZ     - COORDONNEE A LAQUELLE ON VEUT INTERPOLER
*         - F      - CHAMP SOURCE (NA,NPLUS1)
*         - Z      - COORDONNEE VERTICALE (NPLUS1)
*         - NA     - DIMENSION HORIZONTALE DE FF, F
*         - NPLUS1 - DEGRE MAXIMUM DES POLYNOMES DE LAGRANGE + 1
*

      
      IF (NA .EQ. 1) THEN
         FF(1) = POLAGR(ZZ,F,Z,NPLUS1)
         RETURN
      ENDIF
      
      
      DO 10 I=1,NA
         FF(I) = 0.0
 10   CONTINUE
      
      
      DO 100 I =1,NPLUS1
         POLY = 1.0
         
         DO 20 J =1,NPLUS1
            IF (J .ne. I) THEN
               POLY = POLY*(ZZ-Z(J))/(Z(I)-Z(J))
            ENDIF
 20      CONTINUE
         
         DO 30 II=1,NA
            FF(II) = POLY * F(II,I) + FF(II)
 30      CONTINUE
         
 100  CONTINUE
      RETURN
      END
