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
***S/P VALIDE - VERIFIER LA VALEUR DE ICHK
*
      LOGICAL FUNCTION VALIDE(NOM, ICHK, MIN, MAX)
      character (len=*) :: NOM
      integer :: ICHK, MIN, MAX

*
*AUTEUR   - P. SARRAZIN  - AVRIL 82
*REVISION 001  M. LEPINE     -  NOV  94   TRADUCTION RATFOR A fortran
*
*LANGAGE  - fortran
*
*OBJET(VALIDE)
*         - VERIFIER SI ICHK => MIN OU =< MAX
*         - MESSAGE SI ICHK N'EST PAS DANS LES LIMITES
*
*LIBRAIRIES
*         - SOURCE  INTPRSOURCELIB,ID=RMNP   DECK=VALIDE
*         - OBJET   INTPRLIB,ID=RMNP
*
*ARGUMENTS
*   IN    - NOM   - NOM DE LA VARIABLE EMPLOYE PAR LA ROUTINE
*   IN    - ICHK  - VALEUR DE LA VARIABLE POUR VERIFICATION
*   IN    - MIN   - VALEUR MINIMUM DE ICHK
*   IN    - MAX   - VALEUR MAXIMUM DE ICHK
*
*MODULES - SCINT - UVINT
*
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
      IF (ICHK.LT.MIN .OR. ICHK.GT.MAX) THEN
        WRITE(6,600) NOM,ICHK,MIN,MAX
      ENDIF
 600  FORMAT("MAUVAISE VALEUR POUR",A10,"VALEUR=",I10,"MINIMUM=",
     %        I10,"MAXIMUM=",I10)
*
      VALIDE = (ICHK.GE.MIN) .AND. (ICHK.LE.MAX)
*
      RETURN
      END
