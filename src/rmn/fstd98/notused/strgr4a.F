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
#include <ftnmacros.hf>
*** S/P STRGR4A    TRANSFORMER UNE CHAINE DE CAR EN R4A
*
      SUBROUTINE STRGR4A(STRG,R4A,POSDEB,POSFIN)
      CHARACTER * (*) STRG
      INTEGER R4A(0:*)
      INTEGER POSDEB,POSFIN,POS,LNG,MAXCAR,MOT,CARPOS,SHFCNT,C,I,MSK
*
*AUTEUR     M. LEPINE  MARS 85
*
*LANGAGE    RATFOR
*
*OBJET(STRGR4A)
*     TRANSFORMER UNE CHAINE DE CARACTERES EN MOT DE TYPE R4A
*  (4 CAR. DE 8 BITS ALIGNES DANS LA PARTIE DROITE DU MOT).
*
*ARGUMENTS
*  IN     STRG    CHAINE DE CARACTERES A TRANSFORMER
*  OUT    R4A     TABLEAU D'ENTIERS QUI CONTIENDRA LES CAR R4A
*  IN     POSDEB  POSITION DU PREMIER CARACTERE DANS R4A
*  IN     POSFIN  POSITION DU DERNIER CARACTERE DANS R4A
*
**
*

      LNG = LEN(STRG)
      MAXCAR = MIN(LNG,POSFIN-POSDEB+1)
      POS = POSDEB
      DO 23000 I = 1,MAXCAR
         READ(STRG(I:I),'(A1)') C
         C = RSHIFT(C, BITMOT - BITCAR)
         MOT = POS / 4
         CARPOS = MOD(POS,4)
         SHFCNT = (3-CARPOS) * 8

         MSK = LXSHIFT(255, SHFCNT)

         R4A(MOT) = IAND(R4A(MOT),COMPL(MSK))
         R4A(MOT) = IOR(R4A(MOT),LXSHIFT(C,SHFCNT))

         POS = POS + 1

23000 CONTINUE 
      RETURN
      END
