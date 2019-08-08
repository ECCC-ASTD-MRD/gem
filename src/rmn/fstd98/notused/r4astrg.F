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
*.S R4ASTRG ROUTINE
*** S/P R4ASTRG     TRANSFORMER DU R4A EN CARACTERES
*
       SUBROUTINE R4ASTRG(STRG,R4A,POSDEB,POSFIN)
       IMPNONE
       CHARACTER * (*) STRG
       INTEGER R4A(0:*)
       INTEGER POSDEB,POSFIN,POS,LNG,MAXCAR,MOT,CARPOS,SHFCNT,C,I,MSK
*
*AUTEUR     M. LEPINE  MARS 85
*
*LANGAGE    RATFOR
*
*OBJET(R4ASTRG)
*     TRANSFORMER N CARACTERES CONTENUS DANS DES MOT R4A (4 CAR. DE 8 BITS
*  ALIGNES DANS LA PARTIE DROITE DU MOT) EN N CARACTERES CONTENUS DANS UNE
*  CHAINE DE CARACTERES AU SENS FORTRAN.
*
*ARGUMENTS
*  IN     R4A     TABLEAU D'ENTIERS CONTENANT LES CARACTERES DE TYPE R4A
*  OUT    STRG    TABLEAU DE CARACTERES TRANSFORMES
*  IN     POSDEB  POSITION DU PREMIER CARACTERE DANS R4A
*  IN     POSFIN  POSITION DU DERNIER CARACTERE DANS R4A
*
**


       POS = POSDEB
       DO 23000 I = 1,POSFIN-POSDEB+1
          MOT = POS / 4
          CARPOS = MOD(POS,4)
          SHFCNT = (3- CARPOS) * 8
          C = IAND(255, RSHIFT(R4A(MOT), SHFCNT))
          WRITE(STRG(I:I),'(A1)') LSHIFT(C, BITMOT - BITCAR)
          POS = POS + 1
23000  CONTINUE
*
       RETURN
       END
