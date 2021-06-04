!/* EDITFST - Collection of useful routines in C and FORTRAN
! * Copyright (C) 1975-2014  Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!**S/R HOLACAR - TRANSFORMER EN VRAIES CHAINES DE CARACTERES
!                LES CARACTERES PASSES VIA DIRECTIVES READLX (STOCKES DANS DES ENTIERS)

      SUBROUTINE HOLACAR(LABEL, LIST, NL, STRING, NC)
  
      IMPLICIT   NONE 
      INTEGER, intent(IN) ::   NL, LIST(NL), STRING(NL*3), NC
      CHARACTER(len=*), intent(OUT) :: LABEL(NL)
!
!AUTEUR   -   Y. BOURASSA  - AVR 91
!REVISION 001 "      "     - JAN 92 
!Revision 002   M. Lepine - mars 98 - extensions pour fstd98
!Revision 003   M. Lepine - Nov  05 - remplacement de fstabt par qqexit
!Revision 004   M. Valin  - Mars 14 - menage et suppression de l'option entiers 64 bits
!LANGAGE  - FTN77
!
!ARGUMENTS
!SORTIE   - LABEL   - ETIKETTES 
!ENTREE   - LIST    - CHAMP RETOURNEE PAR ARGDOPE. (readlx)
!   "     - NL      - DIMENSION DE LABEL ET LIST (nombre de strings a decoder).
!   "     - STRING  - CHAINE DE CARACTERES A DECODER.(tasse dans des entiers)
!   "     - NC      - NOMBRE DE CARACTERES ALLOUE POUR LABEL. (on suppose <=12 dans le code)
!
      EXTERNAL qqexit
      INTEGER  :: I, J, K, L
      integer, parameter :: NCW = 4  ! stocke dans des entiers a 32 bits par readlx
      character *12 temp12
      INTEGER  X, Y

!     PASSE DE HOLLERITH A CARACTERES
      DO 10 K=1,NL
         L = ishft(LIST(K), -16)             ! position du debut d'extraction dans string
         I = IAND(255, ishft(LIST(K), -8))   ! nombre de caracteres a extraire
         IF(I .GT. NC) THEN
            WRITE(6,*)' LIMITE DE',NC,' CARACTERES DEPASSEE :',I
            CALL qqexit(40)
         ENDIF
         IF(IAND(255, LIST(K)) .NE. 3) THEN  ! ce ne sont pas des caracteres (dixit readlx)
            WRITE(6,*)' ARGUMENT PAS DE TYPE CARACTERE'
            CALL qqexit(41)
         ENDIF
         J = (I+NCW-1)/NCW             ! nombre de mots de 32 bits
         WRITE(temp12, '(3A4)') (STRING(J),J=L,L+J-1)
         LABEL(K) = temp12(1:I)
   10    CONTINUE 

      RETURN
      END 




















