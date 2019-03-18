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
!** S/R DMPDES MET EDITFST EN MODE XPRES SI L'USAGER VEUT TOUT LE FICHIER
!       IMPRIME L'INTERPRETATION DES DIRECTIVES DESIRE SI EN MODE DEBUG
      SUBROUTINE DMPDES
      use ISO_C_BINDING
      use configuration
      IMPLICIT NONE

      include 'excdes.inc'
  
!     AUTEUR - Y. R. BOURASSA AVR 86
!Revision 002   M. Lepine - mars 98  - extensions pour fstd98
!Revision 003   M.Valin   - mai 2014 - utiliser les fonctions des fichiers standard pour
!                                      l'impression des directives
!     LANGUAGE FTN77
!*  
      INTEGER I, J, K, L
  

      ESAIS = .TRUE.   !   ON  NOTE UNE TENTATIVE DE COPIE
      DM1   = .TRUE.
  
!     AUCUNE DIRECTIVE ?
      IF(NREQ .EQ. 0) THEN  ! NREQ = no de la directive courante
         XPRES = .NOT.SCRI  ! si pas de criteres supplementaires, on veut vraiment tout
         IF( DEBUG ) WRITE(6,*)'INFO: ON COPIE TOUT LE FICHIER '
         RETURN
      ENDIF
  
!     IMPRESSION DES REQUETES
!     appeler la routine appropriee des fichiers standard pour ce faire (si DEBUG)
      IF( DEBUG ) call Dump_Request_Table  
      return
!     il faudrait idealement pouvoir demander aux fichiers standard s'il y a des directives en vigueur
      END 
