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
!** S/R SAUVDEZ CONTROLE LA PORTEE DES DIRECTIVES APRES COPIE
      SUBROUTINE SAUVDEZ
      use ISO_C_BINDING
      use configuration
      IMPLICIT NONE 
      include 'excdes.inc'
!
!AUTEURS
!VERSION ORIGINALE -   Y. BOURASSA NOV 89
!REVISION      001 -   Y. BOURASSA AVR 92 ANNULE LE ZAP SI SAUV=0
!REVISION      002 -   M. Valin mai 2014 utilisation des fonctions des fichiers standard pour
!                                        la gestion des requetes
!REVISION      003 -   M. Valin aout 2016 correction d'un bug de logique 
!REVISION      004 -   M. Lepine nov 2016 corrections, remettre le compteur de directives a zero
!                                         ainsi que le compteur d'enregistrement 
!
!LANGUAGE   - FTN77 
!
      EXTERNAL FSTCVT, ZAP
!
!*
      INTEGER  FSTCVT, I, J, N
      integer :: status

!     CERTAINS CRITERES DE SELECTION DOIVENT-ILS RESTER VALIDES ?
!     sauv = 0 : on ne conserve rien, et on annulle le 'ZAP'
!     sauv > 0 : on conserver les criteres sauv, sauv+1, ...
!     editfst numerote ses jeux de criteres a partir de 1
!     les fichiers standard les numerotent a partir de 0
!     d'ou le sauv -1 (compatibilite arriere des directives)
      IF( SAUV .LE. 0) THEN
         NP = 1  ! simuler directive ZAP(-1) en mettant NP (nombre de parametres READLX) a 1
         CALL ZAP( -1 )
         do N=0,MAX_REQUETES-1   ! on annulle TOUTES les requetes
           status = f_requetes_reset(N,0,0,0,0,0,0,0)  ! annuller le jeu de criteres N
         enddo
      ELSE
        do N=0,min(SAUV-2,MAX_REQUETES-1)  ! on annulle editfst 1 a sauv -1  ( 0 a sauv -2 pour le code C)
          status = f_requetes_reset(N,0,0,0,0,0,0,0)  ! annuller le jeu de criteres N
        enddo
      ENDIF
  
      NREQ = 0
      COPIES = 0
      RETURN
      END 
