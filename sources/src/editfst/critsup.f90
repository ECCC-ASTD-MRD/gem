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
!**S/R CRITSUP DEFIFIT LES CRITERES SUPLEMENTAIRES DE SELECTION
      SUBROUTINE CRITSUP(NI, NJ, NK, GRID, IG1, IG2, IG3, IG4)
      use ISO_C_BINDING
      use configuration
      IMPLICIT NONE 
      include 'excdes.inc'
      INTEGER     NI, NJ, NK, GRID, IG1, IG2, IG3, IG4
!
!AUTEUR       AUTEUR YVON BOURASSA . JAN 90
!                      "      "      OCT 90 VERSION QLXINS
!Revision 002   M. Lepine - mars 98 - extensions pour fstd98
!Revision 003   M. Lepine - Nov 2002 - passage du bon type d'arguments pour fstcvt
!Revision 004   M. Valin  - Mai 2014 - remplacement des cdk par un module
!
!LANGAGE      FTN77 
!
!ARGUMENTS
! ENTRE     - NI,NJ,NK - DIMENSIONS DE LA GRILLE
! ENTRE     - GRID     - TYPE DE GRILLE 
! ENTRE     - IG1@AG4  - DESCRIPTEURS DE LA GRILLE
!
!OBJET(CRITSUP)
!              DEFIFIT DES DETAILS SUR LA GRILLE, ILS SERONT CONSIDERES
!              COMME CLES DE RECHERCHE APPLIQUABLES A TOUS LES DIRECTIVES
!              DESIRE & EXCLURE.
!
! NOTE:        LES DETAILS SERONT AUTOMATIQUEMENT DESACTIVES APRES EXECUTION
!              D'UNE DIRECTIVE SEQCOPI OU STDCOPI.
!
!MODULE
      EXTERNAL FSTCVT
!*
      INTEGER  I, FSTCVT
      CHARACTER*2 TV
      CHARACTER*4 NV
      CHARACTER*12 LBL
  
!     DESACTIVAGE DES CRITERES EN FONCRION
      IG4S = -1
      IG3S = -1
      IG2S = -1
      IG3S = -1
      GTYPS= ' '
      NKS  = -1
      NJS  = -1
      NIS  = -1
  
!     RECUPERATION DES CRITERES FOURNIS PAR LA DIRECTIVE
      GO TO (8, 7, 6, 5, 4, 3, 2, 1) NP   ! la valeur de NP vient de READLX
    1 IG4S = IG4
    2 IG3S = IG3
    3 IG2S = IG2
    4 IG1S = IG1
    5 I = FSTCVT(-1, -1, -1, GRID, NV, TV, LBL, GTYPS, .TRUE.) ! traduire GRID en caracteres, resultat dans GTYPS
    6 NKS  = NK
    7 NJS  = NJ
    8 NIS  = NI
!     DESACTIVER LES CRITERES SUPLEMENTAIRES AVEC "CTITSUP(-1)"
      IF(NP.EQ.1 .AND. NI.EQ.-1) THEN
         SCRI = .FALSE.
         IF( DEBUG ) PRINT*,'CRITERES SUPLEMENTAIRES DESACTIVES'
      ELSE
         IF( DEBUG ) THEN
            IF( SCRI ) THEN
               PRINT*,'CRITERES SUPLEMENTAIRES DE SELECTION MODIFIES' 
            ELSE
               PRINT*,'CRITERES SUPLEMENTAIRES DE SELECTION ACTIVES'
            ENDIF
         ENDIF
         SCRI = .TRUE.  ! on a des criteres supplementaires
      ENDIF
  
      RETURN
      END 
