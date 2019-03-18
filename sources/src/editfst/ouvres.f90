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
!** S/R OUVRES OUVRE UN FICHIER SOURCE
      SUBROUTINE OUVRES( DN ) 
      use configuration
      IMPLICIT      NONE
      CHARACTER*(*) DN(120)
  
!ARGUMENTS
! ENTREE DN    -  LISTE DES NOMS DE FICHIER SOURCE
!
!AUTEUR -      Y. BOURASSA SEP 90
!REVISION 001  "      "    OCT 90 VERSION QLXINS, STANDARD 90
!         002  "      "     "  91 LINK FICHIERS SOURCE
!         003  "      "    FEV 92 CHECK SI FICHIER INEXISTANT
!                                 FNOM-OUVRE-LINK VIA FSTNOL
!                                 CHANGE APPEL A FSTVOI
!         004  "      "    NOV 92 liimite a 15 nombre de fichiers
!         005  M. Lepine   Juil 2005, limite du nombre de fichiers a 120
!         006  M. Valin    Mai  2014, remplacement des comecks par un module
!LANGUAGE FTN77
!
!
!MODULES
      EXTERNAL FSTNOL, FSTVOI, FSTRWD, FERMES, FERMED
      INTEGER  FSTNOL, FSTVOI, FSTRWD, I

!     TRAITEMENT DU FICHIER SOURCE
!     si le fichier source etait ouvert comme fichier destination
      IF( DN(1).EQ.ND .AND. OUVD ) CALL FERMED
      IF( OUVS ) THEN
         IF( DN(1).EQ.NS .AND. NFS.EQ.1 ) RETURN
!        si on change de fichier source
         CALL FERMES
      ENDIF
      I    = FSTNOL(SOURCES, DN, NFS, SNOM)
      NFSO = NFS
      SSEQ = INDEX(SNOM,'SEQ') .GT. 0

!     REFERME LE FICHIER SI 'RND' ET VIDE
      IF(SSEQ .OR. (I .GE. 0)) THEN
         OUVS = .TRUE.
         NS   = DN(1)
         IF( VS ) THEN
            IF( SSEQ ) I = FSTRWD( SOURCES ) 
            IF( INDEX(SNOM,'FTN') .GT. 0) THEN
               I = FSTVOI(SOURCES, 'SEQ')
            ELSE
               I = FSTVOI(SOURCES, 'STD')
            ENDIF
            IF( SSEQ ) I = FSTRWD( SOURCES ) 
         ENDIF
      ELSE
         CALL FERMES
      ENDIF
      
      RETURN
      END 





