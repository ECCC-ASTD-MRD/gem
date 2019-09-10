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
!** FUNCTION FSTNOL OUVRE ET LINKE UNE LISTE DE N FICHIERS
!     L'ENCHAINEMENT SE FAIT SEULEMENT SI LES FICHIERS SONT 'RND'
      FUNCTION FSTNOL(LISTE, NOMS, N, OPTIONS) 
      IMPLICIT NONE
      INTEGER  FSTNOL, N, LISTE(N)
      CHARACTER*(*) NOMS(N), OPTIONS
  
!ARGUMENTS
!SORTIE      - FSTNOL  -  NOMBRE D'ENREGISTREMENTS ENCHAINES
!ENT/SORTIE  - LISTE   -  UNITES ASSIGNES AUX FICHIERS SOURCES
!ENTREE      - NOMS    -  NOMS DES FICHIERS SOURCES
!ENT/SORT    - N       -  DIMENSION DE LISTE ET NOMS
!ENTREE      - OPTIONS -  'RND','SEQ','SEQ/FTN'
!
!AUTEURS
!VERSION ORIGINALE  Y. BOURASSA FEV 92
!REVISION      001  "      "    AVR 92 ENLEVE APPEL A FSTABT
!              002  "      "    MAI 92 SI UN DES FICHIERS INEXISTANT
!                                      FCLOS LES FICHIERS ASSIGNES
!              003  M. Lepine   Fev 05 Utilisation optionnelle des fichiers remotes
!
!LANGUAGE FTN77
!
!MODULES
      EXTERNAL FSTLNK, FSTOUV, FNOM, FCLOS
      INTEGER  FSTLNK, FSTOUV, FNOM, FCLOS, I, J, K

      FSTNOL = 0
      IF(N.GT.1 .AND. ((INDEX(OPTIONS,'SEQ') + INDEX(OPTIONS,'SQI') + INDEX(OPTIONS,'FTN')) .NE. 0)) THEN
         PRINT*,'PAS PERMIS D''ENCHAINER DES FICHIERS SEQUENTIELS'
         N = 0
      ELSE
!        TEST SI LES FICHIERS EXISTENT
         DO 20 I=1,N
            IF(FNOM(LISTE(I), NOMS(I), OPTIONS//'+REMOTE', 0) .NE. 0) THEN
               DO 10 J=1,I
                  K = FCLOS ( LISTE(J) )
   10             CONTINUE
               N = 0
               GO TO 30
            ENDIF
   20       CONTINUE
      ENDIF

!     LES FICHIERS EXISTENT
   30 IF(N .GT. 0) THEN
         DO 40 I=1,N
            FSTNOL = FSTNOL + FSTOUV(LISTE(I), OPTIONS)
   40       CONTINUE
         IF(N .GT. 1) I = FSTLNK(LISTE, N)
      ENDIF

      RETURN
      END 




