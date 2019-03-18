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
!** S/R WEOFILE
!     ECRIT UNE MARQUE DE FIN DE FICHIER LOGIQUE DE NIVO M DANS LE FICHIER
!     DESTINATION LA POSITION COURANTE. 
!
      SUBROUTINE WEOFILE(OUPT, MARC, TD)
      use configuration
      IMPLICIT NONE 
  
      INTEGER    OUPT(*), MARC, TD(*)
!
!ARGUMENTS
!  ENTRE    - OUPT    - DN FICHIER DESTINATION
!    "      - MARC    - MARQUE A ECRIRE DANS FICHIER DESTINATION
!    "      - TD      - TYPE FICHIER DESTINATION (SEQ,STD/SEQ,STD/SQI)
!
!LANGUAGE   - FTN77 
!
!AUTEURS
!
!AUTEURS
!VERSION ORIGINALE Y. BOURASSA - AVR 86
!REVISION 001      "     "       OCT 90 ADAPTATION AUX FSTDS.
!         002      "     "       MAR 92 APPEL A LOW2UP POUR TD
!         003      "     "       AVR 93 FIX BUG DECODE DNOM
!         004      "     "       MAI 92 SKIP ABORT SI INTERACTIF
!         005      M. Lepine     Nov 05 remplacement de fstabt par qqexit
!         006      M. Valin      Mai 14 Remplacement des comdecks par un module
!
!MODULES
      EXTERNAL      ARGDIMS, FSTWEO, OUVRED, qqexit, LOW2UP
!
!*
      INTEGER       ARGDIMS, FSTWEO, OUVRED, L, M
      CHARACTER*128 DN

      M = 1
      GO TO(30, 20, 10) NP

!     DECODE TYPE DU FICHIER DESTINATION
   10 IF(TD(1).NE.-1 .AND. OUPT(1).NE.-1) THEN
         WRITE(DNOM, LIN128) (TD(L), L=1,ARGDIMS(3))
         CALL LOW2UP(DNOM, DNOM)
         IF(INDEX(DNOM,'FTN') .GT. 0) THEN
            DNOM = 'STD+SEQ+FTN'
         ELSEIF(INDEX(DNOM,'SEQ').GT.0 .OR. INDEX(DNOM,'SQI').GT.0) THEN
            DNOM = 'STD+SEQ'
         ENDIF
      ENDIF

   20 IF(INDEX(DNOM,'SEQ') .EQ. 0) THEN
         PRINT*,'IMPOSSIBLR DE MARQUER IN FICHIER RND'
         IF( INTERAC ) THEN
            RETURN
         ELSE
            CALL qqexit(90)
         ENDIF
      ENDIF

      IF(NP.GT.1) THEN
         IF(MARC.GT.14 .OR. MARC.LT.1) THEN
            WRITE(6,*)'MARQUE DE FIN DE FICHIER LOGIQUE ',MARC,' INACCEPTABLE, DOIT ETRE >0 ET <15.'
            IF( INTERAC ) THEN
               RETURN
            ELSE
               CALL qqexit(91)
            ENDIF
         ELSE
            M = MARC
         ENDIF
      ENDIF
     
   30 IF(OUPT(1) .NE. -1) THEN
         WRITE(DN,LIN128) (OUPT(L),L=1,ARGDIMS(1))
!        OUVERTURE LE FICHIER DESTINATION
         L = OUVRED( DN )
      ENDIF
  
      IF(.NOT. OUVD) THEN
         PRINT*,'IMPOSSIBLR DE MARQUER FICHIER INCONNU'
         IF( INTERAC ) THEN
            RETURN
         ELSE
            CALL qqexit(92)
         ENDIF
      ENDIF

!     AJOUTE LA MARQUE LOGIQUE
      L = FSTWEO(3, M)
      IF( DIAG ) WRITE(6,*)' MARQUE ',M,' AJOUTEE AU FICHIER ',ND
  
      RETURN
      END 

