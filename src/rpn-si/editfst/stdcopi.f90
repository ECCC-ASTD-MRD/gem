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
!** S/P STDCOPI - APPELEE VIA DIRECTIVE STDCOPI(S,TYPE,D,TYPE)
!                 OUVRE LES FICHIERS AU BESOIN, APPEL A COPYSTX
      SUBROUTINE STDCOPI(INPT, TS, OUPT, TD)
      use configuration
      IMPLICIT NONE 
  
      INTEGER    INPT(*), TS(*), OUPT(*), TD(*)
!
!ARGUMENTS
!  ENTRE    - INPT  - DN  FICHIER SOURCE
!    "      - TS    - TYPE   "       "      (SQI,RND,FTN)
!    "      - OUPT  - DN     "    DESTINATION
!    "      - TD    - TYPE   "         "    (SQI,RND,FTN)
  
!AUTEUR       Y. BOURASSA  - AVR 86 VERSION ORIGINALE         
!REVISION 001 "     "      - FEV 92 ENLEVE EXTRENALS INUTILES 
!         002 "     "      - MAR 92 APPEL A LOW2UP POUT 'TS' ET 'DS'
!         003 "     "      - AVR 92 FIX BUG DECODE DNOM
!         004 "     "      - MAI 92 SKIP ABORT SI INTERACTIF
!         005 M. Lepine      Nov 05 remplacement de fstabt par qqexit
!         006       "        Fev 06 ouverture des fichiers sources R/O
!
!LANGAGE  - FTN77
!
!MODULERS
      EXTERNAL      ARGDIMS, OUVRES, OUVRED, COPYSTX, DMPDES, SAUVDEZ, LOW2UP,  qqexit
      INTEGER       ARGDIMS, OUVRED, I
      CHARACTER*128 DD
 
!     FICHIER SOURCE
      IF(NP.GE.2. .AND. (TS(1).NE.-1 .AND. INPT(1).NE.-1)) THEN
         WRITE(SNOM, LIN128) (TS(I), I=1,ARGDIMS(2)) 
         CALL LOW2UP(SNOM, SNOM)
         IF(INDEX(SNOM,'FTN') .GT. 0) THEN 
            SNOM = 'STD+SEQ+FTN+OLD+R/O' 
         ELSEIF(INDEX(SNOM,'SEQ').GT.0 .or. INDEX(SNOM,'SQI').GT.0) THEN
            SNOM = 'STD+SEQ+OLD+R/O'
         ELSE
            SNOM = 'STD+RND+OLD+R/O'
         ENDIF
      ENDIF
      IF(INPT(1) .NE. -1) THEN  ! on a fourni un nom
         NFS = 1
         WRITE(DD, LIN128) (INPT(I), I=1,ARGDIMS(1))
         CALL OUVRES( DD )
      ENDIF
      IF( .NOT. OUVS ) THEN
         PRINT*,'****  FICHIER SOURCE INCONNU  ****'
         CALL qqexit(80)
      ENDIF

!     FICHIER DESTINATION
      IF(NP.EQ.4 .AND. (TD(1).NE.-1 .AND. OUPT(1).NE.-1)) THEN
         WRITE(DNOM, LIN128) (TD(I), I=1,ARGDIMS(4)) 
         CALL LOW2UP(DNOM, DNOM)
         IF(INDEX(DNOM,'FTN') .GT. 0) THEN
            DNOM = 'STD+SEQ+FTN'
         ELSEIF(INDEX(DNOM,'SEQ') .GT. 0) THEN
            DNOM = 'STD+SEQ'
         ELSEIF(INDEX(DNOM,'SQI') .GT. 0) THEN
            DNOM = 'STD+SEQ'
         ELSE
            DNOM = 'STD+RND'
         ENDIF
      ENDIF
      IF(NP.GE.3 .AND. (OUPT(1) .NE. -1)) THEN
         WRITE(DD,LIN128) (OUPT(I),I=1,ARGDIMS(3))
         I = OUVRED( DD )
      ENDIF
      IF( .NOT. OUVD) THEN
         PRINT*,'****  FICHIER DESTINATION INCONNU  ****'
         CALL qqexit(81)
      ENDIF

!     DEVONS=NOUS METTRE COPYSTX,EN MODE XPRES?
      CALL DMPDES
      IF(LIMITE .NE. 0) THEN  ! nombre  limite d'enregistrements a copier pas encore atteint
         CALL COPYSTX
      ELSE
         WRITE(6,*)'LA LIMITE DES TRANSFERS DEJA ATEINTE'
      ENDIF

!     CONTROLE DE LA PORTEE DES DIRECTIVES
      CALL SAUVDEZ
   
      RETURN
      END 


