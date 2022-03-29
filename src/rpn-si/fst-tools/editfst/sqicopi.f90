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
!** S/R SQICOPI COPIE UN FICHIER SEQUENTIEL SQI DANS UN FICHIER STANDARD
      SUBROUTINE SQICOPI(INPT, OUPT, TD, PR, TO, PS, NB)
      use configuration
      IMPLICIT   NONE
      INTEGER    INPT(*), OUPT(*), TD(*), PR, TO, PS, NB
!
!ARGUMENTS
!  ENTRE    - INPT  - DN  FICHIER SOURCE
!    "      - OUPT  -  "     "    DESTINATION
!    "      - TD    - TYPE   "         "      (SEQ,RND,SEQ+FTN)
!    "      - PR    - EOF LIGIQUE PRECEDANT ENREGISTREMENT RECHERCHE [0]
!    "      - TO    - EOF LIGIQUE TERMINANT ENREGISTREMENT RECHERCHE [0]
!    "      - PS    - EOF LIGIQUE QUI TERMINE LE CYCLE [0]
!    "      - NB    - NOMBRE DE BOUCLE A EXECUTER [MIN = 1] 
!
!AUTEURS
!VERSION ORIGINALE  Y. BOURASSA AVRL 86
!REVISION 001       "      "    NOV  90 ACCEPTE FICHIER SQI COMME FICHIER SOURCE
!         002       "      "    JUIL 91 ACCEPTE -1 POUR LES ARGUMENTS INPT, OUPT, TD
!         003       "      "    MARS 92 FICHIER SOURCE DOIT EXISTER CHANGE TEST FSTINF
!                                       CHANGE ALLEL A OUVRED
!         004       "      "    MAI 92  SKIP ABORT SI EN INTERACTIF
!                                       OUVRE FICHIER SOURCE AVANT DESTINATION
!         005       M. Lepine   Nov 05  rempalcement de fstabt par qqexit
!         006       M. Valin    Mai 14 Remplacement des comdecks par un module
!
!LANGUAGE   - FTN77 
!
!MODULES  
      EXTERNAL      SAUVDEZ, qqexit, COPYSTX, OUVRES, OUVRED, DMPDES
      EXTERNAL      FSTINF,  FSTEOF, ARGDIMS, FERMES, LOW2UP
      INTEGER       FSTINF,  FSTEOF, ARGDIMS, OUVRED, NI, NJ, NK
      INTEGER       N, M, I, J
      INTEGER       PRE, CSD,  POS
      CHARACTER*128 DD
  
!     INITIALISATION
      SNOM = 'STD+SEQ+OLD'
    1 PRE  = 0
      CSD  = 0
      POS  = 0
      N    = -1
      J    = SOURCES(1)
!     UTILISATION DES PARAMETRES PASSES PAR L'APPELEUR
!     DIRECTION A PRENDRE SELON LE NOMBRE D'ARGUMENTS PASSES
      GO TO(50, 50, 50, 40, 30, 20, 10) NP
   10 N    = NB
   20 POS  = PS
   30 CSD  = TO
   40 PRE  = PR
      IF(PRE .GT. MEOF) THEN
         PRINT*,' PRE = ',PRE, '    MAXEOF = ',MEOF
         PRINT*,'*    ON DEMANDE DE COPIER PASSE MEOF   *'
         PRINT*,'*        ABORT DANS SUB. SEQCOPI       *'
         IF( INTERAC ) THEN
            RETURN
         ELSE
            CALL qqexit(70)
         ENDIF
      ENDIF


!     FICHIER DESTINATION
   50 IF(NP.GE.3 .AND. (TD(1).NE.-1 .AND. OUPT(1).NE.-1)) THEN
         WRITE(DNOM, LIN128) (TD(I), I=1,ARGDIMS(3)) 
         CALL LOW2UP(DNOM, DNOM)
         IF(INDEX(DNOM,'FTN') .GT. 0) THEN
            DNOM = 'STD+SEQ+FTN' 
         ELSEIF(INDEX(DNOM,'SEQ').GT.0 .OR. INDEX(DNOM,'SQI').GT.0) THEN
            DNOM = 'STD+SEQ'
         ELSE
            DNOM = 'STD+RND'
         ENDIF
      ENDIF

!     FICHIER SOURCE
      IF(INPT(1) .NE. -1) THEN
         WRITE(DD,LIN128) (INPT(I),I=1,ARGDIMS(1))
         NFS = 1
         CALL OUVRES( DD )
      ENDIF
      IF( .NOT. OUVS ) THEN
         PRINT*,'****  FICHIER SOURCE INCONNU  ****'
         IF( INTERAC ) THEN
            RETURN
         ELSE
            CALL qqexit(71)
         ENDIF
      ENDIF

      IF(NP.GE.2 .AND. (OUPT(1) .NE. -1)) THEN
         WRITE(DD,LIN128) (OUPT(I),I=1,ARGDIMS(2))
         I = OUVRED( DD )
      ENDIF
      IF( .NOT. OUVD) THEN
         PRINT*,'****  FICHIER DESTINATION INCONNU  ****'
         IF( INTERAC ) THEN
            RETURN
         ELSE
            CALL qqexit(72)
         ENDIF
      ENDIF

!     DEVONS=NOUS METTRE COPYSTX,EN MODE XPRES?
      CALL DMPDES
      M    = MEOF
      MEOF = CSD
  
!     BOUCLE DES N REPETITIONS (N=-1 ON BOUCLE JUSQUA LA FIN)
!     SAUTE DES MARQUES DE FIN DE FICHIER LOGIQUES AVANT COPIE
   80 LEOF = 0
  
      IF(PRE .GT. 0) THEN
!        SAUTE AU PROCHAIN EOF NIVEAU PRE
   90    IF(FSTINF(J, NI, NJ, NK, 0, '0', 0, 0, 0, '0', '0') .GE. 0) GOTO 90
         LEOF = FSTEOF(J)
         IF( DIAG ) WRITE(6,*)'RENCONTRE UN EOF #',LEOF
         IF(LEOF.GT.15 .OR. LEOF.LT.1) THEN
            WRITE(6,*) LEOF,' N''EST PAS ACCEPTABLE COMME EOF LOGOQUE'
            CALL qqexit(73)
         ENDIF
         IF(LEOF .LT. PRE) GO TO 90
         IF(LEOF.GE.POS .AND. POS.NE.0) GO TO 110 
      ENDIF
  
!     SI POSSIBLE COPIE JUAQU'AU PROCHAIN EOF DESIGNE
      IF(LEOF.LT.M .AND. (LEOF.LT.CSD .OR. CSD.EQ.0)) THEN
         IF(LIMITE .NE. 0) THEN
            CALL COPYSTX
         ELSE
            WRITE(6,*)'LA LIMITE DES TRANSFERS DEJA ATEINTE'
            GO TO 120
         ENDIF
      ENDIF
  
!     SAUTE DES MARQUES DE FIN DE FICHIER LOGIQUES APRES COPIE
      IF(LEOF .LT. POS) THEN
!        SAUTE AU PROCHAIN EOF NIVEAU POS
  100    IF(FSTINF(J, NI, NJ, NK, 0, '0', 0, 0, 0, '0', '0') .GE. 0) GOTO 100
         LEOF = FSTEOF(J)
         IF( DIAG ) WRITE(6,*)'RENCONTRE UN EOF #',LEOF
         IF(LEOF.GT.15 .OR. LEOF.LT.1) THEN
            WRITE(6,*) LEOF,' N''EST PAS ACCEPTABLE COMME EOF LOGOQUE'
            CALL qqexit(74)
         ENDIF
         IF(LEOF .LT. POS) GO TO 100
      ENDIF
  
  110 IF(LEOF .LT. M) THEN
         N = N-1
         IF(N .NE. 0) GO TO 80
      ENDIF
  
!     CONTROLE DE LA PORTEE DES DIRECTIVES
  120 CALL SAUVDEZ
      MEOF = M
      RETURN
  
!**   COPIE UN FICHIER SEQUENTIEL FTN DANS UN FICHIER STANDARD
      ENTRY SEQCOPI(INPT, OUPT, TD, PR, TO, PS, NB)
      SNOM = 'STD+SEQ+FTN+OLD'
      GO TO 1
  
      END 
