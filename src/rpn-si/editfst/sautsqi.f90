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
!** S/R SAUTSQI
!        POSITIONNE LE FICHIER SEQUENTIEL DN 
!        A UN CERTAIN NIVEAU D'EOF LOGIQUE
      SUBROUTINE SAUTSQI(DN, LEV, NL)
      use configuration
  
      IMPLICIT   NONE
        
      INTEGER    DN(*), LEV, NL
  
!ARGUMENTS
!  ENTRE    - DN    - DATASET NAME DU FICHIER DESTINATION A POSITIONNER
!    "      - LEVEL - NIVEAU DE EOF LIGIQUE RECHERCHE [1]
!    "      - NL    - NOMBRE DE DE CES EOF RECHERCHES 1]
!
!AUTEURS
!VERSION ORIGINALE  C. THIBEAULT, C. CHOUINARD ET M.VALIN (STDSAUT)
!REVISION 001       Y. BOURASSA MARS 86
!         002       "      "    NOV 90 VERSION FTN/SQI
!         003       "      "    FEV 92 LE FICHIER DOIT EXISTER
!         004       "      "    MAR 92 CORRIGE BUG DANS TEST FSTINF
!         005       "      "    MAI 92 SKIP ABORT SI EN INTERACTIF
!         006       M. Lepine   Nov 05 Remplacement de fstabt par qqexit
!         007       M. Valin    Mai 14 Remplacement des comdecks par un module
!
!LANGUAGE   - FTN77 
!
!MODULES  
      EXTERNAL      ARGDIMS, FSTINF, FSTEOF, OUVRED, qqexit 
!
!*
      INTEGER       ARGDIMS, FSTINF, FSTEOF, OUVRED, LEVEL
      INTEGER       I, J, K, L, M
      CHARACTER*128 CLE
      CHARACTER*15  SORTE, TAPE
      CHARACTER*3   DND

      DND   = 'SEQ'
      SORTE = 'STD+SEQ+OLD'
!     ASSIGNATION DES PARAMETRES PAR DEFAUT ET OUVRE FICHIER
    1 IF(NP .GT. 1) THEN
         LEVEL = LEV
      ELSE
         LEVEL = 1
      ENDIF
      IF(NP .EQ. 3) THEN
         I = NL
      ELSE
         I = 1
      ENDIF
  
      WRITE(CLE, LIN128) (DN(M), M=1,ARGDIMS(1))
      IF(CLE.EQ.ND .AND. OUVD) THEN
         IF(INDEX(DNOM, DND) .EQ. 0) THEN
            WRITE(6,*)'PROBLEME AVEC FICHIER D= ', ND
            WRITE(6,*)'DEJA OUVERT AVEC   TYPE= ', DNOM
            WRITE(6,*)'TYPE PRESEMT=            ', SORTE
            IF( INTERAC ) THEN
               RETURN
            ELSE
               CALL qqexit(61)
            ENDIF
         ENDIF
         J = 3
      ELSEIF(CLE.EQ.NS .AND. OUVS) THEN 
         IF(INDEX(SNOM, DND) .EQ. 0) THEN
            WRITE(6,*)'PROBLEME AVEC FICHIER S= ', NS
            WRITE(6,*)'DEJA OUVERT AVEC   TYPE= ', SNOM
            WRITE(6,*)'TYPE PRESEMT=            ', SORTE
            IF( INTERAC ) THEN
               RETURN
            ELSE
               CALL qqexit(62)
            ENDIF
         ENDIF
         J = SOURCES(1)
      ELSE
         DNOM = SORTE
         K    = OUVRED( CLE )
         IF(K .EQ. 0) THEN
            J    = 3
         ELSE
            PRINT*,'FICNIER N''EXISTE PAS'
            IF( INTERAC ) THEN
               RETURN
            ELSE
               CALL qqexit(63)
            ENDIF
         ENDIF
      ENDIF
      IF( DIAG .OR. DEBUG ) THEN
         TAPE = CLE
         WRITE(6,600) I, LEVEL, J, TAPE
  600    FORMAT(' SAUTE',I3,' EOF 'I2,' FICHIER',I3,'=',A15,'...')
      ENDIF       
!     SAUTE AU N..IEME EOF DE NIVEAU LEVEL
   10 IF(FSTINF(J, K, L, M, 0, '0', 0, 0, 0, '0', '0') .GE. 0) GOTO 10
      M = FSTEOF(J) 
      IF(M.LT.1 .OR. M.GT.15) THEN
         PRINT*,' MAUVAISE MARQUE DE FIN DE FICHIER =',M,' RENCONTREE DANS TAPE=',J
         CALL qqexit(64)
      ENDIF
      IF(DIAG .OR. DEBUG) WRITE(6,*)'RENCONTRE UN EOF #',M
      IF(M .LT. LEVEL) GO TO 10
      IF(M .EQ. LEVEL) THEN
         I = I-1
         IF(I .NE. 0) GO TO 10
      ENDIF
      IF(J .EQ. SOURCES(1)) LEOF = M
      RETURN
  
!**   SAUTSEQ OUVRE LE FICHIER 2 (SEQUENTIEL 'SEQ+FTN')
      ENTRY SAUTSEQ(DN, LEV, NL)
      DND   = 'FTN'
      SORTE = 'STD+SEQ+FTN+OLD'
      GO TO 1
  
      END 
