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
!**   FUNCTION OUVRE UN FICHIER DESTINATION
      FUNCTION OUVRED( DN ) 
      use configuration
      IMPLICIT      NONE
      INTEGER       OUVRED
      CHARACTER*(*) DN
  
!ARGUMENTS
!SORTIE OUVRED - >=0 DIMENSION DU FICHIER OUVERT
!                <0  ERREUR DE FNOM PAS OUVERT
!ENTREE DN    -  NOM DU FICHIERUS DESTINATION
!
!AUTEUR -       Y. BOURASSA OCT 91 SEPARATION DE OUVRES
!REVISION  001  "     "     FEV 92 ENLEVE EXTERNAL INUTILE
!          002  "     "     MAR 92 CHANGE S/R EN FUNCTION
!                                  CHANGE APPEL A FNOM
!          003  "     "     MAR 92 DEPLACE LE RETURN QUAND DEJA OUVERT
!          004  "     "     MAI 92 SKIP ABORT EN INTERACTIF
!          005  "     "     DEC 94 SI FICHIER DESTINATION EST SEQUENTIEL
!                                  ON PLACE LE FICHIER A LA FIN AVANT COPIE
!          006  M. Lepine   Fev 05 Utilisation optionnelle des fichiers remotes
!          007  M. Lepine   Nov 05 Remplacement de tous les fstabt par qqexit
!          008  M. Valin    Fev-Avr 14 mode DRYRUN - simplification de la logique d'erreur
!
!LANGUAGE FTN77
!
!
!MODULES
      EXTERNAL FSTOUV, FNOM, FSTNBR, FSTINF, FERMED, qqexit

      INTEGER  FSTOUV, FNOM, FSTNBR, FSTINF, RENDUA, I, J, K
      integer ier

      if(dryrun) then  ! mode dryrun, on fait comme si on avait ouvert le fichier de sortie
        PRINT*,'DRYRUN: le fichier de sortie ne sera pas ouvert'
        ouvred = 0
        ouvd = .true.
        return
      endif
!     SI DEJA OUVERT COMME DESTINATION, TROUVER LE NOMBRE D'ENREGISTREMENTS
      IF(DN.EQ.ND .AND. OUVD) THEN
         IF(INDEX(DNOM,'SEQ') .GT. 0) THEN
            OUVRED = 0
            IF( DEBUG ) PRINT*,'FICHIER SEQUENTIEL ',ND,' DEJA OUVERT'
         ELSE
            OUVRED = FSTNBR( 3 )
            IF( DEBUG ) PRINT*,'FICHIER ',ND,' DEJA OUVERT RANDOM ','TAILLE =',OUVRED
         ENDIF
         RETURN
      ENDIF
  
!     SI DEJA OUVERT COMME SOURCE, ERREUR FATALE
      IF(DN.EQ.NS .AND. OUVS) THEN
         PRINT*,'  **************************************'
         PRINT*,' *              ATTENTION               *'
         PRINT*,'* DESTINATION = SOURCE,  ERREUR FATALE   *'
         PRINT*,' *                ABORT                 *'
         PRINT*,'  **************************************'
         CALL qqexit(54)
      ENDIF

!     FERMER LE FICHIER DESTINATION SI OUVERT
      IF( OUVD ) CALL FERMED

!     RETOURNE OUVRED >= 0 SI OUVERT
      ier = FNOM(3, DN, DNOM//'R/W+REMOTE', 0)
      IF(ier .EQ. 0) THEN
         OUVRED = FSTOUV(3, DNOM)
         if (ouvred .lt. 0) then
            PRINT*,'  **************************************'
            PRINT*,' *              ATTENTION               *'
            PRINT*,'* ERREUR DANS FSTOUV FICHIER INUTILISABLE *'
            PRINT*,' *                ABORT                 *'
            PRINT*,'  **************************************'
            CALL qqexit(55)
         endif
         ND     = DN
         OUVD   = .TRUE.
         DSEQ   = INDEX(DNOM,'SEQ') .NE. 0
         IF( DSEQ ) THEN  ! fichier sequentiel, aller se placer a la fin
            RENDUA = 0
 10         COPIES = FSTINF( 3, I, J, K, -1, ' ', -1, -1, -1, ' ', ' ')
            IF(COPIES .GE. 0) THEN
               RENDUA = RENDUA+1
               GO TO 10
            ENDIF
            IF( DEBUG ) PRINT*,'ENREGISTREMENTS AVANT COPIE =',RENDUA
         ENDIF
         COPIES = 0
      ELSE
         PRINT*,'  **************************************'
         PRINT*,' *              ATTENTION               *'
         PRINT*,'* ERREUR DANS FNOM FICHIER INUTILISABLE  *'
         PRINT*,' *                ABORT                 *'
         PRINT*,'  **************************************'
         CALL qqexit(55)
      ENDIF

      RETURN
      END 
