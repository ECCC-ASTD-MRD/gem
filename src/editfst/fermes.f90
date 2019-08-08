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
!** S/R FERMES FERME  LES FICHIERS SOURCES
      SUBROUTINE FERMES
      use configuration
      IMPLICIT      NONE
!
!AUTEUR
!VERSION ORIGINALE  - Y. BOURASSA NOV 90
!REVISION 001         "      "    MAR 92 VARIABLE NFSO (NOMBRE DE SOURCE OUVERTS)
!                                        CHANGE ALLEL A FATVOI
!         002         "      "    MAI 92 FCLOS SUB.>FUNCTION.
!         003         M.Valin     FEV 14 mode DRYRUN / remplacement des comdecks par un module
!     LANGUAGE FTN77
!
!
!*
      INTEGER, external :: FSTVOI, FSTFRM, FSTRWD, FSTUNL, FSTOPC, FCLOS
      integer :: I, J
  
!     TRAITEMENT DES FICHIERS SOURCES
      IF( OUVS ) THEN
         IF(NFSO .GT. 1) I = FSTUNL( )  ! unlink
         DO 10 J=1,NFSO
            I = FSTFRM( SOURCES(J) )
            I = FCLOS(  SOURCES(J) )
   10       CONTINUE
         OUVS = .FALSE.
         NFSO = 0
      ENDIF
      RETURN
      END
  
!** S/R FERMED FERME  LE FICHIER DESTINATION
      SUBROUTINE FERMED
      use configuration
      IMPLICIT      NONE
!AUTEUR
!VERSION ORIGINALE  - Y. BOURASSA NOV 90
!REVISION 001         "      "    MAR 92 VARIABLE NFSO (NOMBRE DE SOURCE OUVERTS)
!                                        CHANGE ALLEL A FATVOI
!         002         "      "    MAI 92 FCLOS SUB.>FUNCTION.
!         003         M.Valin     FEV 14 mode DRYRUN / remplacement des comdecks par un module
!*
      INTEGER, external :: FSTVOI, FSTFRM, FSTRWD, FSTUNL, FSTOPC, FCLOS
      integer :: I
!     TRAITEMENT DU FICHIER DESTINATION 
      if(dryrun) then  ! dry run, on ne fait rien
        OUVD = .FALSE.
        return
      endif
      IF( OUVD ) THEN    ! fichier destination ouvert
         IF( VD ) THEN   ! voir contenu du fichier destination
            I = FSTOPC('MSGLVL', 'INFORM', .FALSE.)
            IF( DSEQ ) I = FSTRWD(3)  ! fichier sequentiel, on rembobine
            IF( INDEX(DNOM,'FTN') .NE. 0) THEN
               I = FSTVOI(3, 'SEQ')
            ELSE
               I = FSTVOI(3, 'STD')
            ENDIF
            I = FSTOPC('MSGLVL', DEF1b(13), .FALSE.)
         ENDIF
         I = FSTFRM( 3 )
         I = FCLOS(  3 )
         OUVD = .FALSE.
      ENDIF
      RETURN
  
      END 
