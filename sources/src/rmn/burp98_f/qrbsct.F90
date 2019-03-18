!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
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
!DECK QRBSCT
!.S QRBSCT
!**S/P  QRBSCT - INITIALISER LE TABLEAU DE CONVERSIONS
      FUNCTION QRBSCT( TABLEAU,           TABDIM, NELELU)
      IMPLICIT NONE
      INTEGER TABDIM
      INTEGER  QRBSCT, TABLEAU(3,TABDIM), NELELU 

!AUTEUR: J. CAVEEN   FEVRIER 1991
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 j. caveen   sept.   1995 changer tableburp pour table_b_bufr
!REV 003 j. caveen   Avril   1996 ajout du traitement des commentaires '#'
!REV 004 M. Lepine   Nov     1996 changement d'un format de lecture
!REV 005 M. Lepine   Mai     1997 changement d'un autre format de lecture
!
!OBJET( QRBSCT )
!     SOUS-PROGRAMME SERVANT A INITIALISER LE TABLEAU CONTENANT
!     LES VALEURS QUI SERVIRONT A FAIRE LA CONVERSION D'UNITES
!     CMC A BUFR OU L'INVERSE.
!     POUR UN ELEMENT I, 
!     LE NOM DE L'ELEMENT (CODE) SE RETROUVE DANS TABLEAU(1,I)
!     L'EXPOSANT DU FACTEUR D'ECHELLE EST DANS    TABLEAU(2,I)
!     LA VALEUR DE REFERENCE SE RETROUVE DANS     TABLEAU(3,I)
!
!ARGUMENTS
!     TABLEAU  SORTIE  TABLEAU CONTENANT LES VALEURS POUR 
!                      LES CONVERSIONS D'UNITES
!     TABDIM   ENTREE  DEUXIEME DIMENSION DE TABLEAU
!     NELELU   SORTI   NOMBRE D'ELEMENTS LUS
!
!IMPLICITES
#include "codes.cdk"
#include "defi.cdk"
#include "burpopt.cdk"
#include <ftnmacros.hf>
!
!MODULES
      EXTERNAL LONGUEUR, QDFERR, MRBCOV, FNOM, FCLOS
      INTEGER  LONGUEUR, QDFERR, MRBCOV, FNOM, MAXELM, I, IUN, PATHLNG,  OUTFIL, TRAVAIL
      CHARACTER*128 LIGNE, PATH, PATH1 
      DATA OUTFIL   /6/
      DATA RPETITIF / MAXREP*0/
!
!*

!     ON OBTIENT LE NOM DE LA TABLEBURP A OUVRIR
!     SI LA VARIABLE TABLEBURP EST DEFINIE, ON OUVRE LA TABLE 
!     EXPERIMENTALE ET ON EMET UN MESSAGE AVERTISSANT
!     QU'ON N'UTILISE PAS LE FICHIER tableburp OFFICIEL.
      CALL GETENV('MA_TABLEBURP_PERSONNELLE',PATH1)
      PATHLNG = LONGUEUR( PATH1 )
      IF(PATHLNG .GT. 0) THEN
         WRITE(OUTFIL, 700)
         WRITE(OUTFIL, 703)
         WRITE(OUTFIL, 701)
         WRITE(OUTFIL, 703)
         WRITE(OUTFIL, 702)

!        ON MET LE MARQUEUR D'UTILISATION DE LA MAUVAISE
!        TABLEBURP A 1
         BADTBL = 1
      ELSE
         CALL GETENV('AFSISIO', PATH)
         PATHLNG = LONGUEUR( PATH )
         PATH1   = PATH(1:PATHLNG)//'/datafiles/constants/table_b_bufr'
!         PATH1   = PATH(1:PATHLNG)//'/datafiles/constants/tableburp'
      ENDIF

      QRBSCT = -1
      IUN    = 0
      QRBSCT = FNOM(IUN, PATH1, 'FTN+FMT+R/O', 0)
      IF(QRBSCT .NE. 0) THEN
         QRBSCT = QDFERR('QRBSCT', ' ERREUR D''OUVERTURE DU FICHIER TABLEBURP', KAPUT, ERBTAB)
         RETURN
      ENDIF

!     LIRE LE NOMBRE D'ELEMENTS PRESENTS ET LE NOMBRE D'ELEMENTS
!     A CONVERTIR DANS LE FICHIER BUFR_B
      READ(IUN, *) MAXELM, NELELU
        
!      IF(NELELU .GT. TABDIM) THEN
!         QRBSCT = QDFERR('QRBSCT',
!     $ ' TABLEAU POUR LA LECTURE EST TROP PETIT CONSULTER SPECIALISTE',
!     $        ERROR, ERBTAB)
!         RETURN
!      ENDIF

!     LECTURE DE TOUTES LES ENTREES DE LA TABLE BUFRB
      I = 0
 100  READ(IUN,'(A128)', END=300) LIGNE
      IF((LIGNE(1:1) .EQ. '*') .or. (ligne(1:1) .eq. '#')) GO TO 100

!     CONVERSION DES NOM DE VARIABLES A DES ENTIERS DE 16 BITS
      READ(LIGNE(1:6),'(I6)') TRAVAIL
      TRAVAIL =  MRBCOV( TRAVAIL )

!     VERIFIER SI L'ELEMENT EST REPETITIF. SI OUI, ON ALLUME
!     LE BIT CORRESPONDANT AU NO D'ELEMENT DANS RPETITIF
      IF(INDEX('Mm',LIGNE(85:85)) .NE. 0 ) PUTBIT(RPETITIF, 1, TRAVAIL, 1)
      
        
      IF(LIGNE(51:51) .NE. '*') THEN
         I = I+1
         IF (I .GT. TABDIM) THEN
           QRBSCT = QDFERR('QRBSCT', ' TABLEAU POUR LA LECTURE EST TROP PETIT CONSULTER SPECIALISTE', ERROR, ERBTAB)
           RETURN
         ENDIF
         TABLEAU(1,I) = TRAVAIL
         READ(LIGNE(64:66),'(I3)') TABLEAU(2,I)
         READ(LIGNE(67:77),'(I11)') TABLEAU(3,I)
      ENDIF

      GO TO 100
 300  CONTINUE

      IF(NELELU .NE. I) THEN
         QRBSCT = QDFERR('QRBSCT', 'AVERTISSEMENT- NELELU <> NOMBRE D ENTREE DANS TABLEBURP', INFORM, ERELEM)
         NELELU = I
      ENDIF

      CALL FCLOS( IUN )
      QRBSCT = 0
      RETURN


! 600  FORMAT(I3, 1X, I3)
 700  FORMAT('0***********************ATTENTION***********************')
 701  FORMAT(' *   ON N''UTILISE PAS LE FICHIER TABLEBURP OFFICIEL   *')
 702  FORMAT(' *******************************************************')
 703  FORMAT(' *                                                     *')
      END
