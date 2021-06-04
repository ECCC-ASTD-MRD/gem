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
!.S MRBCVT
!**S/PMRBCVT - FAIRE UNE CONVERSION D'UNITES
      FUNCTION MRBCVT( LISTE,  TBLVAL, RVAL, NELE, NVAL, NT, MODE)
      IMPLICIT NONE
      INTEGER  MRBCVT, MRBSCT, MRBTBL, NELE, NVAL, NT, MODE,  &
     &         LISTE(NELE), TBLVAL(NELE, NVAL, NT), ELEUSR, NSLOTS,  &
     &         NELEUSR, TBLUSR(3, NELEUSR), TBLBURP(NSLOTS, NELE)
      REAL     RVAL(NELE, NVAL, NT)
!
!AUTEUR: J. CAVEEN   JANVIER 1991
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( MRBCVT )
!     SOUS-PROGRAMME SERVANT A LA CONVERSION D'UNITES DE REEL A ENTIER
!     POSITIF OU L'INVERSE, SELON LA VALEUR DE MODE.  LA CONVERSION
!     SE FAIT EN CONSULTANT UN TABLEAU DUQUEL ON EXTRAIT UN FACTEUR
!     D'ECHELLE ET UNE VALEUR DE REFERENCE QUI SONT APPLIQUES A LA 
!     VARIABLE A CONVERTIR.
!     LE TABLEAU DE REFERENCE EST INITIALISE PAR LE SOUS-PROGRAMME
!     MRBSCT QUI EST APPELE LORS DE LA PREMIERE EXECUTION DE MRBCVT
!
!ARGUMENTS
!     LISTE    ENTREE    LISTE DES ELEMENTS A CONVERTIR
!     NELE       "       NOMBRE D'ELEMENTS (LONGUEUR DE LISTE)
!     NVAL       "       NOMBRE DE VALEURS PAR ELEMENT
!     NT         "       NOMBRE D'ENSEMBLES NELE*NVAL
!     MODE       "       MODE DE CONVERSION
!                        MODE = 0 DE TBLVAL(CODE BUFR) A RVAL
!                        MODE = 1 DE RVAL A TBLVAL(CODE BUFR)
!     TBLVAL   ENT/SRT   TABLEAU DE VALEURS EN CODE BUFR
!     RVAL               TABLEAU DE VALEURS REELLES
!
!IMPLICITES 
#include "codes.cdk"
#include "defi.cdk"
#include "burpopt.cdk"
#include <ftnmacros.hf>
!
!MODULES
      EXTERNAL QRBSCT, QDFERR, BUFRCHR, QBRPTRI 
      INTEGER  QRBSCT, QDFERR, BUFRCHR

      INTEGER TABLEAU(3,MAXNELE), NELELU
!
!*

!******************************************************************************
!       LA MATRICE TABLEAU CONTIENT:
!          TABLEAU(1,I) - CODE DE L'ELEMENT I
!          TABLEAU(2,I) - FACTEUR A APPLIQUER A L'ELEMENT I
!          TABLEAU(3,I) - VALEUR DE REFERENCE A AJOUTER A L'ELEMENT I
!       LA VARIABLE NELELU INDIQUE LE NOMBRE D'ELEMENTS PRESENT DANS LE
!       FICHIER BUFR
!
!       POUR CODER LA VALEUR D'UN ELEMENT AVANT UN APPEL A MRFPUT, ON FAIT 
!       L'OPERATION SUIVANTE: 
!       ELEMENT(I)_CODE = ELEMENT(I) * TABLEAU(2,I) - TABLEAU(3,I)
!
!       ON NE FAIT AUCUNE CONVERSION LORSQUE QU'UN ELEMENT EST DEJA CODE
!       COMME PAR EXEMPLE POUR LES DIFFERENTS MARQUEURS.
!
!       POUR DECODER UN ELEMENT ON FAIT L'OPERATION INVERSE.  DANS LE CAS 
!       DES ELEMENTS NE REQUERANT AUCUN DECODAGE (E.G. MARQUEURS), ON INSERE
!       DANS LE TABLEAU RVAL LA VALEUR -1.1E30 CE QUI INDIQUE A L'USAGER 
!       QU'IL DOIT CONSULTER LE TABLEAU TBLVAL POUR OBTENIR CET ELEMEMT
!
!
!*****************************************************************************
      LOGICAL PREMIER
      DATA    PREMIER /.TRUE./
      SAVE    TABLEAU, NELELU, PREMIER 
      INTEGER I, J, K, L, REFEREN, ZEROCPL
      REAL    ECHELE

      MRBCVT = -1
      IF( PREMIER ) THEN
          MRBCVT = QRBSCT(TABLEAU,MAXNELE,NELELU)
          IF(MRBCVT .EQ. -ERBTAB) RETURN

!         TRIER LE TABLEAU PAR ELEMENT EN ORDRE CROISSANT
!         (NECESSAIRE POUR LA RECHERCHE BINAIRE QUE L'ON FAIT PLUS TARD)
          CALL QBRPTRI(TABLEAU, 3, NELELU)
          PREMIER = .FALSE.
      ENDIF

      ZEROCPL = COMPL( 0 )

      DO 50 I = 1, NELE
!        TROUVER L'INDEX J POINTANT A L'ELEMENT DANS TABLEAU
         ECHELE  = 1.0
         REFEREN = 0
         J       = BUFRCHR(LISTE(I), TABLEAU, NELELU)
         IF(J .GT. 0) THEN
            ECHELE  = 10.0 ** TABLEAU(2,J)
            REFEREN = TABLEAU(3,J)
!            print *,'Debug+ ECHELE=',ECHELE,' TABLEAU(2,J)=',TABLEAU(2,J)
!            print *,'Debug+ REFEREN=',REFEREN

!           CONVERSION DE CODE BUFR A UNITES CMC
            IF(MODE .EQ. 0) THEN
               DO 20 L = 1, NT
                  DO 10 K = 1, NVAL
                     IF(TBLVAL(I,K,L) .EQ. ZEROCPL) THEN
                        RVAL(I,K,L) = MANQUE
                     ELSE
                        IF(TBLVAL(I,K,L) .LT. 0) TBLVAL(I,K,L) = TBLVAL(I,K,L) + 1
                        RVAL(I,K,L) = FLOAT(TBLVAL(I,K,L) + REFEREN)/ECHELE
                     ENDIF
10                   CONTINUE
20                CONTINUE 
            ELSE
!              CODAGE BUFR
               DO 40 L = 1, NT
                  DO 30 K = 1, NVAL
                     IF(RVAL(I,K,L) .EQ. MANQUE) THEN
                        TBLVAL(I,K,L) = ZEROCPL
                     ELSE
                        TBLVAL(I,K,L) = NINT(RVAL(I,K,L) * ECHELE) - REFEREN
                        IF(TBLVAL(I,K,L) .LT. 0) TBLVAL(I,K,L) = TBLVAL(I,K,L) - 1
                     ENDIF
30                  CONTINUE
40               CONTINUE
            ENDIF
         ENDIF
50       CONTINUE

      MRBCVT = 0
      RETURN

!**S/P - MRBSCT - INITIALISER LE TABLEAU DE CONVERSION DE L'USAGER
      ENTRY MRBSCT(TBLUSR, NELEUSR)
!
!OBJET( MRBSCT ) 
!     SOUS-PROGRAMME SERVANT A AJOUTER AU TABLEAU DE
!     CONVERSION STANDARD, UNE LISTE D'ELEMENTS QUE
!     L'USAGER A DEFINIE LUI-MEME.
!     SI LE TABLEAU STANDARD N'A PAS ETE INITIALISE,
!     ON APPELLE QRBSCT POUR EN FAIRE L'INITIALISATION.
!     ON AJOUTE LE TABLEAU DE L'USAGER A LA FIN.
!
      MRBSCT = -1
      IF( PREMIER ) THEN
          MRBSCT = QRBSCT(TABLEAU,MAXNELE,NELELU)
          IF (MRBSCT .EQ. -ERBTAB) RETURN

!         TRIER LE TABLEAU PAR ELEMENT EN ORDRE CROISSANT
!         (NECESSAIRE POUR LA RECHERCHE BINAIRE QUE L'ON FAIT PLUS TARD)
          CALL QBRPTRI(TABLEAU, 3, NELELU)
          PREMIER = .FALSE.
      ENDIF

!     AJOUTER LE TABLEAU DE L'USAGER A LA FIN DU TABLEAU STANDARD
!     APRES L'AVOIR TRIE
      CALL QBRPTRI(TBLUSR, 3, NELEUSR)

      DO 70 J = 1, NELEUSR
         NELELU = NELELU + 1
         DO 60 I = 1, 3
             TABLEAU(I,NELELU) = TBLUSR(I,J)
60           CONTINUE
70       CONTINUE
      MRBSCT = 0

      RETURN

!**S/P - MRBTBL - REMPLIR UN TABLEAU A PARTIR DE TABLEBURP
      ENTRY MRBTBL(TBLBURP, NSLOTS, NELE)
!
!OBJET( MRBTBL )
!     SOUS-PROGRAMME SERVANT A REMPLIR LE TABLEAU TBLBURP
!     A PARTIR DES DESCRIPTIONS D'ELEMENTS TROUVEES DANS
!     LE FICHIER TABLEBURP.  POUR CHAQUE ELEMENT, 
!     ON RETOURNE:
!        - FACTEUR D'ECHELLE
!        - VALEUR DE REFERENCE
!        - SI L'ELEMENT EST CONVERTISSABLE OU NON
!          0 - non convertissable
!          1 - convertissable
!
!ARGUMENTS:
!     NELE      ENTREE        - NOMBRE D'ELEMENTS A TRAITER
!     TBLBURP   ENTREE CONTIENT LES CODES D'ELEMENTS
!        "      SORTIE CONTIENT LES PARAMETRES DE CHAQUE ELEMENT
!
!     ARANGEMENT DE TBLBURP:
!
!     ----------------------------------------------------------
!     | code elem 16 bits | echelle | reference | convertissable |
!     |                   |         |           |                |
!              .               .          .             .
!              .               .          .             .
!              .               .          .             .
!
!*


      MRBTBL = -1

!     S'ASSURER QUE LA VALEUR DE NSLOTS EST BONNE 
      IF(NSLOTS .NE. NCELLMAX) THEN
         MRBTBL = QDFERR('MRBTBL', 'DIMENSION NCELL INCORRECTE', WARNIN, ERRCELL)
         RETURN
      ENDIF

      IF( PREMIER ) THEN
          MRBTBL = QRBSCT(TABLEAU,MAXNELE,NELELU)
          IF (MRBTBL .EQ. -ERBTAB) RETURN

!         TRIER LE TABLEAU PAR ELEMENT EN ORDRE CROISSANT
!         (NECESSAIRE POUR LA RECHERCHE BINAIRE QUE L'ON FAIT PLUS TARD)
          CALL QBRPTRI(TABLEAU, 3, NELELU)
          PREMIER = .FALSE.
      ENDIF

!     TRAITER CHAQUE ELEMENT DU TABLEAU TBLBURP
      DO 80 I = 1, NELE
!        TROUVER L'INDEX J POINTANT A L'ELEMENT DANS TABLEAU
         J = BUFRCHR(TBLBURP(1,I), TABLEAU, NELELU)
         IF(J .GT. 0) THEN
            TBLBURP(2,I) = TABLEAU(2,J)
            TBLBURP(3,I) = TABLEAU(3,J)
            TBLBURP(4,I) = 1
         ELSE
            TBLBURP(4,I) = 0
         ENDIF
80       CONTINUE

      MRBTBL = 0

      RETURN
      END
