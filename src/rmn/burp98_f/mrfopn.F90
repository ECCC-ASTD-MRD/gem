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
!.S MRFOPN
!**S/P MRFOPN - OUVRIR UN FICHIER BURP
!
      FUNCTION MRFOPN( IUN, INMODE)
      IMPLICIT NONE
      INTEGER  MRFOPN, IUN
      CHARACTER*(*)         INMODE
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!
!OBJET( MRFOPN )
!     INITIALISER LES DESCRIPTEURS DE CLEFS ET CREER UN FICHIER
!     BURP OU OUVRIR UN FICHIER BURP DEJA EXISTANT.
!     MRFOPN RETOURNE LE NOMBRE D'ENREGISTREMENTS ACTIFS
!     CONTENUS DANS LE FICHIER.
!                                                                       
!ARGUMENTS
!
!     IUN     ENTREE  NUMERO DU FICHIER A OUVRIR
!     INMODE    "     MODE D'OUVERTURE (READ,CREATE,APPEND)
!                                                         
!IMPLICITES
#include "defi.cdk"
#include "codes.cdk"
#include "burpopt.cdk"
#include "bpl.cdk"
#include "enforc8.cdk"
!
!MODULES 
      EXTERNAL XDFOPN, XDFCLE, XDFSTA, QDFERR, MRFOPR, MRFOPC
      external longueur,genvdt8
      integer longueur
      INTEGER  XDFOPN, XDFCLE, XDFSTA, QDFERR, MRFOPR, MRFOPC, IOUT,  &
     &         NOMBRE, STAT,   PRII,   PRI(2,NPRITOT), AUX(2,NAUXTOT),  &
     &         IER,    AUXX
      logical initdone
      integer lng
      character * 128 STRICT8
      CHARACTER*4 APPL, VERSN
      CHARACTER*6 MODE
      DATA        IOUT /6/
      DATA        BADTBL /0/
      data initdone /.false./
      SAVE initdone
!
!*

!     VERIFIER SI LE FICHIER EST OUVERT DANS UN MODE PERMIS
      if (.not. initdone) then
!         call getenv('ENFORC8',STRICT8)
!         lng = longueur(STRICT8)
!         if (lng .gt. 0) then
!            enforc8 = .true.
!         else
!            enforc8 = .false.
!         endif
         call genvdt8(enforc8)
         initdone = .true.
      endif
      MRFOPN = -1
      MODE   = INMODE
      IF(INDEX(MODE, 'WRITE').NE.0 .OR. INDEX(MODE, 'R-W').NE.0) THEN
         MRFOPN = QDFERR('MRFOPN',  &
     &        'SEULS LES MODES READ, CREATE ET APPEND SONT PERMIS',  &
     &        ERROR, ERFMOD)
         RETURN
      ENDIF

!     INITIALISER LES CLEF SI MODE = 'CREATE'
      IF(INDEX(MODE, 'CREATE') .NE. 0) THEN
         IER = 0

!        DEFINITION DES CLEFS PRIMAIRES
         IER=IER+XDFCLE('STI1', BPSTI1, LSTI1, 33, PRI(1,1),PRI(2,1))
         IER=IER+XDFCLE('STI2', BPSTI2, LSTI2, 33, PRI(1,2),PRI(2,2))
         IER=IER+XDFCLE('STI3', BPSTI3, LSTI3, 33, PRI(1,3),PRI(2,3))
         IER=IER+XDFCLE('STI4', BPSTI4, LSTI4, 33, PRI(1,4),PRI(2,4))
         IER=IER+XDFCLE('STI5', BPSTI5, LSTI5, 33, PRI(1,5),PRI(2,5))
         IER=IER+XDFCLE('STI6', BPSTI6, LSTI6, 33, PRI(1,6),PRI(2,6))
         IER=IER+XDFCLE('STI7', BPSTI7, LSTI7, 33, PRI(1,7),PRI(2,7))
         IER=IER+XDFCLE('STI8', BPSTI8, LSTI8, 33, PRI(1,8),PRI(2,8))
         IER=IER+XDFCLE('STI9', BPSTI9, LSTI9, 33, PRI(1,9),PRI(2,9))
         IER=IER+XDFCLE('FLGS', BPFLGS, LFLGS, 00, PRI(1,10),PRI(2,10))
         IER=IER+XDFCLE('LATI', BPLATI, LLATI, 00, PRI(1,11),PRI(2,11))
         IER=IER+XDFCLE('LONG', BPLONG, LLONG, 00, PRI(1,12),PRI(2,12))
         IER=IER+XDFCLE('DATE', BPDATE, LDATE, 00, PRI(1,13),PRI(2,13))
         IER=IER+XDFCLE('DX',   BPDX,   LDX,   00, PRI(1,14),PRI(2,14))
         IER=IER+XDFCLE('IDTP', BPIDTP, LIDTP, 00, PRI(1,15),PRI(2,15))
         IER=IER+XDFCLE('DY',   BPDY,   LDY,   00, PRI(1,16),PRI(2,16))
         IER=IER+XDFCLE('HEUR', BPHEUR, LHEUR, 00, PRI(1,17),PRI(2,17))
         IER=IER+XDFCLE('MIN',  BPMIN,  LMIN,  00, PRI(1,18),PRI(2,18))
         
!        DEFINITION DES CLEFS AUXILIAIRES 
         IER=IER+XDFCLE('NBLK', BPNBLK, LNBLK, 00, AUX(1,1),AUX(2,1))
         IER=IER+XDFCLE('OARS', BPOARS, LOARS, 00, AUX(1,2),AUX(2,2))
         IER=IER+XDFCLE('ELEV', BPELEV, LELEV, 00, AUX(1,3),AUX(2,3))
         IER=IER+XDFCLE('DRCV', BPDRCV, LDRCV, 00, AUX(1,4),AUX(2,4))
         IER=IER+XDFCLE('RUNN', BPRUNN, LRUNN, 00, AUX(1,5),AUX(2,5))
         IF(IER .LT. 0) RETURN
      ELSE
         IF(INDEX(MODE, 'APPEND') .NE. 0) MODE='R-W'
      ENDIF

!     OUVERTURE DU FICHIER
      NOMBRE = XDFOPN(IUN, MODE, PRI, NPRITOT, AUX, NAUXTOT, 'BRP0')
      IF(NOMBRE .LT. 0) THEN
         MRFOPN = NOMBRE
         RETURN
      ENDIF

!     OBTENIR LES INFORMATIONS CONCERNANT LE FICHIER.
      MRFOPN = XDFSTA(IUN, STAT, 0, PRII, 0, AUXX, 0, VERSN, APPL)
      IF((INDEX(VERSN,'XDF').EQ.0) .OR. ((INDEX(APPL,'BRP0').EQ.0) .AND. (INDEX(APPL,'bRp0') .EQ. 0))) THEN
         MRFOPN = QDFERR('MRFOPN', 'LE FICHIER N''EST PAS UN FICHIER RAPPORT', ERROR, ERFRAP)
         RETURN
      ENDIF

!     S'ASSURER QUE LE FICHIER A ETE CREE EN UTILISANT LA BONNE 
!     TABLE BURP
      IF(INDEX(APPL,'bRp0') .NE. 0) THEN
         IF (MESSNIV .LE. WARNIN) THEN
           WRITE(0, 700)
           WRITE(0, 703)
         ENDIF
         MRFOPN = QDFERR('MRFOPN', 'FICHIER CREE AVEC TABLEBURP NON-OFFICIELLE', WARNIN, NONOFTB)
         IF (MESSNIV .LE. WARNIN) THEN
           WRITE(0, 703)
           WRITE(0, 702)
         ENDIF
      ENDIF

      IF(MESSNIV .LE. INFORM) THEN
         IF(INDEX(MODE,'CREATE') .NE. 0) WRITE(IOUT, 300) IUN
         WRITE(IOUT, 400) IUN
      ENDIF
      MRFOPN = NOMBRE
        
 300  FORMAT(/' UNITE = ',I3,' FICHIER RAPPORT EST CREE')
 400  FORMAT(/' UNITE = ',I3,' FICHIER RAPPORT EST OUVERT')
 700  FORMAT(' ***********************ATTENTION***********************')
 702  FORMAT(' *******************************************************')
 703  FORMAT(' *')

      RETURN
      END
