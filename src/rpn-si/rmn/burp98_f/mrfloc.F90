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
!.S MRFLOC
!**S/P MRFLOC - TROUVER UN RAPPORT DANS UN FICHIER BURP
! 
      FUNCTION MRFLOC(IUN, HANDLE, STNID, IDTYP, LAT, LON, DATEin, TEMPS, SUP, NSUP)
      IMPLICIT NONE
      INTEGER  MRFLOC, IUN, HANDLE, IDTYP, LAT, LON, DATEin, TEMPS, NSUP,  SUP(*)
      CHARACTER*(*) STNID
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 j. caveen   sept.   1995 ajout d'un appel a mrfprm pour produire 
!                                 un message plus explicite
!REV 003 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!REV 004 M. Lepine   Avr     2000 appel a char2rah au lieu de read et hrjust
!REV 005 M. Lepine   Jan     2003 date est un argument d'entree seulement
!
!OBJET( MRFLOC )
!     TROUVER LE POINTEUR DU RAPPORT DECRIT PAR LES PARAMETRES STNID,
!     IDTYP,LAT,LON,DATE ET LE CONTENU DU TABLEAU SUP.  LA RECHERCHE SE
!     FAIT A PARTIR DE L'ENREGISTREMENT DONT LE POINTEUR EST HANDLE.
!     SI UN ELEMENT DE STNID = '*', CET ELEMENT SERA IGNORE POUR LA
!     RECHERCHE.  UNE VALEUR DE -1 POUR LES AUTRES ARGUMENTS A LE
!     MEME EFFET.  SI LA VALEUR DE HANDLE EST ZERO, ON AFFECTUE LA
!     RECHERCHE A PARTIR DU DEBUT DU FICHIER.
!       
!ARGUMENTS
!     IUN     ENTREE   NUMERO D'UNITE DU FICHIER
!     HANDLE    "      POINTEUR A L'ENREGISTREMENT D'OU PART LA RECHERCHE
!     STNID     "      IDENTIFICATEUR DE LA STATION
!     IDTYP     "      TYPE DE RAPPORT
!     LAT       "      LATITUDE DE LA STATION
!     LON       "      LONGITUDE DE LA STATION
!     DATE      "      DATE DE VALIDITE DU RAPPORT
!     TEMPS     "      HEURE DE L'OBSERVATION
!     SUP       "      TABLEAU DE CLEFS DE RECHERCHES SUPPLEMENTAIRES
!     NSUP      "      NOMBRE DE CLEFS SUPPLEMENTAIRES
!
!
!IMPLICITES
#include "codes.cdk"
#include "defi.cdk"
#include "burpopt.cdk"
#include "enforc8.cdk"
!
!MODULES 
      EXTERNAL QDFERR, XDFLOC, CHAR2RAH, mrfprm
      INTEGER  QDFERR, XDFLOC, mrfprm
!*
      CHARACTER*9  ISTNID
      INTEGER IIDTYP, ILAT, ILON, IDATE, ITEMPS, ISUP(1), INSUP
      INTEGER IDX,IDY,IFLGS,ILNGR,irien

      INTEGER PRI(NPRITOT), NPRI, I
      integer annee, mois, AA, MM, JJ, date

      date = DATEin
      MRFLOC = -1
      NPRI   = NPRIDEF

!     POUR LA VERSION 1990, LES CLEFS SUPPLEMENTAIRES NE SONT PAS PERMISES
      IF(NSUP .GT. NPRISUP) THEN
         MRFLOC = QDFERR('MRFLOC','IL Y A TROP DE CLEFS PRIMAIRES', ' SUPPLEMENTAIRES', WARNIN, ERCLEF)
         NSUP = NPRISUP
      ENDIF

!     COMPOSER LES CLEFS A PARTIR DU STNID
      ISTNID = STNID
      DO 10 I=1,9
         IF(ISTNID(I:I) .NE. '*') THEN
!            READ(ISTNID(I:I),'(A1)') PRI(I)
!            PRI(I) = HRJUST(PRI(I), 1)
            CALL CHAR2RAH(ISTNID(I:I),PRI(I),1)
         ELSE
            PRI(I) = -1
         ENDIF
 10      CONTINUE

!     COMPOSER LE RESTE DES CLEFS DE RECHERCHE
      PRI(10) = -1
      PRI(11) = LAT
      PRI(12) = LON
      if ((enforc8) .and. (date .ne. -1)) then
         if (date .lt. 999999) then
         MRFLOC = QDFERR('MRFLOC', 'LA DATE DOIT ETRE EN FORMAT AAAAMMJJ', ERFATAL,ERRDAT)
         endif
      endif
      if (date .gt. 999999) then
         annee = date/10000
         AA = mod((date/10000),100)
         MM = (((annee - 1900) /100) * 12) + mod((date/100),100)
         JJ = mod(date,100)
         date = (AA * 10000) + (MM * 100) + JJ
      endif
      PRI(13) = DATE
      PRI(14) = -1
      PRI(15) = IDTYP
      PRI(16) = -1
      IF(TEMPS .EQ. -1) THEN
         PRI(17) = -1
      ELSE
         PRI(17) = TEMPS/100
      ENDIF
      PRI(18) = -1

!     INCLURE LES CLEFS SUPPLEMENTAIRES
      IF(NSUP .GT. 0) THEN
         DO 20 I=1,NSUP
            PRI(NPRIDEF+I) = SUP(I)
 20      CONTINUE
         NPRI = NPRI + NSUP
      ENDIF

!     TROUVER L'ENREGISTREMENT
      MRFLOC = XDFLOC(IUN, HANDLE, PRI, NPRI)
      IF (MESSNIV .LE. INFORM) THEN
         IF(MRFLOC .LT. 0) THEN
            WRITE(6,1000) STNID, IDTYP, LAT, LON, DATEin, TEMPS
         ELSE
            IRIEN = MRFLOC
            INSUP = 0
            IRIEN = MRFPRM(IRIEN,ISTNID, IIDTYP, ILAT, ILON, IDX,IDY, IDATE, ITEMPS,IFLGS, ISUP, INSUP,ILNGR)
            WRITE(6,1100) ISTNID, IIDTYP, ILAT, ILON, IDX,IDY, IDATE, ITEMPS,IFLGS,ILNGR
         ENDIF
      ENDIF

 1000 FORMAT(' MRFLOC- INEXISTANT - STNID=',A9,' IDTYP=',I3, ' LAT=',I5,' LON=',I5,' DATE=',I8,' TEMPS=',I4)
 1100 FORMAT(' MRFLOC- TROUVE - STNID=',A9,' IDTYP=',I3, ' LAT=',I5,' LON=',I5,' DX=',i4,' DY=',i4,' DATE=',I8, ' TEMPS=',I4,' FLGS=',i8,' LNGR=',i6)
      RETURN
      END
