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
!.S MRBUPD
!**S/P MRBUPD - DONNER UNE VALEUR AUX CLEFS D'UN RAPPORT
!
      FUNCTION MRBUPD(IUN, BUF, TEMPS, FLGS, STNID, IDTYP, LATI, LONG, DX, DY, ELEV, DRCV, DATEin, OARS, RUN, SUP, NSUP, XAUX, NXAUX)
      IMPLICIT NONE
      INTEGER  MRBUPD, NSUP, NXAUX, IUN, BUF(*), TEMPS, FLGS, IDTYP,  &
     &         LATI,   LONG, ELEV,  DX,  DY,     DRCV,  DATEin, OARS,  &
     &         RUN,    SUP(*), XAUX(*)
      CHARACTER*(*) STNID
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!REV 003 M. Lepine   Avr     2000 appel a char2rah au lieu de read et hrjust
!
!OBJET( MRBUPD )
!     MISE A JOUR DE L'ENTETE D'UN RAPPORT. SEULES LES CLEFS QUI
!     N'ONT PAS POUR VALEUR -1 SERONT MISE A JOUR.  IL EN VA DE MEME
!     POUR CHAQUE CARACTERE DE STNID S'IL EST DIFFERENT DE '*'.
!
!ARGUMENTS
!     IUN     ENTREE  NUMERO D'UNITE ASSOCIE AU FICHIER
!     TYPREC    "     TYPE D'ENREGISTREMENT
!     TEMPS     "     DIFF DE TEMPS ENTRE T VALIDITE ET T SYNOPTIQUE
!     FLGS      "     MARQUEURS GLOBAUX
!     STNID     "     IDENTIFICATEUR DE LA STATION
!     IDTYP     "     TYPE DE RAPPORT
!     LATI      "     LATITUDE DE LA STATION EN CENTIDEGRES
!     LONG      "     LONGITUDE DE LA STATION EN CENTIDEGRES
!     DX        "     DIMENSION X D'UNE BOITE
!     DY        "     DIMENSION Y D'UNE BOITE
!     ELEV      "     ALTITUDE DE LA STATION EN METRES
!     DRCV      "     DELAI DE RECEPTION
!     DATE      "     DATE SYNOPTIQUE DE VALIDITE (AAMMJJHH)
!     OARS      "     RESERVE POUR ANALYSE OBJECTIVE
!     RUN       "     IDENTIFICATEUR DE LA PASSE OPERATIONNELLE
!     SUP       "     CLEFS PRIMAIRES SUPPLEMENTAIRES (AUCUNE POUR
!                     LA VERSION 1990)
!     NSUP      "     NOMBRE DE CLEFS PRIMAIRES SUPPLEMENTAIRES (DOIT
!                     ETRE ZERO POUR LA VERSION 1990)
!     XAUX      "     CLEFS AUXILIAIRES SUPPLEMENTAIRES (=0 VRSN 1990)
!     NXAUX     "     NOMBRE DE CLEFS AUXILIAIRES SUPPLEMENTAIRES(=0)
!     BUF       "     VECTEUR QUI CONTIENDRA LES ENREGISTREMENTS
!
!IMPLICITES
#include "defi.cdk"
#include "codes.cdk"
#include "enforc8.cdk"
!
!MODULES 
      integer getbuf8
      external getbuf8
      EXTERNAL XDFUPD, CHAR2RAH, QDFERR
      INTEGER  XDFUPD, QDFERR, KLPRIM(NPRITOT), NKLPRIM, I, NKLAUX, TYPREC, KLAUX(NAUXTOT)
      integer AA, MM, JJ, annee, date
      CHARACTER*9 ISTNID
!
!*
      date = DATEin
      MRBUPD  = -1
      NKLPRIM = NPRIDEF
      NKLAUX  = NAUXDEF
      TYPREC  = 1

!     POUR LA VERSION 1990, NSUP ET NXAUX DOIVENT ETRE EGAL A ZERO
      IF(NSUP .GT. NPRISUP) THEN
         MRBUPD = QDFERR('MRBUPD','IL Y A TROP DE CLEFS PRIMAIRES', ' SUPPLEMENTAIRES', WARNIN, ERCLEF)
         NSUP = NPRISUP
      ENDIF
      IF(NXAUX .GT. NAUXSUP) THEN
         MRBUPD = QDFERR('MRBUPD', 'IL Y A TROP DE CLEFS AUXILIAIRES', ' SUPPLEMENTAIRES', WARNIN, ERCLEF)
         NXAUX = NAUXSUP
      ENDIF

!     TRANSFORMER CHAQUE CARACTERE DE STNID EN UN ENTIER
      ISTNID = STNID
      DO 10 I = 1,9
         IF(ISTNID(I:I).NE.'*') THEN
!            READ(ISTNID(I:I),'(A1)') KLPRIM(I)
!            KLPRIM(I) = HRJUST(KLPRIM(I),1)
            CALL CHAR2RAH(ISTNID(I:I),KLPRIM(I),1)            
         ELSE
            KLPRIM(I) = -1
         ENDIF
 10      CONTINUE

!     COMPOSER LES AUTRES CLEFS PRIMAIRES
      KLPRIM(11) = LATI
      KLPRIM(12) = LONG
      KLPRIM(10) = FLGS
      if ((enforc8) .and. (date .ne. -1)) then
         if (date .lt. 999999) then
         MRBUPD = QDFERR('MRBUPD', 'LA DATE DOIT ETRE EN FORMAT AAAAMMJJ', ERFATAL,ERRDAT)
         endif
      endif
      if (date .gt. 999999) then
         annee = date/10000
         AA = mod((date/10000),100)
         MM = (((annee - 1900) /100) * 12) + mod((date/100),100)
         JJ = mod(date,100)
         date = (AA * 10000) + (MM * 100) + JJ
      endif
      KLPRIM(13) = DATE
      KLPRIM(14) = DX
      KLPRIM(15) = IDTYP
      KLPRIM(16) = DY
      IF(TEMPS .EQ. -1) THEN
         KLPRIM(17) = -1
         KLPRIM(18) = -1
      ELSE
         KLPRIM(17) = TEMPS/100
         KLPRIM(18) = MOD(TEMPS,100)
      ENDIF

!     AJOUTER LES CLEFS PRIMAIRES SUPPLEMENTAIRES
      IF(NSUP .GT. 0) THEN
         DO 20 I = 1,NSUP
            KLPRIM(NPRIDEF+I) = SUP(I)
 20         CONTINUE
         NKLPRIM = NKLPRIM + NSUP
      ENDIF

!     COMPOSER LES CLEFS AUXILIAIRES
      KLAUX(1) = getbuf8(BUF)
      KLAUX(2) = OARS
      KLAUX(3) = ELEV
      KLAUX(4) = DRCV
      KLAUX(5) = RUN

!     AJOUTER LES CLEFS AUXILIAIRES SUPPLEMENTAIRES
      IF(NXAUX .GT. 0) THEN
         DO 30 I = 1,NXAUX
            KLAUX(NAUXDEF+I) = XAUX(I)
 30         CONTINUE
         NKLAUX = NKLAUX + NXAUX
      ENDIF

!     INITIALISER LE TOUT
      MRBUPD = XDFUPD(IUN, BUF, TYPREC, KLPRIM, NKLPRIM, KLAUX, NKLAUX)

      RETURN
      END
