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
!.S MRBINI
!**S/P MRBINI - INITIALISER L'ENTETE D'UN RAPPORT
!
      FUNCTION MRBINI(IUN, BUF, TEMPS, FLGS, STNID, IDTYP, LATI, LONG, DX, DY, ELEV, IDRCV, DATEin, OARS, RUN, SUP, NSUP, XAUX, NXAUX)
      IMPLICIT NONE
      INTEGER  FLGS, DX, NXAUX, IDTYP, LONG,   &
     &         DATEin, RUN, BUF(*), XAUX(*),   &
     &         ELEV, DY, TEMPS, IDRCV, LATI, OARS, IUN, SUP(*), MRBINI,   &
     &         NSUP
      CHARACTER*(*) STNID
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!REV 003 M. Lepine   Avr     2000 appel a char2rah au lieu de read et hrjust
!
!OBJET( MRBINI )
!     INITIALISER L'ENTETE D'UN RAPPORT.  AVANT DE METTRE QUOI QUE 
!     CE SOIT DANS UN RAPPORT, ON INITIALISE LES DIFFERENTES CLEFS
!     PRIMAIRES ET AUXILIAIRES.
!
!ARGUMENTS
!     IUN     ENTREE   NUMERO D'UNITE ASSOCIE AU FICHIER
!     TYPREC    "      TYPE D'ENREGISTREMENT
!     IDELT     "      DIFF DE TEMPS ENTRE T VALIDITE ET T SYNOPTIQUE
!     FLGS      "      MARQUEURS GLOBAUX
!     STNID     "      IDENTIFICATEUR DE LA STATION
!     IDTYP     "      TYPE DE RAPPORT
!     LATI      "      LATITUDE DE LA STATION EN CENTIDEGRES
!     LONG      "      LONGITUDE DE LA STATION EN CENTIDEGRES
!     DX        "      DIMENSION X D'UNE BOITE
!     DY        "      DIMENSION Y D'UNE BOITE
!     ELEV      "      ALTITUDE DE LA STATION EN METRES
!     IDRCV     "      DELAI DE RECEPTION
!     DATEin    "      DATE SYNOPTIQUE DE VALIDITE (AAMMJJHH)
!     OARS      "      RESERVE POUR ANALYSE OBJECTIVE
!     RUN       "      IDENTIFICATEUR DE LA PASSE OPERATIONNELLE
!     SUP       "      CLEFS PRIMAIRES SUPPLEMENTAIRES 
!                      (AUCUNE POUR LA VERSION 1990)
!     NSUP      "      NOMBRE DE CLEFS PRIMAIRES SUPPLEMENTAIRES 
!                      (DOIT ETRE ZERO POUR LA VERSION 1990)
!     XAUX      "      CLEFS AUXILIAIRES SUPPLEMENTAIRES (=0 VRSN 1990)
!     NXAUX     "      NOMBRE DE CLEFS AUXILIAIRES SUPPLEMENTAIRES(=0)
!     BUF       "      VECTEUR QUI CONTIENDRA LES ENREGISTREMENTS
!
!IMPLICITES
#include "defi.cdk"
#include "codes.cdk"
#include "enforc8.cdk"
!
!MODULES 
      EXTERNAL XDFINI, CHAR2RAH, QDFERR
      INTEGER  XDFINI, QDFERR, KLPRIM(NPRITOT), KLAUX(NAUXTOT), NKLAUX, TYPREC, NKLPRIM, I
      integer AA, MM, JJ, annee, date
      CHARACTER*9 ISTNID
!
!*
      date = DATEin
      MRBINI  = -1
      NKLPRIM = NPRIDEF
      NKLAUX  = NAUXDEF
      TYPREC  = 1

!     POUR LA VERSION 1990, NSUP ET NXAUX DOIVENT ETRE EGAL A ZERO
      IF(NSUP .GT. NPRISUP) THEN
         MRBINI = QDFERR('MRBINI', 'IL Y A TROP DE CLEFS PRIMAIRES SUPPLEMENTAIRES', WARNIN, ERCLEF)
         NSUP = NPRISUP
      ENDIF
      IF(NXAUX .GT. NAUXSUP) THEN
         MRBINI = QDFERR('MRBINI', 'IL Y A TROP DE CLEFS AUXILIAIRES SUPPLEMENTAIRES', WARNIN, ERCLEF)
         NXAUX = NAUXSUP
      ENDIF


!     TRANSFORMER CHAQUE CARACTERE DE STNID EN UN ENTIER
      ISTNID = STNID
      DO 10 I = 1,9
!         READ(ISTNID(I:I),'(A1)') KLPRIM(I)
!         KLPRIM(I) = HRJUST(KLPRIM(I),1)
         CALL CHAR2RAH(ISTNID(I:I),KLPRIM(I),1)
10       CONTINUE

!     COMPOSER LES AUTRES CLEFS PRIMAIRES
      KLPRIM(10) = FLGS
      KLPRIM(11) = LATI
      KLPRIM(12) = LONG
      if (enforc8) then
         if (date .lt. 999999) then
         MRBINI = QDFERR('MRBINI', 'LA DATE DOIT ETRE EN FORMAT AAAAMMJJ', ERFATAL,ERRDAT)
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
20          CONTINUE           
         NKLPRIM = NKLPRIM + NSUP
      ENDIF

!     COMPOSER LES CLEFS AUXILIAIRES
      KLAUX(1) = 0
      KLAUX(2) = OARS
      KLAUX(3) = ELEV
      KLAUX(4) = IDRCV
      KLAUX(5) = RUN

!     AJOUTER LES CLEFS AUXILIAIRES SUPPLEMENTAIRES
      IF(NXAUX .GT. 0) THEN
         DO 30 I = 1,NXAUX
            KLAUX(NAUXDEF+I) = XAUX(I)
30          CONTINUE           
         NKLAUX = NKLAUX + NXAUX
      ENDIF

!     INITIALISER LE TOUT
      MRBINI = XDFINI(IUN, BUF, TYPREC, KLPRIM, NKLPRIM, KLAUX, NKLAUX)

!     INITIALISER BUF(8) AU NOMBRE DE BLOCS (=0)
!      BUF(8) = 0

!     METTRE LE BIT DE DEBUT DES BLOCS DANS BUF(9)
!      BUF(9) = 0
      call buf89a0(buf)

      RETURN
      END
