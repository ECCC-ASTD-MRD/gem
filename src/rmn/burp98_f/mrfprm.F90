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
!.S MRFPRM
!**S/P MRFPRM - OBTENIR LES PARAMETRES PRINCIPAUX D'UN RAPPORT
      FUNCTION MRFPRM(HANDLE, STNID, IDTYP, LAT, LON, DX, DY, DATE, TEMPS,  FLGS,  SUP,   NSUP,     LONENR)
      IMPLICIT NONE
      INTEGER  MRFPRM, NSUP, IDTYP, LAT, DX, DATE, SUP(*), LONENR, HANDLE, FLGS, TEMPS, LON, DY 
      CHARACTER*9 STNID
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!REV 003 M. Lepine   Avr     2000 appel a rah2char au lieu de write et hljust
!
!OBJET( MRFPRM )
!     ALLER CHERCHER LES PARAMETRES PRINCIPAUX DE L'ENREGISTREMENT
!     DONT LE POINTEUR EST HANDLE.
!
!ARGUMENTS
!
!     HANDLE   ENTREE  POINTEUR A L'ENREGISTREMENT
!     NSUP        "    NOMBRE DE DESCRIPTEURS SUPPLEMENTAIRES 
!                      (DOIT ETRE EGAL A ZERO POUR VERSION 1990)
!     STNID    SORTIE  IDENTIFICATEUR DE LA STATION
!     IDTYP       "    TYPE DE RAPPORT
!     LAT         "    LATITUDE EN CENTIDEGRES PAR RAPPORT AU POLE SUD
!     LON         "    LONGITUDE EN CENTIDEGRES (0-35999)
!     DX          "    DIMENSION X D'UNE BOITE
!     DY          "    DIMENSION Y D'UNE BOITE
!     DATE        "    DATE DE VALIDITE (AAMMJJHH)
!     TEMPS       "    HEURE DE L'OBSERVATION
!     FLGS        "    MARQUEURS GLOBAUX
!     SUP         "    LISTE DE DESCRIPTEURS SUPPLEMENTAIRES (AUCUN)
!     LONENR      "    LONGUEUR DE L'ENREGISTREMENT EN MOTS HOTE
!
!IMPLICITES
#include "codes.cdk"
#include "defi.cdk"
#include <ftnmacros.hf>
#include "enforc8.cdk"
!
!MODULES 
      EXTERNAL QDFERR, XDFPRM, RAH2CHAR
      INTEGER  QDFERR, XDFPRM, ADDR, LNGR, TYPREC, PRI(NPRITOT), NPRI,   I
      integer annee, mois, AA, MM, JJ
!
!*
      NPRI = NPRIDEF

!     POUR LA CUVEE 90, NSUP DOIT ETRE EGAL A ZERO
      IF(NSUP .GT. NPRISUP) THEN
         MRFPRM = QDFERR('MRFPRM', 'IL Y A TROP DE CLEFS PRIMAIRES SUPPLEMENTAIRES', WARNIN, ERCLEF)
         NSUP = NPRISUP
      ENDIF

!     POUR LES VERSIONS SUBSEQUENTES, NSUP PEUT ETRE PLUS GRAND QUE ZERO.
!     ON AJOUTE ALORS LES CLEFS PRIMAIRES SUPPLEMENTAIRES AU VECTEUR PRI.
      IF(NSUP .GT. 0) NPRI = NPRI + NSUP

!     OBTENIR LES CLEFS PRIMAIRES
      MRFPRM = XDFPRM(HANDLE, ADDR, LNGR, TYPREC, PRI, NPRI)
      IF(MRFPRM .LT. 0) RETURN

!     CALCUL DE LA LONGUEUR
      LONENR = LNGR * UNITES/BITMOT

!     DECOMPOSER LES CLEFS PRIMAIRES
      DO 10 I = 1,9
!         PRI(I) = HLJUST(PRI(I),1)
!         WRITE(STNID(I:I),'(A1)') PRI(I)
         CALL RAH2CHAR(STNID(I:I),PRI(I),1)
 10      CONTINUE

      FLGS  = PRI(10)
      LAT   = PRI(11)
      LON   = PRI(12)
      DATE  = PRI(13)
      if ((mod((date/100),100) .gt. 12) .or. (enforc8)) then
!
!     retourner la date en format AAAAMMJJ
!
        AA = mod((date/10000),100)
        MM = mod((date/100),100)
        JJ = mod(date,100)
        annee = 1900 + AA + (((MM-1)/12)*100)
        mois = 1 + mod(MM-1,12)
        date = (annee * 10000) + (mois * 100) + JJ
      endif
      DX    = PRI(14)
      IDTYP = PRI(15)
      DY    = PRI(16)
      TEMPS = (PRI(17)*100) + PRI(18)
      
!     OBTENIR LES CLEFS PRIMAIRES SUPPLEMENTAIRES
      IF(NSUP .GT. 0) THEN
         DO 20 I = 1,NSUP
            SUP(I) = PRI(NPRIDEF+I)
 20         CONTINUE
      ENDIF

      MRFPRM = 0

      RETURN
      END
