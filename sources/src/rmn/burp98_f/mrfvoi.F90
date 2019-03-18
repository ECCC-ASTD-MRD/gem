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
!.S MRFVOI
!**S/P MRFVOI - VOIR LE CONTENU D'UN FICHER BURP
!
      FUNCTION MRFVOI( IUN )
      IMPLICIT NONE
      INTEGER  MRFVOI, IUN
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!REV 002 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!REV 003 M. Lepine   Avr     2000 appel a rah2char au lieu de write et hljust
!
!OBJET( MRFVOI )
!     FONCTION SERVANT A OBTENIR LES PARAMETRES DESCRIPTEURS DE CHAQUE
!     ENREGISTREMENT DU FICHIER BURP DONT LE NUMERO D'UNITE EST IUN.
!     MRFVOI RENSEIGNE AUSSI SUR LE TYPE D'ENREGISTREMENT, SA LONGUEUR
!     AINSI QUE SON ADRESSE EN UNITES DE 64 BITS.
!                                                                       
!ARGUMENTS
!
!     IUN   ENTREE   NUMERO D'UNITE DU FICHIER
!
!IMPLICIT
#include "defi.cdk"
#include "enforc8.cdk"
!
!MODULES 
      EXTERNAL XDFSTA, XDFPRM, QDFIND, RAH2CHAR, XDFLOC, XDFOPN, XDFCLS, QQQFNOM
      INTEGER  XDFSTA,  XDFPRM, QDFIND, TEMPS, LONG, STAT(12),   &
     &         QQQFNOM, INDFIC, HANDLE, TYPREC, NLIGN, ADDR, IOUT, J,  &
     &         NKLPRIM, NCPRM,  XDFLOC, PRIDEF(2,NPRIDEF),   XDFCLS,  &
     &         NBPAGES, LENGR,  XDFOPN, AUXDEF(2,NAUXDEF),       &
     &         KLPRIM(NPRIDEF), CPRM(NPRIDEF)
      CHARACTER*9  STNID
      CHARACTER*4  VERSION, APPLIC
      CHARACTER*50 NOMFIC,  TYPFIC
      LOGICAL      ONFERME 
      integer date, mois, annee, AA, MM, JJ
!
!*

      DATA IOUT /6/
      DATA CPRM /NPRIDEF*-1/

      NCPRM   = NPRIDEF
      ONFERME = .FALSE.
      MRFVOI  = -1
      INDFIC  = QDFIND( IUN )
      IF(INDFIC .GT. MAXFIL) THEN
         ONFERME = .TRUE.
         MRFVOI  = XDFOPN(IUN, 'READ', PRIDEF, NPRIDEF, AUXDEF, NAUXDEF, 'BURP')
         IF(MRFVOI .LT. 0) RETURN
      ENDIF

!     OBTENIR LE NOM DU FICHIER
      MRFVOI = QQQFNOM(IUN, NOMFIC, TYPFIC, LENGR)

!     OBTENIR LES STATISTIQUES SE RAPPORTANT AU FICHIER
      MRFVOI  = XDFSTA(IUN, STAT, 12, PRIDEF, NPRIDEF, AUXDEF, NAUXDEF, VERSION, APPLIC)
      NKLPRIM = STAT(7)

!     OBTENIR LES PARAMETRES DESCRIPTEURS
      HANDLE  = 0
      HANDLE  = XDFLOC(IUN, HANDLE, CPRM, NCPRM)
      NLIGN   = 60
      NBPAGES = 1
 10   IF(HANDLE .GE. 0) THEN
         MRFVOI = XDFPRM(HANDLE, ADDR, LONG, TYPREC, KLPRIM, NKLPRIM)
!        ENREGISTREMENT VALIDE?
         IF(TYPREC .NE. 255) THEN
!           EXTRAIRE LES CLEFS PRIMAIRES
            DO 20 J = 1,9
!               KLPRIM(J) = HLJUST(KLPRIM(J), 1)
!               WRITE(STNID(J:J),'(A1)') KLPRIM(J)
               CALL RAH2CHAR(STNID(J:J),KLPRIM(J),1)
 20            CONTINUE
!           COMMENCER NOUVELLE PAGE?
            IF(NLIGN .EQ. 60) THEN 
               WRITE(IOUT, 400) IUN, NOMFIC, NBPAGES
               WRITE(IOUT, 500)
               NLIGN   = 0
               NBPAGES = NBPAGES + 1
            ENDIF

            TEMPS = (KLPRIM(17)*100) + KLPRIM(18)
            date = klprim(13)
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
            WRITE(IOUT, 600) STNID, KLPRIM(11), KLPRIM(12), KLPRIM(14),  &
     &                       KLPRIM(16), KLPRIM(10), date, TEMPS,  &
     &                       KLPRIM(15), LONG,  ADDR
            NLIGN = NLIGN + 1
         ENDIF
         HANDLE = XDFLOC(IUN, HANDLE, CPRM, NCPRM)
         GO TO 10
      ENDIF

!     FERMER LE FICHIER SI ON L'A OUVERT

      IF( ONFERME ) MRFVOI = XDFCLS( IUN )

!     ECRIRE LES STATISTTIQUES SE RAPPORTANT AU FICHIER
      IF(NLIGN .GT. 46) WRITE(IOUT,400) IUN, NOMFIC, NBPAGES
      WRITE(IOUT,800)
      WRITE(IOUT,900)' TAILLE DU FICHIER                   ',STAT(1)
      WRITE(IOUT,900)' NOMBRE DE REECRITURES               ',STAT(2)
      WRITE(IOUT,900)' NOMBRE D''EXTENSIONS                ',STAT(3)
      WRITE(IOUT,900)' NOMBRE D''EFFACEMENTS               ',STAT(11)
      WRITE(IOUT,900)' NOMBRE D''ENREGISTREMENTS VALIDES   ',STAT(12)
      WRITE(IOUT,900)' TAILLE DU PLUS GROS ENREGISTREMENT  ',STAT(6)
      WRITE(IOUT,700)
      MRFVOI = 0

      RETURN

 400  FORMAT('1  MRFVOI  UNITE  ',I3,'  NOM ',A,T86,'  PAGE ',I3)
 500  FORMAT('0  STATION   LATI   LONG     DX     DY   FLGS(HEX)   DATE'  &
     &       ,'   TEMPS   IDTYP   LONGUEUR  ADRESSE '/)
#if defined (NEC)
 600  FORMAT(' ',A9,1X,4(I6,1X),4X,Z6,1X,I8,3X,I4,3X,I3,3X,I8,1X,I10)
#else
 600  FORMAT(' ',A9,1X,4(I6,1X),4X,Z6.6,1X,I8,3X,I4,3X,I3,3X,I8,1X,I10)
#endif
 700  FORMAT(/' N.B. DIMENSIONS ET ADRESSES EN UNITES DE 64 BITS'//)
 800  FORMAT('0 STATISTIQUES'//)
 900  FORMAT(A38, I10)

      END
