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

!.S QRBNBDT
!**S/P QRBNBDT - DETERMINER NBITS ET DATYP POUR UN BLOC DE DONNEES
!
      FUNCTION QRBNBDT( NBIT, DATYP, TBLVAL, TBLDIM)
      IMPLICIT NONE
      INTEGER  QRBNBDT, NBIT, DATYP, TBLDIM, TBLVAL(*)
!
!AUTEUR  J. CAVEEN   DECEMBRE 1992
!REV 001 Y. BOURASSA MARS     1995 RATFOR @ FTN77
!REV 002 J. Caveen   NOVEMBRE 1995 BUG FIX: DATYP = 0
!REV 003 J. Caveen   Avril 1995 reduction de 1 bit
!REV 004 M. Lepine   Juin 1995 bug fix, reduction de 1 bit
!
!OBJET( QRBNBDT )
!     FONCTION SERVANT A DETERMINER LE NOMBRE DE BITS REQUIS AINSI
!     QUE LE DATYP A UTILISER POUR INSERER TBLVAL DANS UN BLOC DE
!     DONNEES.
!     METHODE: - TROUVER MIN,MAX ET MAX(ABS(MIN),MAX)) DE TBLVAL
!                VOIR SI ENTIERS SIGNES REQUIS,
!                CALCULER NOMBRE DE BITS REQUIS
!                ON NE RETOURNE JAMAIS MOINS DE PRECISION QUE
!                CE QUE L'USAGER A DEMANDE.
!
!     LA FONCTION RETOURNE UN CODE D'ERREUR SI LA PRECISION REQUISE
!     EST IRREALISABLE.  SINON, ON RETOURNE ZERO.
!                                                                       
!
!ARGUMENTS
!     NBIT     ENT/SRT  NOMBRE DE BIT A CONSERVER PAR VALEUR
!     DATYP    ENT/SRT  TYPE DE DONNES POUR LA COMPACTION
!     TBLVAL   ENTREE   TABLEAU DES VALEURS A ECRIRE (NELE*NVAL*NT)
!     TBLDIM      "     NOMBRE D'ELEMENTS DANS TBLVAL
!
!IMPLICITE
#include "codes.cdk"
!
!MODULES 
      EXTERNAL QDFERR
      INTEGER  QDFERR, I, TBLMAX, TBLMIN, ERREUR, SUPVAL(32)
!
!     SUPVAL:TABLEAU DE REFERENCE POUR TROUVER LE NOMBRE DE BITS REQUIS
      DATA SUPVAL /  1,          2,           4,           8,  &
     &              16,         32,          64,         128,  &
     &             256,        512,        1024,        2048,  &
     &            4096,       8192,       16384,       32768,  &
     &           65536,     131072,      262144,      524288,  &
     &         1048576,    2097152,     4194304,     8388608,  &
     &        16777216,   33554432,    67108864,   134217728,  &
     &       268435456,  536870912,  1073741824,  2147483647  / 

      logical nomanke
!
!*
      ERREUR = 0
      IF(NBIT .LE. 0) NBIT = 1

!     ON RETOURNE SI DATYP ET NBIT == MODE TRANSPARENT
!     Ou si datatyp == 0
      IF(((DATYP.EQ.2) .AND. (NBIT.EQ.32)) .OR. (DATYP .EQ. 0)) THEN
         QRBNBDT = 0
         RETURN
      ENDIF

 
#if defined (NEC)       
!
!     Verifier si datyp est plus grand que 5 --> pas permis sur NEC
!
      if(datyp .ge. 6) then
         qrbnbdt = QDFERR('QRBNBDT', 'DATYP=6,7,8 OU 9, PAS PERMIS SUR NEC',ERFATAL, EREDAT)
         return
      endif
#endif

!
!     SI datyp >=6 (real, real*8, complex, complex*8)
!     on met le nombre de bits a 32
!
      if(datyp .ge. 6) then
         nbit = 32
         qrbnbdt = 0
         return
      endif

!     SI DATYP = CARACTERES, NBIT = 8
      IF(DATYP.EQ.3 .OR. DATYP.EQ.5) THEN
         NBIT    = 8
         QRBNBDT = 0
         RETURN
      ENDIF

!     TROUVER LE MIN ET LE MAX DU CHAMP
      TBLMAX = 0
      TBLMIN = 0
      DO 10 I = 1,TBLDIM
         TBLMAX = MAX(TBLMAX, TBLVAL(I))       
         TBLMIN = MIN(TBLMIN, TBLVAL(I))       
 10      CONTINUE
        
!     DETERMINER LA VALEUR ABS MAXIMALE ET SI ENTIERS SIGNES SONT REQUIS
      IF(TBLMIN .LT. -1) THEN
         TBLMAX  = MAX(TBLMAX, ABS( TBLMIN ))
         DATYP   = 4
      ENDIF   

!
!     Verifier si il y a des valeurs manquantes et mettre indicateur 
!     a vrai si necessaire. On devra alors calculer 1 bit de plus
!
!      nomanke = .true.
!      do i = 1, tbldim
!         if(tblval(i) .eq. -1) then
!            nomanke = .false.
!            goto 15
!         endif
!      enddo
!
! 15   continue

!     DETERMINER LE NOMBRE DE BITS A UTILISER
!     si la valeur maximale occupe tous les bits requis, ou si
!     il y ades valeurs manquantes (=-1), on rajoute 1 bit
!
      IF(TBLMAX .GE. SUPVAL(NBIT)) THEN
         DO 20 I = NBIT+1, 32
            IF(TBLMAX .LT. SUPVAL(I)) THEN
              NBIT = I
              IF(TBLMAX .LT. SUPVAL(I) - 1) NBIT=NBIT-1
              GO TO 100
            ENDIF
 20         CONTINUE
         ERREUR = QDFERR('QDFNBDT', 'ON CODE AVEC NBIT=32 ET DATYP=2', WARNIN, ERCMPR)
      ENDIF
     
 100  CONTINUE

!     S'ASSURER QUE LES PARAMETRES SONT VALABLES:
!     SI NBIT = 32 ==> DATYP = 2
!     SI DATYP = 4, ON ALLOUE UN BIT DE PLUS POUR LE SIGNE
      IF(DATYP .EQ. 4) THEN
         NBIT = NBIT + 1
         IF(NBIT .GT. 31) THEN
            NBIT   = 32
            DATYP  = 2
            ERREUR = QDFERR('QDFNBDT', ' ON CODE VALEURS <0  AVEC NBIT=32 ET DATYP=2', WARNIN, ERCMPR)
         ENDIF
      ELSE
         NBIT = MIN(NBIT, 32)
      ENDIF

      QRBNBDT = ERREUR

      RETURN
      END
