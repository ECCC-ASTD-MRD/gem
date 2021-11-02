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
!.S MRBPRM
!**S/P MRBPRM - EXTRAIRE LES PARAMETRES DESCRIPTEURS D'UN BLOC
!
      FUNCTION MRBPRM(BUF,  BKNO, NELE, NVAL, NT, BFAM, BDESC, BTYP, NBIT, BIT0, DATYP)
      IMPLICIT NONE
      INTEGER  MRBPRM, BUF(*), BKNO, NELE, NVAL, NT, BFAM, BDESC, BTYP, NBIT,   BIT0,   DATYP
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!rev 002 j. caveen   sept    1995 elimine bdsec et bfam passe a 12 bits
!
!OBJET( MRBPRM )
!     FONCTION SERVANT A RETOURNER LES PARAMETRES DESCRIPTEURS
!     DU BLOC DE DONNES DONT LE NUMERO D'ORDRE EST BKNO.
!                                                                       
!ARGUMENTS
!     BUF     ENTREE  VECTEUR CONTENANT LE RAPPORT
!     BKNO       "    NUMERO D'ORDRE DU BLOC
!     NELE    SORTIE  NOMBRE D'ELEMENTS METEOROLOGIQUES
!     NVAL       "    NOMBRE DE VALEURS PAR ELEMENTS
!     NT         "    NOMBRE DE NELE*NVAL VALEURS
!     BFAM       "    FAMILLE DU BLOC (12 bits)
!     BDESC      "    DESCRIPTEUR DU BLOC (mis a zero)
!     BTYP       "    TYPE DU BLOC
!     NBIT       "    NOMBRE DE BITS CONSERVES PAR VALEUR
!     BIT0       "    PREMIER BIT DU TABLEAU DE VALEURS
!     DATYP      "    TYPE DE COMPACTION
!
!IMPLICITES
#include "defi.cdk"
#include "bpl.cdk"
#include <ftnmacros.hf>
!
!MODULE 
      EXTERNAL XDFXTR
      INTEGER  XDFXTR
!
!*
      integer  ENTETE(DIMENT), BITPOS
      integer famdesc

      MRBPRM = -1

!     EXTRAIRE L'ENTETE DU BLOC
      BITPOS = (BKNO-1)*NBENTB
      MRBPRM = XDFXTR(BUF, ENTETE, BITPOS, DIMENT, BITMOT, 0)

!     EXTRAIRE DU BLOC LES DIFFERENTS PARAMETRES
      famdesc= GETBIT(ENTETE, BPFMDSC,LFMDSC)
      BTYP   = GETBIT(ENTETE, BPTYP,   LTYP)
      NBIT   = GETBIT(ENTETE, BPNBIT,  LNBIT) + 1
      BIT0   = GETBIT(ENTETE, BPBIT0,  LBIT0)
      DATYP  = GETBIT(ENTETE, BPDATYP, LDATYP)
      NELE   = GETBIT(ENTETE, BPNELE,  LNELE)

      IF(NELE .GE. GRONELE) THEN
         NELE = GETBIT(ENTETE, BPNELE16, LNELE16)
         NVAL = GETBIT(ENTETE, BPNVAL16, LNVAL16)
         NT   = GETBIT(ENTETE, BPNT16,   LNT16)
      ELSE
         NVAL = GETBIT(ENTETE, BPNVAL, LNVAL)
         NT   = GETBIT(ENTETE, BPNT,   LNT)
      ENDIF

!
!     construire bfam a partir de famdesc 
!     (interchange 6 bits du bas avec 6 bits du haut)
!
      bfam = LSHIFT(IAND(famdesc,RMASK(6)),6)
      bfam = IOR(bfam,(IAND(RSHIFT(famdesc,6),RMASK(6))))

      bdesc = 0

      MRBPRM = 0

      RETURN
      END
