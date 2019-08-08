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
!.S MRBPRML
!**S/P MRBPRML - EXTRAIRE LES PARAMETRES DESCRIPTEURS DE TOUS LES BLOCS
!
      FUNCTION MRBPRML( BUF,    INBKNO, TBLPRM, NPRM, INBLOCS)
      IMPLICIT NONE

      INTEGER  MRBPRML, BUF(*), INBKNO,  NPRM, INBLOCS, TBLPRM(NPRM, INBLOCS)
!
!AUTEUR  J. CAVEEN   OCTOBRE 1990
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!rev 002 j. caveen   sept    1995 elimine bdsec et bfam passe a 12 bits
!
!OBJET( MRBPRML )
!     FONCTION SERVANT A RETOURNER DANS LE TABLEAU TBLPRM
!     LES PARAMETRES DESCRIPTEURS DES INBLOCS BLOCS A PARTIR 
!     DU BLOC SUIVANT LE BLOC NUMERO BKNO.
!                                                                       
!ARGUMENTS
!     BUF        ENTREE    VECTEUR CONTENANT LE RAPPORT
!     INBKNO        "      NUMERO D'ORDRE DU PREMIER BLOC
!     NPRM          "      NOMBRE DE PARAMETRES A EXTRAIRE (DIM 1 DE TBLPRM)
!     INBLOCS       "      NOMBRE DE BLOCS DONT ON VEUT LES PARAMETRES
!     TBLPRM     SORTIE    TABLEAU CONTENANT LES PARAMETRES DES INBLOCS
!
!     STRUCTURE DE TBLPRM(NPRM,INBLOCS)
!     TBLPRM(1,I) - NUMERO DU BLOC I
!     TBLPRM(2,I) - NOMBRE D'ELEMENTS DANS LE BLOC I  (NELE)
!     TBLPRM(3,I) - NOMBRE DE VALEURS  PAR ELEMENT    (NVAL)
!     TBLPRM(4,I) - NOMBRE DE PAS DE TEMPS            (NT)
!     TBLPRM(5,I) - FAMILLE DU BLOC                   (BFAM) (12 bits)
!     TBLPRM(6,I) - DESCRIPTEUR DE BLOC               (BDESC) (mis a zero)
!     TBLPRM(7,I) - TYPE DU BLOC                      (BTYP)
!     TBLPRM(8,I) - NOMBRE DE BITS PAR ELEMENT        (NBIT)
!     TBLPRM(9,I) - NUMERO DU PREMIER BIT             (BIT0)
!     TBLPRM(10,I)- TYPE DE DONNEES POUR COMPACTION   (DATYP)
!
!IMPLICITES
#include "defi.cdk"
#include "bpl.cdk"
#include "codes.cdk"
#include <ftnmacros.hf>
!
!MODULES 
      EXTERNAL XDFXTR, QDFERR, getbuf8
      INTEGER  XDFXTR, QDFERR, getbuf8
!*

      integer ENTETE(DIMENT), BITPOS, NBLOCS, BKNO, NOBL
      integer famdesc
!
!*
      MRBPRML = -1

!     S'ASSURER QUE LES DIMENSIONS DU TABLEAU SONT ADEQUATES
      IF(NPRM .NE. NPRMMAX) THEN
         MRBPRML = QDFERR('MRBPRML', 'DIMENSIONS DE TBLPRM INCORRECTES', ERROR, ERNPRM)
         RETURN
      ENDIF

!     BLOC DE DEPART
      BKNO = MAX(INBKNO, 0)

!     NOMBRE DE BLOCS A EXTRAIRE
      NBLOCS = MIN(INBLOCS, getbuf8(buf))

!     EXTRAIRE TOUTES LES ENTETES DE BLOCS
      DO 10 NOBL = 1,NBLOCS
!        ADRESSE DU BLOC
         BITPOS = BKNO*NBENTB
!        EXTRACTION DE L'ENTETE DU BLOC
         MRBPRML = XDFXTR(BUF, ENTETE, BITPOS, DIMENT, BITMOT, 0)
         TBLPRM(1,NOBL)  = BKNO + 1
         famdesc         = GETBIT(ENTETE, BPFMDSC,   LFMDSC)
         TBLPRM(6,NOBL)  = 0
         TBLPRM(7,NOBL)  = GETBIT(ENTETE, BPTYP,   LTYP)
         TBLPRM(8,NOBL)  = GETBIT(ENTETE, BPNBIT,  LNBIT) + 1
         TBLPRM(9,NOBL)  = GETBIT(ENTETE, BPBIT0,  LBIT0)
         TBLPRM(10,NOBL) = GETBIT(ENTETE, BPDATYP, LDATYP)
         TBLPRM(2,NOBL)  = GETBIT(ENTETE, BPNELE,  LNELE)
         IF(TBLPRM(2,NOBL) .GE. GRONELE) THEN
            TBLPRM(2,NOBL) = GETBIT(ENTETE, BPNELE16, LNELE16)
            TBLPRM(3,NOBL) = GETBIT(ENTETE, BPNVAL16, LNVAL16)
            TBLPRM(4,NOBL) = GETBIT(ENTETE, BPNT16,   LNT16)
         ELSE
            TBLPRM(3,NOBL) = GETBIT(ENTETE, BPNVAL, LNVAL)
            TBLPRM(4,NOBL) = GETBIT(ENTETE, BPNT,   LNT)
         ENDIF

!     
!     construire bfam a partir de famdesc 
!     (interchange 6 bits du bas avec 6 bits du haut)
!     
         TBLPRM(5,NOBL) = LSHIFT(IAND(famdesc,RMASK(6)),6)
         TBLPRM(5,NOBL) = IOR(TBLPRM(5,NOBL),(IAND(RSHIFT(famdesc,6),RMASK(6))))
         
         BKNO = BKNO + 1
10       CONTINUE 

        MRBPRML = NBLOCS

        RETURN
        END
