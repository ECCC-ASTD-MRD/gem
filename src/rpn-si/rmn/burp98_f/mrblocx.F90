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
 
!.S MRBLOCX
!**S/P MRBLOCX - TROUVER UN BLOC DANS UN RAPPORT
!
      FUNCTION MRBLOCX( BUF,   BFAM, BDESC, BKNAT, BKTYP, BKSTP, BLKNO)
      IMPLICIT NONE
      INTEGER  MRBLOCX, BUF(*),BFAM, BDESC, BKNAT, BKTYP, BKSTP, BLKNO
!
!AUTEUR  J. CAVEEN     JANVIER 1992
!REV 001 YVON BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( MRBLOCX )
!     FONCTION SERVANT A BATIR UNE CLEF DE RECHERCHE BTYP A
!     PARTIR DE BKNAT, BKTYP ET BKSTP POUR L'APPEL SUBSEQUENT A MRBLOC
!     POUR CHAQUE CLEF D'ENTREE, ON TRANSPOSE LA VALEUR DU BIT
!     28, 29 OU 30 RESPECTIVEMENT DANS INBTYP (CES BITS NE SONT ALLUMES
!     QUE SI LES CLEFS D'ENTREE SON MISE A -1)
!                                                                       
!ARGUMENTS
!     BUF     ENTREE  VECTEUR C ONTENANT LE RAPPORT
!     BFAM       "    FAMILLE DU BLOC RECHERCHE
!     BDESC      "    DESCRIPTION DU BLOC RECHERCHE
!     BKNAT      "    PORTION NATURE DU BTYP DE BLOC RECHERCHE
!     BKTYP      "    PORTION TYPE DU BTYP DE BLOC RECHERCHE
!     BKSTP      "    PORTION SOUS-TYPE DU BTYP DE BLOC RECHERCHE
!     BLKNO      "    BLOC D'OU PART LA RECHERCHE
!
!IMPLICITES
#include "masques.cdk"
#include "defi.cdk"
#include "bpl.cdk"
#include "burpopt.cdk"
#include "codes.cdk"
#include <ftnmacros.hf>
!
!MODULE
      EXTERNAL MRBLOC
      INTEGER  MRBLOC, INBTYP 
!
!*

!     CONSTRUIRE LA CLEF BTYP
      INBTYP = 0
      INBTYP = IAND(BKSTP, IOR(SGBKSTP, BKSTPMSK))
      INBTYP = IOR(INBTYP, IOR(IAND(SGBKTYP, BKTYP), LSHIFT(IAND(BKTYP,BKTYPMSK), BPBKTYP)))
      INBTYP = IOR(INBTYP,IOR(IAND(SGBKNAT, BKNAT), LSHIFT(IAND(BKNAT, BKNATMSK), BPBKNAT)))

      MRBLOCX = MRBLOC(BUF, BFAM, BDESC, INBTYP, BLKNO)

      RETURN
      END
