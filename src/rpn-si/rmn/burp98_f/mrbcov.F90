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
!.S MRBCOV
!**S/P MRBCOV - RETOURNER LA VALEUR D'UN ELEMENT EN SEIZE BITS
      FUNCTION MRBCOV( ELEM )
      IMPLICIT NONE
      INTEGER  MRBCOV, ELEM

!AUTEUR: J. CAVEEN   FEVRIER 1991
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( MRBCOV )
!     FONCTION RETOURNANT LA VALEUR D'UN ELEMENT DE TELLE SORTE
!     QU'IL PUISSE TENIR EN SEIZE BITS.
!     POUR UN ELEMENT AYANT LE FORMAT DECIMAL  ABBCCC, (A,B,C DE 0 A 9)
!     ON RETOURNE UN ENTIER CONTENANT A SUR DEUX BITS, BB SUR SIX BITS
!     ET CCC SUR HUIT BITS
!
!IMPLICITES
#include <ftnmacros.hf>
#include "defi.cdk"
#include "masques.cdk"
!
!ARGUMENT
!     ELEM    ENTREE    VALEUR DECIMALE DE L'ELEMENT
!
      INTEGER IELEM, IIBIT, VIBIT, VIIIBIT
!
!*
      IELEM   = ELEM
      IIBIT   = IELEM/100000
      IELEM   = MOD(IELEM, 100000)
      VIBIT   = IELEM/1000
      VIIIBIT = MOD(IELEM, 1000)
      MRBCOV  = VIIIBIT
      MRBCOV  = IOR(MRBCOV, LSHIFT(IAND(VIBIT, MSK6B), BP6B))
      MRBCOV  = IOR(MRBCOV, LSHIFT(IAND(IIBIT, MSK2B), BP2B))

      RETURN
      END
