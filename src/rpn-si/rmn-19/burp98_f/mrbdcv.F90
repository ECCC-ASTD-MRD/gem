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

!.S MRBDCV
!**S/P MRBDCV - RETOURNER LA VALEUR DECODEE D'UN ELEMENT
      FUNCTION MRBDCV( ELEM )
      IMPLICIT NONE
      INTEGER  MRBDCV, ELEM

!AUTEUR:  J. CAVEEN   FEVRIER 1991
!REV 001  Y. BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( MRBDCV )
!     FONCTION RETOURNANT LA VALEUR DECODEE D'UN ELEMENT QUI A ETE CODE
!     DE TELLE SORTE QU'IL PUISSE TENIR EN SEIZE BITS.
!
!     POUR UN ELEMENT, ON RETOURNE SA VALEUR SOUS FORMAT DECIMAL  ABBCCC,
!                                                        (A,B,C DE 0 A 9)
!     OU A    PROVIENT DES BITS 14 ET 15 DE L'ELEMENT
!        BB       "    DES BITS 8 A 13 DE L'ELEMENT
!        CCC      "    DES BITS 0 A 7  DE L'ELEMENT
!
!IMPLICITES
#include <ftnmacros.hf>
#include "defi.cdk"
#include "masques.cdk"
!
!ARGUMENT
!     ELEM    ENTREE    VALEUR CODEE DE L'ELEMENT
!
      INTEGER IELEM, IIBIT, VIBIT, VIIIBIT
!
!*
      IELEM   = ELEM
      VIIIBIT = IAND(IELEM, MSK8B)
      VIBIT   = IAND(RSHIFT(IELEM,BP6B), MSK6B)
      IIBIT   = IAND(RSHIFT(IELEM,BP2B), MSK2B)
      MRBDCV  = IIBIT * 100000 + VIBIT * 1000 + VIIIBIT

      RETURN
      END
