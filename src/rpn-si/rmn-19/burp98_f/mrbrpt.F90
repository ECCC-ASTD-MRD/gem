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

!.S MRBRPT
!**S/P  MRBRPT - VERIFIER SI UN ELEMENT EST REPETITIF OU NON
      FUNCTION MRBRPT( ELEMENT )
      IMPLICIT NONE
      INTEGER  MRBRPT, ELEMENT
!
!AUTEUR: J. CAVEEN    JANVIER 1992
!REV001  Y. BOURASSA  MARS    1995 RATFOR @ FTN77
!
!OBJET( MRBRPT )
!     FONCTION SERVANT A VERIFIER SI UN ELEMENT EST REPETITIF OU NON.
!     LA FONCTION RETOURNE:
!        1 - ELEMENT REPETITIF
!        0 - ELEMENT NON REPETITIF
!       <0 - CODE D'ELEMENT NON VALABLE
!            (PLUS PETIT QUE UN OU PLUS GRAND QUE MAXREP)
!
!ARGUMENT
!     ELEMENT  ENTREE  CODE DE L'ELEMENT A VERIFIER
!
!IMPLICITES
#include "defi.cdk"
#include "burpopt.cdk"
#include "codes.cdk"
#include <ftnmacros.hf>
!
!MODULE
      EXTERNAL QDFERR
      INTEGER  QDFERR
!
!*
      IF(ELEMENT.LT.1 .OR. ELEMENT .GT. MAXREP*BITMOT) THEN
         MRBRPT = QDFERR('MRBRPT', 'NOM D''ELEMENT NON VALIDE', WARNIN, ERELEM)
      ELSE
         MRBRPT = GETBIT(RPETITIF, ELEMENT, 1)
      ENDIF

      RETURN
      END
