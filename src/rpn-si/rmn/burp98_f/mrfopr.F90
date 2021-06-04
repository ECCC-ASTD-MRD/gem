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
!.S MRFOPR
!**S/P MRFOPR - INITIALISER UNE OPTION DE TYPE REEL
!
      FUNCTION MRFOPR(OPTNOM, OPVALR)
      IMPLICIT NONE
      INTEGER  MRFOPR
      CHARACTER*(*)   OPTNOM
      REAL                    OPVALR
!
!AUTEUR  J. CAVEEN   AVRIL  1991
!REV 001 Y. BOURASSA MARS   1995 RATFOR @ FTN77
!
!OBJET( MRFOPR )
!     FONCTION SERVANT A INITIALISER UNE OPTION DE TYPE REEL
!     LA VALEUR DE L'OPTION EST CONSERVEE DANS LE COMMON BURPUSR
!
!ARGUMENTS
!     OPTNOM   ENTREE  NOM DE L'OPTION A INITIALISER
!     OPVALR     "     VALEUR A DONNER A L'OPTION
!
!IMPLICITES
#include "defi.cdk"
#include "burpopt.cdk"
#include "codes.cdk"
!
!MODULES
      EXTERNAL QDFERR
      INTEGER  QDFERR
      DATA     MANQUE /RMANQUE/
!
!*
      MRFOPR = -1
      IF(INDEX(OPTNOM, 'MISSING') .NE. 0) THEN
         MANQUE = OPVALR
         MRFOPR = 0
      ELSE
         MRFOPR = QDFERR('MRFOPR', 'NOM D''OPTION INCONNU', ERROR, EROPTN)
      ENDIF

      RETURN
      END
