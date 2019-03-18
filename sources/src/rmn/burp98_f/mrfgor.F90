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
!.S MRFGOR
!**S/P MRFGOR - OBTENIR LA VALEUR D'UNE OPTION DE TYPE REEL
!
      FUNCTION MRFGOR(OPTNOM, OPVALR)
      IMPLICIT NONE
      INTEGER  MRFGOR
      CHARACTER*(*)  OPTNOM
      REAL                   OPVALR
!
!AUTEUR  J. CAVEEN   AVRIL 1991
!REV 001 Y. BOURASSA MARS  1995 RATFOR @ FTN77
!
!OBJET( MRFGOR )
!     FONCTION SERVANT A OBTENIR LA VALEUR D'UNE OPTION DE TYPE REEL
!     LA VALEUR DE L'OPTION EST CONSERVEE DANS LE COMMON BURPUSR
!
!ARGUMENTS
!     OPTNOM     ENTREE     NOM DE L'OPTION A INITIALISER
!     OPVALR     SORTIE     VALEUR  DONNEE A L'OPTION
!
!IMPLICITES
#include "defi.cdk"
#include "burpopt.cdk"
#include "codes.cdk"
!
!MODULE 
      EXTERNAL QDFERR
      INTEGER  QDFERR
!
!* 
      MRFGOR = -1
      IF(INDEX(OPTNOM, 'MISSING') .NE. 0) THEN
         OPVALR =  MANQUE 
         MRFGOR = 0
      ELSE
         MRFGOR = QDFERR('MRFGOR', 'NOM D''OPTION INCONNU', ERROR, EROPTN)
      ENDIF

      RETURN
      END
