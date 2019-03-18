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
!.S MRFMXL
!**S/P MRFMXL - OBTENIR LES PARAMETRES PRINCIPAUX D'UN RAPPORT
!
      FUNCTION MRFMXL( IUN )
      IMPLICIT NONE
      INTEGER  MRFMXL, IUN
!
!AUTEUR  J. CAVEEN   FEVRIER 1991
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( MRFMXL )
!     RETOURNER LA LONGUEUR DE L'ENREGISTREMENT LE PLUS LONG
!     CONTENU DANS LE FICHIER IUN
!
!ARGUMENT
!     IUN      ENTREE         NUMERO D'UNITE DU FICHIER
!
!IMPLICITES
#include "defi.cdk"
#include <ftnmacros.hf>
!
!MODULES 
      EXTERNAL XDFSTA
      INTEGER  XDFSTA, STAT(6), PRIDEF, AUXDEF
      CHARACTER*4 VERSN, APPL
!
!*

      MRFMXL = XDFSTA(IUN, STAT, 6, PRIDEF, 0, AUXDEF, 0, VERSN, APPL)

      IF(MRFMXL .LT. 0) RETURN

!     OBTENIR LA LONGUEUR DE STAT(6) ET LA CONVERTIR EN MOTS-HOTE
      MRFMXL = STAT(6) * UNITES/BITMOT

      RETURN
      END
