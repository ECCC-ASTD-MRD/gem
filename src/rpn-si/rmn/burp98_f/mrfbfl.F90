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
!**S/P MRFBFL - OBTENIR LA DIMENSION DU BUFFER NECESSAIRE POUR STOCKER
!               LE PLUS LONG RAPPORT CONTENU DANS LE FICHIER
!
      FUNCTION MRFBFL( IUN )
      IMPLICIT NONE
      INTEGER  MRFBFL, IUN
!
!AUTEUR  J. BLEZIUS   FEVRIER 2013
!
!OBJET( MRFMXL )
!     RETOURNER LA LONGUEUR DE L'ENREGISTREMENT LE PLUS LONG
!     CONTENU DANS LE FICHIER IUN + LE NOMBLE D'ELEMENTS EXTRA
!     NECESSAIRE A LA STRUCTURE INTERNE DU BUFFER
!
!ARGUMENT
!     IUN      ENTREE         NUMERO D'UNITE DU FICHIER
!
!IMPLICITES
#include <ftnmacros.hf>
!
!MODULES 
      EXTERNAL MRFMXL
      INTEGER  MRFMXL
!
!*

      MRFBFL = MRFMXL(IUN) + 10

      RETURN
      END
