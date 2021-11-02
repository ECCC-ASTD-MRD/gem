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
!.S    MRFNBR
!**S/P MRFNBR -  OBTENIR LE NOMBRE D'ENREGISTREMENTS DANS UN FICHIER
!
      FUNCTION MRFNBR( IUN )
      IMPLICIT NONE
      INTEGER  MRFNBR, IUN
!
!AUTEUR  J. CAVEEN   MAI   1991
!REV 001 Y. BOURASSA MARS  1995 RATFOR @ FTN77
!
!OBJET( MRFNBR )
!     FONCTION RETOURNANT LE NOMBRE D'ENREGISTREMENTS ACTIFS CONTENUS 
!     DANS UN FICHIER RAPPORT
!                                                                       
!ARGUMENTS
!
!     IUN      ENTREE  NUMERO DU FICHIER A OUVRIR
!                                                         
!IMPLICITES
#include "defi.cdk"
!
!MODULES 
      EXTERNAL XDFSTA
      INTEGER  XDFSTA, STAT(NSTATMAX), PRII, AUXX
      CHARACTER*4 APPL, VERSN
!
!*

      MRFNBR = -1

!     OBTENIR LES INFORMATIONS CONCERNANT LE FICHIER.
      MRFNBR = XDFSTA(IUN, STAT, NSTATMAX, PRII, 0, AUXX, 0, VERSN, APPL)
        
      IF(MRFNBR .GE. 0) MRFNBR = STAT(12)
        
      RETURN
      END
