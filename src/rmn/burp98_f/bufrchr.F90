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

!**S/P BUFRCHR  - TROUVER L'INDEX POINTANT A UN ELEMENT DANS UN TABLEAU
!
      FUNCTION BUFRCHR(ELEM, TABLEAU, NELE)
      IMPLICIT NONE
      INTEGER BUFRCHR, ELEM, NELE, TABLEAU(3, NELE)
!
!AUTEUR  J. CAVEEN   JANVIER 1991
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( BUFRCHR )
!     FAIRE UNE RECHERCHE BINAIRE SUR LES NOMS D'ELEMENTS CONTENUS DANS
!     LA MATRICE TABLEAU ET RETOURNER L'INDEX DE L'ELEMENT DE TABLEAU
!     CONTENANT LA VARIABLE ELEM.
!
!ARGUMENTS
!     ELEM      ENTREE  NOM DE L'ELEMENT A CHERCHER
!     NELE         "    NOMBRE D'ELEMENTS DANS TABLEAU
!     TABLEAU      "    MATRICE OU ON RECHERCHE ELEM
!
      INTEGER FIN, DEBUT, MILIEU
!
!*

      BUFRCHR = -1
      DEBUT   = 0
      FIN     = NELE + 1

  10  MILIEU  = (DEBUT + FIN)/2

      IF(DEBUT .NE. MILIEU) THEN
         IF(ELEM .NE. TABLEAU(1, MILIEU)) THEN
            IF(ELEM .GT. TABLEAU(1, MILIEU)) THEN
               DEBUT = MILIEU
            ELSE
               FIN   = MILIEU 
	    ENDIF
            GO TO 10
         ENDIF
         BUFRCHR = MILIEU
      ENDIF

      RETURN
      END

