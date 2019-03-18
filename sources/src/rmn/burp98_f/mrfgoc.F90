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

!**S/P MRFGOC - INITIALISER LA VALEUR D'UNE OPTION DE TYPE CARACTERE
!
        FUNCTION MRFGOC(OPTNOM,OPVALC)
        implicit none
        INTEGER MRFGOC
        CHARACTER*(*) OPTNOM,OPVALC
        
!
!AUTEUR           J. CAVEEN AVRIL 1991
!Revision         james caveen - mai 1995 - convertir a fortran
!Revision  002    Mario Lepine - mai 2000 - remplace -1 par variable notused
!
!OBJET(MRFGOC)
!     FONCTION SERVANT A OBTENIR LA VALEUR D'UNE OPTION DE TYPE CARACTERE
!     LA VALEUR DE L'OPTION EST CONSERVEE DANS LE COMMON XDFTLR
!
!ARGUMENTS
!
!     OPTNOM     ENTREE     NOM DE L'OPTION A INITIALISER
!     OPVALR     SORTIE     VALEUR  DONNEE A L'OPTION
!
!
!MODULES 
        INTEGER XDFGOP
        EXTERNAL XDFGOP
        integer notused
!*

        MRFGOC = XDFGOP(OPTNOM,OPVALC,notused)

        RETURN
        END
