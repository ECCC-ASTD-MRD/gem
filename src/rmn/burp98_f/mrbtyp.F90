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

!.S MRBTYP
!**S/P MRBTYP - CONVERTIR BKNAT, BKTYP ET BKSTP A BTYP OU L'INVERSE
!
      FUNCTION MRBTYP(BKNAT, BKTYP, BKSTP, BTYP)
      IMPLICIT NONE
      INTEGER  MRBTYP, BTYP, BKNAT, BKTYP, BKSTP
!
!AUTEUR  J. CAVEEN   JANVIER 1992
!REV 001 Y. BOURASSA MARS    1995 RATFOR @ FTN77
!
!OBJET( MRBTYP )
!     FONCTION SERVANT A BATIR UNE CLEF DE RECHERCHE BTYP A
!     PARTIR DE BKNAT, BKTYP ET BKSTP OU A EXTRAIRE 
!     BKNAT, BKTYP ET BKSTP DE BTYP
!                                                                       
!ARGUMENTS
!     BTYP    ENT/SRT  CLEF COMPOSITE INDIQUANT LE TYPE DE BLOC
!     BKNAT      "     PORTION NATURE DU BTYP DE BLOC RECHERCHE
!     BKTYP      "     PORTION TYPE DU BTYP DE BLOC RECHERCHE
!     BKSTP      "     PORTION SOUS-TYPE DU BTYP DE BLOC RECHERCHE
!
!                      BTYP=0 - DE BKNAT, BKTYP, BKSTP -> BTYP
!                               FONCTION RETOURNE BTYP
!                      BTYP>0 - DE BTYP -> BKNAT, BKTYP, BKSTP
!                               FONCTION RETOURNE 0
!MODULE 
      INTEGER  QDFERR
      EXTERNAL QDFERR
!
!IMPLICITES
#include "codes.cdk"
#include "bpl.cdk"
#include "masques.cdk"
#include <ftnmacros.hf>
!
!*
      IF(BTYP .LT. -1) THEN
         MRBTYP = QDFERR('MRBTYP', 'VALEUR DE BTYP INVALIDE', WARNIN, ERBTYP)
         RETURN
      ENDIF
       
!     CONSTRUIRE LA CLEF BTYP
      IF(BTYP .EQ. -1) THEN
         MRBTYP = 0
         MRBTYP = IAND(BKSTP, BKSTPMSK)
         MRBTYP = IOR(MRBTYP, LSHIFT(IAND(BKTYP,BKTYPMSK),BPBKTYP))
         MRBTYP = IOR(MRBTYP, LSHIFT(IAND(BKNAT,BKNATMSK),BPBKNAT))
      ELSE 
         BKSTP  = IAND(BTYP,BKSTPMSK)
         BKTYP  = IAND(RSHIFT(BTYP, BPBKTYP), BKTYPMSK)
         BKNAT  = IAND(RSHIFT(BTYP, BPBKNAT), BKNATMSK)
         MRBTYP = 0
      ENDIF
  
      RETURN
      END
