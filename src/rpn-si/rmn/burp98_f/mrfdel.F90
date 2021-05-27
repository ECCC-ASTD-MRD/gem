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
!.S MRFDEL
!**S/P MRFDEL - EFFACER UN ENREGISTREMENT D'UN FICHIER BURP
!
      FUNCTION MRFDEL( HANDLE )
      IMPLICIT NONE
      INTEGER  MRFDEL, HANDLE
!
!AUTEUR  J. CAVEEN   NOVEMBRE 1990
!REV 001 Y. BOURASSA MARS     1995 RATFOR @ FTN77
!rev 002 j. caveen   sept     1995 ajout d' appel a mrfprm pour
!                                  produire un message plus explicite
!REV 003 M. Lepine   sept    1997 nouveau format de date AAAAMMJJ (an 2000)
!
!OBJET( MRFDEL )
!     EFFACER L'ENREGISTREMENT CONTENU DANS UN FICHIER BURP DONT
!     LE POINTEUR EST HANDLE.
!                                                                       
!ARGUMENT
!     HANDLE  ENTREE  POINTEUR A L'ENREGISTREMENT A EFFACER
!                                                         
!IMPLICITES
#include "defi.cdk"
#include "burpopt.cdk"
#include "codes.cdk"
!
!MODULE 
      EXTERNAL XDFDEL, mrfprm
      INTEGER  XDFDEL,mrfprm
!
!*

      integer IOUT
      DATA     IOUT /6/
!
      CHARACTER*9  ISTNID
      INTEGER IIDTYP, ILAT, ILON, IDATE, ITEMPS, ISUP(1), INSUP
      INTEGER IDX,IDY,IFLGS,ILNGR,irien

      MRFDEL = -1

!
!     Obtenir les parametres descripteurs de l'enregistrement
!
      IF (MESSNIV .LE. INFORM) THEN
         IRIEN = handle
         INSUP = 0
         IRIEN = MRFPRM(IRIEN,ISTNID, IIDTYP, ILAT, ILON, IDX,IDY, IDATE, ITEMPS,IFLGS, ISUP, INSUP,ILNGR)
      ENDIF

!     EFFACER L'ENREGISTREMENT
      MRFDEL = XDFDEL( HANDLE )

      IF(MESSNIV.LE.INFORM .AND. MRFDEL.GE.0) THEN
         WRITE(6,1100) ISTNID, IIDTYP, ILAT, ILON, IDX,IDY, IDATE, ITEMPS,IFLGS,ILNGR
      ENDIF

 1100 FORMAT(' MRFDEL- EFFACE - STNID=',A9,' IDTYP=',I3, ' LAT=',I5,' LON=',I5,' DX=',i4,' DY=',i4,' DATE=',I8, ' TEMPS=',I4,' FLGS=',i8,' LNGR=',i6)

      RETURN
      END
