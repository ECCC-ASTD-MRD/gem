!/* EDITFST - Collection of useful routines in C and FORTRAN
! * Copyright (C) 1975-2014  Environnement Canada
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
!** S/R SETPER ETABLIR UNE DATE (DEBUT/FIN DE PERIODE)
      SUBROUTINE SETPER(DN, ECART, DUREE, DELTA)
      use configuration
      IMPLICIT NONE 
  
      INTEGER    DN(*), ECART, DUREE, DELTA
!
!ARGUMENTS
!  ENTRE    - DN    - A) 'OPRUN'
!                     C) 'YYJJJZZ'
!                     D) -(CMCSTAMP)
!    "      - ECART - NOMBRE D'HEURES PAR LEQUEL IF FAUT MODIFIER
!                     LA DATE QUI VIENT DU FICHIER "DN" (SI PRESENT)
!    "      - DUREE - DUREE DE LA PERIODE (SI PRESENT)
!    "      - DELTA - INTERVALE EN HEURES ENTRE LES CAS (SI PRESENT)
!
!AUTEURS
!         Y BOURASSA NOV  90
!REV 001  "     "    JUIL 91 ACCEPTE LES DATESTAMP CMC *
!REV 002  "     "    FEV  92 APPEL A LOW2UP AVANT IOPDATM
!Rev 003  M. Lepine  Nov  05 Remplacement de fstabt par qqexit
!Rev 004  M. Valin   Mai 14 Remplacement des comdecks par un module
!
!LANGUAGE   - FTN77 
!
!MODULES  
      EXTERNAL      IOPDATM, ARGDIMS, INCDAT, DATMGP, qqexit, LOW2UP
!
!*
      INTEGER       IOPDATM, ARGDIMS, DTG(14), I, K
      EQUIVALENCE   (K, DTG(14))
      CHARACTER*128 C

         print *,'PERIODE',jours
!     ETABLIR LE DATESTAMP DU CAS OU DU DEBUT DE LA PERIODE 
      IF(DN(1) .LT. 0) THEN   ! CMC date-time stamp numerique
         K = -DN(1)
      ELSE                    ! chaine de caracteres
         WRITE(C,LIN128) (DN(I), I=1,ARGDIMS(1))
         CALL LOW2UP(C, C)
         K = IOPDATM( C )     ! aller chercher le date-time stamp numerique associe
         IF(K .EQ. 10101011) THEN
            PRINT*,'DATE DE IOPDATM INACCEPTABLE'
            CALL qqexit(67) 
         ENDIF
      ENDIF
  
!     SI ON DOIT MODIFIER LA DATE QUI A ETE PASSEE
      IF(NP .GT. 1) CALL INCDAT(K, K, ECART)
      JOURS(1) = K
      JOURS(2) = K           ! date1 @ date1
      JOURS(3) = 0           ! DELTA 0 (pas de DELTA)
      JOURS(4) = 1           ! periode valide presente si jours(4) .ne. 0
  
!     SI ON A SEULEMENT DONNE UNE DATE (pas de duree ni de delta)
      IF(NP .LT. 3) THEN   ! plus rien a faire
         IF( DEBUG ) THEN
            CALL DATMGP( DTG )
            WRITE(6,600) (DTG(I),I=7,13), JOURS(1)
  600       FORMAT(' * CAS DU ',7A4,'  JULH = ',I8)
         ENDIF
      ELSE   !  NP >= 3,  on a @ date2 et peut-etre un delta
         JOURS(4) = -1 
         IF(DUREE .GE. 0) THEN
            call INCDAT(JOURS(2), JOURS(1), DUREE)  ! JOURS(2) = JOURS(1) + DUREE
         ELSE   ! les dates sont a l'envers, DUREE < 0
            JOURS(2) = JOURS(1)
            call INCDAT(JOURS(1), JOURS(2), DUREE)  ! JOURS(1) = JOURS(2) + DUREE
         ENDIF
         IF(NP .GT. 3) THEN
            JOURS(3) = ABS( DELTA ) 
         ELSE
            JOURS(3) = 0                      ! pas de DELTA
         ENDIF
         IF( DEBUG ) THEN
            CALL DATMGP( DTG )
            WRITE(6,601) (DTG(I),I=7,13)
  601       FORMAT(' * DEBUT *  ', 7A4) 
            CALL INCDAT(K, K, DUREE)
            CALL DATMGP( DTG )
            WRITE(6,602) (DTG(I),I=7,13)
  602       FORMAT(' *  FIN  *  ',7A4)
            WRITE(6,*) JOURS(1), ' @ ', JOURS(2), ' DELTA ', JOURS(3) 
         ENDIF
      ENDIF   ! NP 
  
      RETURN
      END 
