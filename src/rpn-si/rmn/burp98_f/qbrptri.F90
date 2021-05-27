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

!.S QBRPTRI
!**S/P QBRPTRI - TRIER LE TABLEAU D'ELEMENTS LUS EN ORDRE CROISSANT
      SUBROUTINE QBRPTRI(TABLEAU, NI, NJ)
      IMPLICIT   NONE
      INTEGER    NI, NJ,  TABLEAU(NI, NJ)
!
!AUTEUR: J. CAVEEN   MARS  1991
!REV 001 Y. BOURASSA MARS  1995 RATFOR @ FTN77
!
!OBJET( QBRPTRI )
!     SOUS-PROGRAMME SERVANT A TRIER UN TABLEAU A DEUX DIMENSIONS
!     PAR ORDRE CROISSANT DU PREMIER ELEMENT I.E: TABLEAU(1,J)
!     ON UTILISE LA METHODE DE TRI DE SHELL
!
! ARGUMENTS
!
!     TABLEAU       ENT-SRT      TABLEAU A TRIER
!     NI            ENTREE       PREMIERE DIMENSION DE TABLEAU
!     NJ            ENTREE       DEUXIEME DIMENSION DE TABLEAU
!      
      INTEGER I, J, K, M, NUM
! 
!*

      M = NJ
 10   IF(M .GT. 1) THEN
         M = (M+2)/3
         DO 40 I = M+1,NJ
            DO 30 J = I,M+1,-M
               IF(TABLEAU(1,J-M) .LT. TABLEAU(1,J)) GO TO 40
              
               DO 20 K = 1,NI
                  NUM            = TABLEAU(K,J)
                  TABLEAU(K,J)   = TABLEAU(K,J-M)
                  TABLEAU(K,J-M) = NUM
 20               CONTINUE

 30            CONTINUE

 40         CONTINUE
         GO TO 10
      ENDIF

      RETURN
      END
