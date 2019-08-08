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
!**S/R ez_CARTAUV - compute the components of the winds in the rotated system of
!                coordinates from the winds in the rotated cartesian space
!
      SUBROUTINE ez_CARTAUV(U, V, UVCART, LON, LAT, NI, NJ)
      implicit none
      INTEGER NI, NJ 
      REAL    UVCART(3,NI*NJ), U(NI,NJ), V(NI,NJ), LON(NI,nj), LAT(ni,NJ)
!
!author michel roch - april 90
!
!arguments
!    out    U       - rotated component  U of the wind
!           V       - rotated component  V of the wind
!    in     xyz     - rotated winds in cartesian space
!           LON     - longitudes of the grid in the rotated system of coordinates
!           LAT     - latitudes of the grid in the rotated system of coordinates
!           NI      - E-W DIMENSION of the grid
!           NJ      - N-S dimension of the grid
!
!*

      INTEGER I, J, K 
      REAL*8    A, B, C, D, E, F, DAR

      DAR = ACOS(-1.)/180.
      K   = 0
 
      DO 20 J=1,NJ
         DO 10 I=1,NI
            K      = K+1
            A      = COS(DAR*LON(I,j))
            B      = SIN(DAR*LON(I,j))
            E      = COS(DAR*LAT(i,J))
            F      = SIN(DAR*LAT(i,J))
            U(I,J) = (UVCART(2,K)*A) - (UVCART(1,K)*B)
            C      = (UVCART(1,K)*A) + (UVCART(2,K)*B)
            D      = SQRT(C**2 + UVCART(3,K)**2 )
            V(I,J) = SIGN(D, (UVCART(3,K)*E)-(C*F))
 10         CONTINUE
 20      CONTINUE

      RETURN
      END 
