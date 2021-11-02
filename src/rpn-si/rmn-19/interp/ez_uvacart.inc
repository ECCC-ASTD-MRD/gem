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
!**S/R UVACART  - compute the winds in the cartesian space from
!                 the components

!
      SUBROUTINE EZ_UVACART( XYZ, U, V, LON, LAT, NI, NJ)
      implicit none
      INTEGER NI, NJ 
      REAL    U(NI,NJ), V(NI,NJ), XYZ(3,NI*NJ), LON(NI,NJ), LAT(NI,NJ)
!
!author michel roch - april 90
!
!arguments
!    out    xyz   - unrotated winds in cartesian space
!    IN     U     - unrotated component  U of the wind
!           V     - unrotated component  V of the wind
!           LON   - longitudes of the grid in the unrotated system of coordinates
!           LAT   - latitudes  of the grid in the unrotated system of coordinates
!           NI    - E-W DIMENSION of the grid
!           NJ    - N-S dimension of the grid
!
!*

      INTEGER I, J, K 
      REAL*8    A, B, C, D, DAR

      DAR = ACOS(-1.)/180.
      K   = 0
 
      DO 20 J=1,NJ
         DO 10 I=1,NI
            K        = K+1
            A        = SIN(DAR*LON(I,J))
            B        = COS(DAR*LON(I,J))
            C        = SIN(DAR*LAT(I,J))
            D        = COS(DAR*LAT(I,J))
            XYZ(1,K) = -(U(I,J)*A) - (V(I,J)*B*C)
            XYZ(2,K) =  (U(I,J)*B) - (V(I,J)*A*C)
            XYZ(3,K) =   V(I,J)*D
 10         CONTINUE
 20      CONTINUE

      RETURN
      END 
