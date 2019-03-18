!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r  verder - compute vertical derivative of input field
!                 with respect to field F_wlnph
!

!
      subroutine verder (F_der, F_infield, F_wlnph, F_con1, F_con2, &
                                   Minx,Maxx,Miny,Maxy, F_Nk,F_i0,F_in,F_j0,F_jn)
!
      use glb_ld
      implicit none
#include <arch_specific.hf>
!
      integer Minx,Maxx,Miny,Maxy, F_Nk,F_i0,F_in,F_j0,F_jn
      real    F_con1, F_con2
      real    F_der(Minx:Maxx,Miny:Maxy,F_Nk),F_infield (Minx:Maxx,Miny:Maxy,F_Nk), &
              F_wlnph(Minx:Maxx,Miny:Maxy,F_Nk)
!
!author
!     alain patoine - after ddsig (efr)
!
!revision
! v2_00 - Lee V.            - initial MPI version (from verder v1_03)
! v3_21 - Lee V.            - Output Optimization
! v4_12 - Lee V.            - OpenMP racing condition:need omp end parallel
!
!object
!     See above id
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_der        O    - derivative of the put field with respect to log of
!                     hydrostatic pressure
! F_infield    I    - input field on the eta levels of the model
! F_wlnph      I    - log of hydrostatic pressure
! F_con1       I    - used for boundary conditions
! F_con2       I    - used for boundary conditions
! F_i0         I    - starting point of calculation on W-E axis
! F_in         I    - ending point of calculation on W-E axis
! F_j0         I    - starting point of calculation on N-S axis
! F_jn         I    - ending point of calculation on N-S axis
!

!
!*
      integer i, j, k
!

!$omp parallel private(k,i)
!$omp do
      do j=F_j0,F_jn
      do k=2,F_nk
      do i=F_i0,F_in
         F_der(i,j,k) =   ( F_infield (i,j,k) - F_infield (i,j,k-1) ) &
                        / ( F_wlnph(i,j,k) - F_wlnph(i,j,k-1) )
      end do
      end do
      end do
!$omp enddo
!$omp end parallel
!
      do j=F_j0,F_jn
      do i=F_i0,F_in
         F_der(i,j,1) =  F_der(i,j,2)
      end do
      end do

!
!$omp parallel private(k,i)
!$omp do
      do j=F_j0,F_jn
      do k=2,F_nk-1
      do i=F_i0,F_in
         F_der(i,j,k)=((F_wlnph(i,j,k+1)-F_wlnph(i,j,k))   * F_der(i,j,k) &
                      +(F_wlnph(i,j,k)  -F_wlnph(i,j,k-1)) * F_der(i,j,k+1)) &
                      /(F_wlnph(i,j,k+1)-F_wlnph(i,j,k-1))
      end do
      end do
      end do
!$omp enddo
!$omp end parallel
!
!$omp parallel private(i)
!$omp do
      do j=F_j0,F_jn
      do i=F_i0,F_in
        F_der(i,j,1)    =   F_con1 * F_der(i,j,1) &
                          + (1.0 - F_con1) * F_der(i,j,2)
        F_der(i,j,F_nk) =   F_con2 * F_der(i,j,F_nk) &
                          + (1.0 - F_con2) * F_der(i,j,F_nk-1)
      end do
      end do
!$omp enddo
!$omp end parallel
!
      return
      end
