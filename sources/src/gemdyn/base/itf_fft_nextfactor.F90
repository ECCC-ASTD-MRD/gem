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

!**s/r itf_fft_nextfactor - calcul du prochain entier > 8 qui se factorise en 2, 3, et 5
!
      subroutine itf_fft_nextfactor2 ( F_n, F_ndown )
      implicit none
#include <arch_specific.hf>
!
      integer F_n, F_ndown
!
!author
!     Jean Cote - 1990
!
!revision
! v2_00 - Lee V.            - initial MPI version (from ngfft v1_03)
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
!  F_n         I/O   - en sortie le plus petit entier >= n qui factorise
!
!
      integer n
      parameter ( n = 3 )
      integer k( n ) , m, i, j, keep
      data m / 8 /
      data k / 2 , 3 , 5 /
!
!     ---------------------------------------------------------------
!
      if ( F_n <= m ) F_n = m + 1
      keep= F_n
! going up
      F_n = F_n - 1
  10  F_n = F_n + 1
      i = F_n

  20  do 30 j=1,n
         if( mod(i,k(j)) == 0 ) go to 40
  30  continue

      go to 10
  40  i = i/k(j)

      if( i /= 1 ) go to 20

! going down
      F_ndown = keep + 1
 100  F_ndown = F_ndown - 1
      i = F_ndown

 200  do 300 j=1,n
         if( mod(i,k(j)) == 0 ) go to 400
 300  continue

      go to 100
 400  i = i/k(j)

      if( i /= 1 ) go to 200
!
!     ---------------------------------------------------------------
!
      return
      end
