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

!**s/r itf_fft_sin - multiple fast real sine transform
!
      subroutine itf_fft_sin ( F_a, F_inc, F_jump, F_lot )
      use fft
      implicit none
#include <arch_specific.hf>

      integer F_inc, F_jump, F_lot
      real*8  F_a(*)
!
!author
!      Jean Cote      Sept 1999
!
!revision
! v4_50 - Desgagne M.   - update standards
!
!object
!     real sine transform of length n-1
!     the transform is its own inverse
!     self inverting implementation of Numerical Recipes pp. 508-512
!     a     is the array containing input and output data
!     inc   is the increment within each data 'vector'
!          (e.g. inc=1 for consecutively stored data)
!     jump  is the increment between the start of each data vector
!     lot   is the number of data vectors
!     definition of transform:
!     -------------------------
!     r(k) = sqrt(2/n)*sum(i=1,...,n-1)(a(i)*sin(i*k*pi/n))
!
!     Note for 'a' stored as a(n1,n2) then
!
!        for a transform along the first dimension
!
!           inc   = 1
!           jump  = n1
!
!        for a transform along the second dimension
!
!           inc   = n1
!           jump  = 1
!
!     The subroutine itf_fft_set must have been called to set-up
!     the commons in fft.cdk


      integer i, j, k, is, j0, jlot
      real*8  ai, as, ya, ys, xnor, w(511*Fft_nstore)

      real*8 zero, half, one, two, four
      parameter( zero = 0.0d0 )
      parameter( half = 0.5d0 )
      parameter( one  = 1.0d0 )
      parameter( two  = 2.0d0 )
      parameter( four = 4.0d0 )

      integer ija, ijw
      ija(i,j) = 1 + (j0+j-1)*F_jump + (i-1)*F_inc
      ijw(i,j) = j + i*511
!     __________________________________________________________________
!
!
      xnor = sqrt( half * Fft_n )

      do 100 j0=0,F_lot-1,511
         jlot  = min( 511, F_lot - j0 )
         do j=1,jlot
            w( ijw(0,j) ) = zero
         end do

         do i=1,Fft_n-Fft_m-1
            is = Fft_n - i
            do j=1,jlot
               ai= F_a( ija(i ,j) )
               as= F_a( ija(is,j) )
               ya=  ( as - ai ) * xnor
               ys=  two * Fft_ssin( i ) * ( as + ai ) * xnor
               w( ijw(i ,j) ) = ys + ya
               w( ijw(is,j) ) = ys - ya
            end do
         end do

         if ( Fft_n == 2 * Fft_m ) then
            do j=1,jlot
               w( ijw(Fft_m,j) ) = four * F_a( ija(Fft_m,j) ) * xnor
            end do
         end if

         call ffft8( w, 511, 1, jlot, -1 )

         do j=1,jlot
            F_a( ija(1,j) ) = half * w( ijw(0,j) )
         end do

         do k = 2,Fft_n-2,2
            do j=1,jlot
               F_a( ija(k+1,j) ) = w( ijw(k,j)   ) + F_a( ija(k-1,j) )
               F_a( ija(k,  j) ) = w( ijw(k+1,j) )
            end do
         end do

         if ( Fft_n /= 2 * Fft_m ) then
            do j=1,jlot
               F_a( ija(Fft_n-1,j) ) = w( ijw(Fft_n,j) )
            end do
         end if

  100 continue
!     __________________________________________________________________
!
      return
      end
