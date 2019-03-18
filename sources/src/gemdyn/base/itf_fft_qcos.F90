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
!
!**s/r itf_fft_qcos - multiple fast real shifted cosine transform
!
      subroutine itf_fft_qcos ( F_a, F_inc, F_jump, F_lot, F_isign )
      use fft
      implicit none
#include <arch_specific.hf>
!
      integer F_inc, F_jump, F_lot, F_isign
      real*8  F_a(*)
!
!author
!      Jean Cote      Sept 1999
!
!revision
! v3_30 - McTaggart-Cowan   - call qcfft8_vec if defined (NEC)
! v3_30 - Qaddouri/Lee      - initialize work vector to zero (LAM bug)
! v4_50 - Desgagne M.       - update standards
!
!object
!     real shifted cosine transform of length n
!     implementation inspired by Numerical Recipes pp. 513
!     but with a different normalization
!     a     is the array containing input and output data
!     inc   is the increment within each data 'vector'
!          (e.g. inc=1 for consecutively stored data)
!     jump  is the increment between the start of each data vector
!     lot   is the number of data vectors
!     isign = +1 for transform from fourier to gridpoint
!           = -1 for transform from gridpoint to fourier
!     definition of transform:
!     -------------------------
!     isign=+1: r(i) = sum(k=0,...,n-1)(a(k)*cos((i+1/2)*k*pi/n))
!     isign=-1: r(k) = sum(i=0,...,n-1)(a(i)*cos((i+1/2)*k*pi/n))
!                      * ((2-delta(k,0))/n)
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
!


      integer i, j, k, is, k1, kk, j0, jlot
      real*8  ai, as, ya, ys, c, s, rr, w(511*Fft_nstore)
      real*8 zero, half, one, two, four
      parameter( zero=0.0, half=0.5, one=1.0, two=2.0, four=4.0 )

      integer ija, ijw
      ija(i,j) = 1 + (j0+j-1)*F_jump + i*F_inc
      ijw(i,j) = j + i*511
!
!-----------------------------------------------------------------------
!
      w(:)=0.0
!
      do 100 j0=0,F_lot-1,511
         jlot  = min( 511, F_lot - j0 )

         if ( F_isign == -1 ) then
!
!     transform from gridpoint to Fourier
!
            do i = 0 , Fft_m-1
               is = Fft_n - i - 1
               do j=1,jlot
                  ai = F_a( ija(i ,j) )
                  as = F_a( ija(is,j) )
                  ys = ai + as
                  ya =  two *  Fft_qsin( i ) * ( as - ai )
                  w( ijw(i ,j) ) = ys + ya
                  w( ijw(is,j) ) = ys - ya
               end do
            end do

            if ( Fft_n /= 2 * Fft_m )  then
               do j=1,jlot
                  w( ijw(Fft_m,j) ) = two * F_a( ija(Fft_m,j) )
               end do
            end if

            call ffft8( w, 511, 1, jlot, -1 )

            do j=1,jlot
               F_a( ija(0,j) ) = w( ijw(0,j) ) * half
            end do

            do k = 1 , Fft_m
               kk = 2*k
               k1 = kk + 1
               if ( k < Fft_m .or. Fft_n /= 2 * Fft_m ) then
                  c = Fft_ccos( k )
                  s = Fft_ssin( k )
                  do j=1,jlot
                     F_a(ija(kk-1,j)) = -s*w( ijw(kk,j) )+c*w( ijw(k1,j) )
                     F_a(ija(kk  ,j)) =  c*w( ijw(kk,j) )+s*w( ijw(k1,j) )
                  end do
               else
                  do j=1,jlot
                     F_a( ija(kk-1,j) ) = -w( ijw(kk,j) )
                  end do
               end if

            end do
            if ( Fft_n == 2 * Fft_m )  then
               do j=1,jlot
                  F_a( ija(Fft_n-1,j) ) = F_a( ija(Fft_n-1,j) ) * half
               end do
            end if

            do k = 2*Fft_m-3 , 1 , -2
            do j = 1,jlot
               F_a( ija(k,j) ) = F_a( ija(k,j) ) + F_a( ija(k+2,j) )
            end do
            end do
!
         elseif ( F_isign == +1 ) then
!
!     transform from Fourier to gridpoint
!
            do j=1,jlot
               w( ijw(0,j) ) = F_a( ija(0,j) )
               w( ijw(1,j) ) = zero
            end do


            do k = 2 , Fft_n-1 , 2
            do j = 1 , jlot
               w( ijw(k,j) ) = F_a( ija(k,j) ) * half
            end do
            end do

            do k = 3 , 2*Fft_m-1 , 2
            do j = 1 , jlot
              w(ijw(k,j)) = ( F_a( ija(k-2,j) ) - F_a( ija(k,j) ) ) * half
            end do
            end do
            if ( Fft_n == 2 * Fft_m )  then
               c = one
            else
               c = half
            end if
            do j=1,jlot
               w( ijw(2*Fft_m+1,j) ) = F_a( ija(2*Fft_m-1,j) ) * c
            end do

            do k = 1 , Fft_m

               if ( k < Fft_m .or. Fft_n /= 2 * Fft_m ) then
                  c = Fft_ccos( k )
                  s = Fft_ssin( k )
               else
                  c = zero
                  s = one
               end if

               kk = 2*k
               k1 = kk + 1

               do j=1,jlot
                  rr = w( ijw(kk,j) )
                  w( ijw(kk,j) ) =  c * rr - s * w( ijw(k1,j) )
                  w( ijw(k1,j) ) =  s * rr + c * w( ijw(k1,j) )
               end do

            end do

            call ffft8( w, 511, 1, jlot, +1 )

            do i = 0 , Fft_m-1

               is = Fft_n - i - 1

               do j=1,jlot
                  ys = ( w( ijw(i ,j) ) + w( ijw(is,j) ) ) * half
                  ya = ( w( ijw(is,j) ) - w( ijw(i ,j) ) ) / &
                       ( four * Fft_qsin( i ) )
                  F_a( ija(i ,j) ) = ys + ya
                  F_a( ija(is,j) ) = ys - ya

               end do

            end do
            if ( Fft_n /= 2 * Fft_m )  then
               do j=1,jlot
                  F_a( ija(Fft_m,j) ) = w( ijw(Fft_m,j) )
               end do
            end if

         end if

  100 continue
!
!-----------------------------------------------------------------------
!
      return
      end
