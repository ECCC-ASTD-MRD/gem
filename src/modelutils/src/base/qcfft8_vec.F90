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

!     subroutine 'qcfft8_vec' - multiple fast real shifted cosine transform

!     real shifted cosine transform of length n
!     implementation inspired by Numerical Recipes pp. 513
!     but with a different normalization

!     created: Sept 99 by j.cote, rpn
!author McTaggart-Cowan Dec 2006 for vectorization
!
!revision
! v3_30 - Qaddouri/Lee      - initialize work vector to zero (LAM bug)
!
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
!     The subroutine SETSCQR must have been called to set-up
!     the commons COMFFT8 and COMFFT8X
!
!-----------------------------------------------------------------------

subroutine qcfft8_vec(a, inc, jump, lot, isign)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
   include "rmnlib_basics.inc"

   ! Input/output variables
   integer :: inc, jump, lot, isign
   real(REAL64) :: a(*)

   ! Comdec variable declarations
   integer :: n, m, nstore
   pointer ( ptss,ssin(n-m-1) ), ( ptcc,ccos(n-m-1) )
   pointer ( ptqs,qsin(0:m-1) )
   common    / comfft8x / ptss, ptcc, ptqs, n, m, nstore

   ! External functions
   real(REAL64) ::  ssin, ccos, qsin

   ! Internal parameters
   integer, parameter :: MAX_VEC=511
   real(REAL64), parameter :: zero=0.0, half=0.5, one=1.0,  &
        two=2.0, four=4.0

   ! Internal variables and parameters
   integer :: i, j, k, is, k1, kk, j0, jlot, ija, maxVec, mMax
   real(REAL64) :: ai, as, ya, ys, c, s, rr
   real(REAL64), dimension(0:nstore-1,MAX_VEC) :: w2,fftData
   logical :: nIs2m
   !-----------------------------------------------------------------------

   ! Statement functions
   ija(i,j) = 1 + (j0+j-1)*jump + i*inc

   ! Vector length setup
   if (n == 2*m) then
      nIs2m = .true.
   else
      nIs2m = .false.
   endif
   w2=zero

   ! Begin main loop for a large number of vectors (>MAX_VEC)
   LONGDO: do j0=0,lot-1,MAX_VEC
      jlot  = min( MAX_VEC, lot - j0 )

      ! Fill local data array
      do j=1,jlot
         do i=0,nstore-1
            fftData(i,j) = a(ija(i,j))
         enddo
      enddo

      if ( isign .eq. -1 ) then

         ! --- Transform from gridpoint to Fourier ---
         do j=1,jlot
!VDIR NOSYNC
            do i=0,m-1
               is = n-i-1
               ai = fftData(i,j)
               as = fftData(is,j)
               ys = ai + as
               ya = two * qsin(i) * (as - ai)
               w2(i,j) = ys + ya
               w2(n-i-1,j) = ys - ya
            enddo
         enddo
         if (.not.nIs2m) w2(m,1:jlot) = two * fftData(m,1:jlot)

         ! Call external FFT calculator
         maxVec = MAX_VEC
         call ffft8( w2(0:nstore-1,1:jlot), 1, nstore, jlot, -1 )

         ! Prepare results for output
         fftData(0,1:jlot) = half * w2(0,1:jlot)
         if (nIs2m) then
            mMax = m-1
         else
            mMax = m
         endif
         do j=1,jlot
            do k=1,mMax
               kk = 2*k
               k1 = kk+1
               c = ccos(k)
               s = ssin(k)
               fftData(kk-1,j) = -s*w2(kk,j) + c*w2(k1,j)
               fftData(kk,j) = c*w2(kk,j) + s*w2(k1,j)
            enddo
         enddo
         if (nIs2m) then
!VDIR NOSYNC
            do j=1,jlot
               fftData(2*m-1,j) = -w2(2*m,j)
               fftData(n-1,j) = fftData(n-1,j) * half
            enddo
         endif
         do j=1,jlot
            do k=2*m-3,1,-2
               fftData(k,j) = fftData(k,j) + fftData(k+2,j)
            enddo
         enddo

      elseif ( isign .eq. +1 ) then

         ! --- Transform from Fourier to gridpoint ---
         do j=1,jlot
            w2(0,j) = fftData(0,j)
            w2(1,j) = zero
         enddo
         do j=1,jlot
            do k=2,n-1,2
               w2(k,j) = fftData(k,j) * half
            enddo
            do k=3,2*m-1,2
               w2(k,j) = (fftData(k-2,j) - fftData(k,j) ) * half
            enddo
         enddo
         if (nIs2m)  then
            c = one
         else
            c = half
         endif
         w2(2*m+1,1:jlot) = fftData(2*m-1,1:jlot) * c
         do j=1,jlot
            do k=1,m
               if ( k .lt. m .or. n .ne. 2 * m ) then
                  c = ccos(k)
                  s = ssin(k)
               else
                  c = zero
                  s = one
               endif
               kk = 2*k
               k1 = kk + 1
               rr = w2(kk,j)
               w2(kk,j) = c * rr - s * w2(k1,j)
               w2(k1,j) = s * rr + c * w2(k1,j)
            enddo
         enddo

         ! Call external FFT calculator
         maxVec = MAX_VEC
         call ffft8( w2(0:nstore-1,1:jlot), 1, nstore, jlot, +1 )

         ! Prepare inverted results for output
         do j=1,jlot
!VDIR NOSYNC
            do i=0,m-1
               is = n-i-1
               ys = ( w2(i,j) + w2(is,j) ) * half
               ya = ( w2(is,j) - w2(i,j) ) / ( four * qsin( i ) )
               fftData(i,j) = ys + ya
               fftData(is,j) = ys - ya
            enddo
         enddo
         if ( .not.nIs2m )  then
            do j=1,jlot
               fftData(m,j) = w2(m,j)
            enddo
         endif

         ! End of main branch for fwd/reverse transformations
      endif

      ! Empty local array into output vector
!VDIR NOSYNC
      do j=1,jlot
         do i=0,nstore-1
            a(ija(i,j)) = fftData(i,j)
         enddo
      enddo

      ! End of main loop for number of vectors
   enddo LONGDO
   !-----------------------------------------------------------------------
   return
end subroutine qcfft8_vec
