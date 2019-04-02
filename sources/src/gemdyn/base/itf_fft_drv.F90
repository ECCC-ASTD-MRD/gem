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

!**s/r itf_fft_drv
!

subroutine itf_fft_drv( F_vec, stride, jump, num, direction)
   ! itf_fft_drv -- legacy-compatible FFT driver

   ! Perform (num) in-place Fourier transforms of the array F_vec,
   ! where elements of each transform are separated by (stride)
   ! REAL*8 units and adjacent arrays are separated by a stride of
   ! (jump).

   ! The transform type is given by the global variable Fft_type_S
   ! in module fft, and the transform direction is specified as a
   ! parameter.

   use fftw3 ! Incorporate the FFTW module
   use gem_fft, only : Fft_type_S, Fft_n, fftw_thread_initialized
   implicit none

   integer, intent(in) :: stride, jump, num, direction

   ! The array is passed as a deferred-size array, so there is no guarantee of
   ! alignment.  In fact, in common use this driver often sees a subset of another
   ! array.
   real*8, intent(inout) ::  F_vec(*)
   real*8 :: norm_factor

   integer :: ii, jj ! Loop variables for normalization

   ! Variables necessary for FFTW
   type(C_PTR) my_plan
   integer fftw_kind(1)
   integer fftw_n(1)

   if (Fft_type_S == 'SIN') then
      ! DST selected, periodic about the points
      ! 0 and Fft_n, which are outside the array
      ! and implicitly zero.  This transform is
      ! its own inverse.
      fftw_kind(1) = FFTW_RODFT00
      norm_factor = 1.0/sqrt(2.0d0*Fft_n)
      fftw_n(1) = Fft_n-1
   elseif (Fft_type_S == 'QCOS' .and. direction == -1) then
      ! DCT selected.  The DCT implied is even about
      ! points 0.5 and Fft_n+0.5, which lie between
      ! the first (1) and last (Fft_n) elements of
      ! the array and the first out-of-bound points.
      fftw_kind(1) = FFTW_REDFT10
      norm_factor = 1.0/(1.0d0*Fft_n)
      fftw_n(1) = Fft_n
   elseif (Fft_type_S == 'QCOS' .and. direction == 1) then
      ! iDCT selected.  This transform has an implied
      ! even boundary condition about point 0.5 and
      ! an implied odd boundary condition about (Fft_n+0.5),
      ! which lie one half-point outside of the supplied
      ! domain.  Additionally, because fftw and the legacy
      ! library treat the k=0 wavenumber differently,
      ! we need a pre-normalization step to double the
      ! term.
      fftw_kind(1) = FFTW_REDFT01
      norm_factor = 0.5
      fftw_n(1) = Fft_n
      do jj=1,num
         F_vec(1+(jj-1)*jump) = F_vec(1+(jj-1)*jump)*2
      end do
   elseif (Fft_type_S == 'PERIODIC') then
      ! Perform a real-to-complex periodic transform.  The native
      ! fftw routines are NOT drop-in replacements here.

      ! FFTW transforms a real input array into a complex output
      ! array, such that the real and imaginary components of each
      ! element are contiguous in memory, whereas the librmn routines
      ! operate exclusively on a REAL*8 array, and the real and
      ! imaginary components are separated by (stride).

      ! These cases are handled by a separate wrapper function, since
      ! temporary storage is necessary.
      call dft_wrapper(F_vec,stride,jump,num,direction)
      return
   else
      ! Do not use gem_error; it could cause a deadlock if OpenMP is used and only
      ! one thread encounters an error
      stop 'Invalid transform type or direction'
   end if

   ! Check for FFTW thread initialization, and initialize if necessary
   if (.not. fftw_thread_initialized) then
!$omp critical(fftw_lock)
      ! Re-check the initialization flag -- it's possible for a thread
      ! to have waited for the lock while another thread was initializing
      ! FFTW, and we don't need to repeat the init_threads call.  The
      ! omp critical section implies synchronization.
      if (.not. fftw_thread_initialized) then
         ii = fftw_init_threads()
         if (ii .eq. 0) then
            ! The caution about deadlocks is even stronger here.  Since this
            ! is in a critical section, only one thread will ever reach
            ! this point in the event of an error.
            stop 'Error in FFTW thread initialization'
         end if
         fftw_thread_initialized = .true.
      end if
!$omp end critical(fftw_lock)
   end if

   ! Perform a real transform, using the fftw_plan_many_r2r interface.
   ! FFTW plan creation is not thread-safe, but ESTIMATE should be fast
!$omp critical(fftw_lock)
   call fftw_plan_with_nthreads(1)
   my_plan = fftw_plan_many_r2r(1, fftw_n, num, & ! Rank 1, transform size, number of transforms
                                F_vec, fftw_n, stride, jump, & ! Input array, physical size, element separation, array separation
                                F_vec, fftw_n, stride, jump, & ! Same for output
                                fftw_kind, FFTW_ESTIMATE)
!$omp end critical(fftw_lock)

   if (.not. c_associated(my_plan) .and. num > 0) then
      ! Do not use gem_error; it could cause a deadlock if OpenMP is used and only
      ! one thread encounters an error
      stop 'FFTW returned a null plan'
   else
      ! FFTW plan execution _is_ thread-safe
      call fftw_execute_r2r(my_plan,F_vec,F_vec)
   end if

   ! Each transform has a built-in normalization factor to apply
   do ii=0,fftw_n(1)
      do jj=1,num
         F_vec(1+(jj-1)*jump + ii*stride) = F_vec(1+(jj-1)*jump+ii*stride)*norm_factor
      end do
   end do

   ! The QCOS transform has special treatment of the DC (0 wavenumber)
   ! term: the librmn routine calculates a term of half the magnitude
   ! compared to FFTW, so this needs correction in a separate pass.
   if (Fft_type_S == 'QCOS' .and. direction == -1) then
      do jj=1,num
         F_vec(1+(jj-1)*jump) = F_vec(1+(jj-1)*jump)/2
      end do
   end if

   ! FFTW plan destruction is again not thread-safe
!$omp critical(fftw_lock)
   if (c_associated(my_plan)) then
      call fftw_destroy_plan(my_plan)
   end if
!$omp end critical(fftw_lock)

end subroutine

subroutine dft_wrapper(F_vec,stride,jump,num,direction)
   ! dft_wrapper -- perform an "in-place" real-to-complex fft
   ! (or complex-to-real ifft), using the FFTW library
   ! to perform the transform itself

   ! Code using the fft in rmnlib treats complex numbers generated
   ! by the FFT as a pair of real values.  However, unlike the
   ! complex datatype that places these values at adjacent memory
   ! memory locations, rmnlib separates real and complex values
   ! along the array dimension.  That is, for a transform along
   ! dimension 3, F(1,1,1) is a real value, then F(1,1,2) is an
   ! imaginary value.  Modeling this inside FFTW demands either
   ! a temporary array or the use of the 'guru' planner routine,
   ! which allows separate arrays for the real and imaginary
   ! components.  This routine uses the latter.

   use fftw3 ! FFTW routines
   use gem_fft, only : Fft_n
   use iso_c_binding
   implicit none

   ! Paramters
   integer, intent(in) :: stride, &   ! Logical distance between adjacent f_vec elements
                          jump,   &   ! Logical distance between adjacent f_vec arrays (usually 1)
                          num,    &   ! Number of transforms to perform
                          direction   ! -1 for real-to-spectral, 1 for spectral-to-real
   real*8, intent(inout) ::  F_vec(*) ! Input/output array

   ! Local variables
   type(fftw_iodim), dimension(1) :: dims, &    ! Stride and number information for the transform
                                     batch_dims ! Corresponding information for the batching
   type(C_PTR) :: my_plan ! FFTW plan
   integer :: ii, jj ! Loop variables, for normalization

   if (direction == -1) then ! Real-to-spectral transform
      ! Stride information
      dims(1)%n = Fft_n
      dims(1)%is = stride ! Input real-values are separated by (stride) units
      dims(1)%os = 2*stride ! Output real-values have twice the separation

      ! Batching information
      batch_dims(1)%n = num
      batch_dims(1)%is = jump
      batch_dims(1)%os = jump

      ! FFTW plan creation is not thread-safe
!$omp critical(fftw_lock)
      call fftw_plan_with_nthreads(1)
      my_plan = fftw_plan_guru_split_dft_r2c(1,dims,1,batch_dims, & ! Dimensioning information
                                             F_vec,F_vec,F_vec(1+stride), & ! Input, output real, output imag
                                             FFTW_ESTIMATE+FFTW_PRESERVE_INPUT) ! FFTW parameters
!$omp end critical(fftw_lock)

      if (.not. c_associated(my_plan) .and. num > 0) then
         ! Do not use gem_error; it could cause a deadlock if OpenMP is used and only
         ! one thread encounters an error
         stop 'FFTW returned a null plan'
      else
         ! Execute the plan, which is thread-safe
         call fftw_execute_split_dft_r2c(my_plan,F_vec,F_vec,F_vec(1+stride))
      end if

      ! Normalize the output
      do ii=0,Fft_n
         do jj=1,num
            F_vec(1+(jj-1)*jump + ii*stride) = F_vec(1+(jj-1)*jump+ii*stride)/Fft_n
         end do
      end do

!$omp critical(fftw_lock)
      ! Clean up any associated memory
      call fftw_destroy_plan(my_plan)
!$omp end critical(fftw_lock)

   elseif (direction == 1) then ! Spectral-to-real transform
      ! Stride information
      dims(1)%n = Fft_n
      dims(1)%is = 2*stride ! real-to-real (and imag-to-imag) distance is 2*stride on input
      dims(1)%os = stride ! Output real-to-real distance is 1*stride

      ! Batching information
      batch_dims(1)%n = num
      batch_dims(1)%is = jump
      batch_dims(1)%os = jump

!$omp critical(fftw_lock)
      call fftw_plan_with_nthreads(1)
      my_plan = fftw_plan_guru_split_dft_c2r(1,dims,1,batch_dims, & ! Dimensioning information
                                             F_vec,F_vec(1+stride),F_vec, & ! Input real, input imag, output
                                             FFTW_ESTIMATE+FFTW_PRESERVE_INPUT) ! FFTW parameters
!$omp end critical(fftw_lock)

      if (.not. c_associated(my_plan) .and. num > 0) then
         stop 'FFTW returned a null plan'
      else
         call fftw_execute_split_dft_c2r(my_plan,F_vec,F_vec(1+stride),F_vec)
      end if

!$omp critical(fftw_lock)
      call fftw_destroy_plan(my_plan)
!$omp end critical(fftw_lock)
   else
      stop 'Invalid periodic transform direction'
   end if

end subroutine
