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

! Update: Dec 2018 -- Christopher Subich

! With the introduction of the fftw library as the "back end" for the Fourier
! Transforms (and related DCT/DST transforms), we have need of a new interface.
! "Legacy" routines itf_fft_set and itf_fft_drv still exist outside this module
! to provide for strict backwards compatibility, but the driver routines located
! here take better advantage of the FFTW library and should be used for all new
! or refactored code.

! The core conceit of FFTW is that it splits the transform into a planning and
! execution stage.  By optimizing (and timing) various planning approaches, the
! library purports to achieve faster execution than via a naive algorithm.
! Implicitly, the library assumes that transforms will be called several times
! on the same set of arrays, justifying the time spent in the planning stages.

! This splitting approach has been adopted (although somewhat less strictly) by
! other libraries like the MKL, so we present an abstracted form of it here to
! allow for future portability.

! To allow for the greatest generality, this module implements transforms for an
! entire (logical) three-dimensional array at once.   That is, one call will
! transform (for example) f(:,:,:) along the third dimension.  Thanks to the
! FFTW 'guru' interface, this is possible with a single plan.

! One lingering question with trigonometric transforms is normalization.  This
! module takes a strictly numerical view, and it leaves normalization to the
! caller.  We provide a get_norm_factor function to provide a multiplicative
! normalization constant to make a forward->backward transform pair an identity,
! but that's all.

module gem_fft
   use iso_c_binding
   implicit none
   public
   save

   ! Necessary variables for the legacy interface
   character(len=12) :: Fft_type_S ! One of 'FOURIER', 'SIN', or 'QCOS'
   integer :: Fft_n

   ! FFTW variables
   logical :: fftw_thread_initialized = .false.

   ! Parameters to turn "-1" and "+1" arguments into descriptive names
   integer, parameter :: DFT_FORWARD = -1
   integer, parameter :: DFT_BACKWARD = +1

   ! DFT descriptor type, to hold bookkeeping information
   ! about a set of generated FFTW plans.  To move from FFTW
   ! to any other library while maintaining a compatible
   ! interface, this type would probably change.
   type dft_descriptor
      ! The size of the input/output arrays used when planning
      ! the transform; this provides sanity and bounds checking
      integer, dimension(3) :: tsize
      ! Opaque array containing the FFTW plans (type void * in C)
      type(C_PTR) :: plan
   end type dft_descriptor

contains

   function get_dft_norm_factor(length, f_type)
      ! Get the multiplicative normalization factor corresponding to
      ! a particular forward/transformation pair.  This factor must
      ! be applied once during a forward->backward transformation
      ! pair to make the combination an identity.

      implicit none
      real*8 :: get_dft_norm_factor ! Result
      integer, intent(in) :: length ! Transform length
      character(len=*), intent(in) :: f_type ! transform type

      if (f_type == 'SIN') then
         get_dft_norm_factor = 0.5d0/(length+1)
      elseif (f_type == 'QCOS') then
         get_dft_norm_factor = 0.5d0/length
      elseif (f_type == 'PERIODIC') then
         get_dft_norm_factor = 1.0d0/length
      else
         call gem_error(-1,'get_dft_norm_factor','Invalid transform type')
      end if

      return
   end function

   subroutine make_r2r_dft_plan(tdesc, F_in, F_out, tdim, f_type, direction)
      ! Create an FFTW plan corresponding to a DCT or DST, storing the plan
      ! and associated data in a descriptor type.

      ! Usage: F_in and F_out are three-dimensional assumed-shape arrays,
      ! corresponding to the logical input/output arrays.  The resulting
      ! plan corresponds to transforming the entire F_in array along
      ! the tdim dimension.  If F_in and F_out are the same array, the
      ! transform is in-place, otherwise the transform is out-of-place.
      ! If F_in and F_out are not distinct but overlap in some way,
      ! the behaviour is undefined and probably erroneous.

      ! When the region to transform is a subset of a larger allocated
      ! block, the caller should supply an array slice.  This subroutine
      ! uses Fortran intrinsics to infer the array (and transform) size
      ! and pointer math to infer the relevant strides.

      ! The input and output arrays must have compatible logical shapes,
      ! although their physical layouts (strides) might differ.

      ! This routine should be called from outside an OpenMP 'parallel'
      ! region because FFTW invokes its own '#pragma omp parallel for'
      ! to internally parallelize; calling that from within an existing
      ! parallel region results in only a single thread being used.
      ! However, OpenMP critical sections are established for protection
      ! regardless.

      use fftw3 ! Incorporate the FFTW module
      use iso_c_binding
      use omp_lib, only : omp_get_max_threads
      implicit none

      integer, intent(in) :: tdim, &   ! Dimension to transform, 1-3
                             direction ! DFT_FORWARD or _BACKWARD

      ! Input and output arrays.  These are intent(inout) because the
      ! FFTW planning process performs measurements, overwriting any
      ! data already present.
      real*8, intent(inout), target :: F_in(:,:,:), F_out(:,:,:)

      character(len=*), intent(in) :: f_type ! Transform type

      type(dft_descriptor), intent(inout) :: tdesc

      ! Local variables:
      integer :: ierr ! Error from FFTW initialization
      integer :: istrides(3) ! Input array strides
      integer :: ostrides(3) ! Output array strides

      ! Variables for FFTW planning

      ! This module uses the 'guru' FFTW interface.  In this interface,
      ! a transform is described by a length, kind, and stride (between
      ! elements), and the guru-planner will plan an array of such
      ! transforms.  Each array is also described by a length and stride,
      ! between the kth element of the (i,j,...) subarray.

      ! This neatly encapsulates the kind of transforms needed for the FFT
      ! solver.
      type(fftw_iodim), dimension(1) :: dims       ! Stride information of the transform
      type(fftw_iodim), dimension(2) :: batch_dims ! ... and for the batching
      integer :: fftw_kind(1) ! Type of transform

      ! Pointers for pointer math -- see the planning section for discussion
      real*8, dimension(:), pointer, contiguous :: in_ptr, out_ptr

      ! The FFTW library has several different planning options, which offer
      ! a planning time/execution time tradeoff.  More time speant measuring
      ! performance in the planning stage (ideally) translates to a faster-
      ! executing FFT plan.

      ! There are several possible options here:
      ! FFTW_ESTIMATE performs no advance planning and takes a generic ``best'' plan
      ! based on system architecture.  This option should guarantee bit-for-bit
      ! compatibility of results between different runs on the same system.

      ! FFTW_MEASURE is the documentation-recommended default, which times several
      ! possible plans to explicitly measure the fastest.  In preliminary testing,
      ! this provided a 10% speedup in CPU time for the transform itself.

      ! FFTW_PATIENT tries a wider variety of plans.  This is unlikely to provide
      ! a noticeable further speedup in most cases, but the FFTW documentation
      ! does suggest that this option will disable OpenMP parallelization for
      ! the inner loops if it does not improve performance.

      ! As of March 2019, we use FFTW_ESTIMATE to provide for exact reproducibility
      ! of results between runs.
      integer, parameter :: FFTW_PLAN_TYPE = FFTW_ESTIMATE


   ! Check for FFTW thread initialization, and initialize if necessary
      if (.not. fftw_thread_initialized) then
!$omp critical(fftw_lock)
      ! Re-check the initialization flag, since omp critical implies
      ! memory synchronization
         if (.not. fftw_thread_initialized) then
            ierr = fftw_init_threads()
            if (ierr .eq. 0) then
               call gem_error(-1,'make_r2r_dft_plan','Error in FFTW thread initialization')
            end if
            fftw_thread_initialized = .true.
         end if
!$omp end critical(fftw_lock)
      end if

      ! Check for compatible input and output arrays
      if (.not. all(shape(F_in) == shape(F_out))) then
         call gem_error(-1,'make_r2r_dft_plan','Incompatible array shapes')
      end if

      ! Compute array strides
      istrides = 1 ! Sensible default
      ostrides = 1

      ! Pointer math: the stride, in terms of REAL*8 units, corresponds
      ! to the difference in memory location along the relevant dimension
      ! (in bytes) divided by the size of a double (8 bytes).  We only
      ! compute the stride along dimensions where the array is a non-
      ! singleton, to ensure that all implied memory addresses are in-
      ! bounds.
      if (size(F_in,1) > 1) then
         istrides(1) = (loc(F_in(2,1,1)) - loc(F_in(1,1,1)))/8
         ostrides(1) = (loc(F_out(2,1,1)) - loc(F_out(1,1,1)))/8
      end if
      if (size(F_in,2) > 1) then
         istrides(2) = (loc(F_in(1,2,1)) - loc(F_in(1,1,1)))/8
         ostrides(2) = (loc(F_out(1,2,1)) - loc(F_out(1,1,1)))/8
      end if
      if (size(F_in,3) > 1) then
         istrides(3) = (loc(F_in(1,1,2)) - loc(F_in(1,1,1)))/8
         ostrides(3) = (loc(F_out(1,1,2)) - loc(F_out(1,1,1)))/8
      end if

      ! With this stride information, we can condense the information
      ! for FFTW:

      ! The transform length is the size of the array over the right dimension
      dims%n = size(F_in,tdim)
      ! And the strides also correspond
      dims%is = istrides(tdim)
      dims%os = ostrides(tdim)

      ! An equivalent data structure is used to batch transforms together
      if (tdim == 1) then
         batch_dims(1)%n = size(F_in,2)
         batch_dims(1)%is = istrides(2)
         batch_dims(1)%os = ostrides(2)
         batch_dims(2)%n = size(F_in,3)
         batch_dims(2)%is = istrides(3)
         batch_dims(2)%os = ostrides(3)
      elseif (tdim == 2) then
         batch_dims(1)%n = size(F_in,1)
         batch_dims(1)%is = istrides(1)
         batch_dims(1)%os = ostrides(1)
         batch_dims(2)%n = size(F_in,3)
         batch_dims(2)%is = istrides(3)
         batch_dims(2)%os = ostrides(3)
      elseif (tdim == 3) then
         batch_dims(1)%n = size(F_in,1)
         batch_dims(1)%is = istrides(1)
         batch_dims(1)%os = ostrides(1)
         batch_dims(2)%n = size(F_in,2)
         batch_dims(2)%is = istrides(2)
         batch_dims(2)%os = ostrides(2)
      else
         call gem_error(-1,'make_r2r_dft_plan','Invalid transform dimension')
      end if

      ! Select the FFTW transfrom type from the input type and direction
      if (f_type == 'SIN') then
         fftw_kind(1) = FFTW_RODFT00
      elseif (f_type == 'QCOS' .and. direction == DFT_BACKWARD) then
         fftw_kind(1) = FFTW_REDFT01
      elseif (f_type == 'QCOS' .and. direction == DFT_FORWARD) then
         fftw_kind(1) = FFTW_REDFT10
      else
         call gem_error(-1,'make_r2r_dft_plan','Unknown transform type or direction')
      end if

      ! Begin setting the transform descriptor
      tdesc%tsize = shape(F_in)

      if (size(F_in,tdim) == 0) then
         ! The array slice is empty along the transform dimension,
         ! so there is in fact no transform.
         tdesc%plan = C_NULL_PTR
         return
      end if

      ! The FFTW Fortran interface is very handy, but it's slightly too restrictive.
      ! The C binding simply needs the address of the first member of the input
      ! and output arrays (double *), &(F_in[1][1][1]) here, and that allows FFTW
      ! to operate with general element strides and inter-array distances.

      ! Here in Fortran, we have that same functionality through the assumed-shape
      ! arrays used as parameters, WITHOUT the restriction that the arrays be
      ! contiguous (to allow for slices).  However, the FFTW-Fortran binding
      ! takes dimension(*) arrays (assumed-size) arrays as parameters, which
      ! must be one-dimensional, congiguous pieces of memory.

      ! However, with the magic of iso_c_binding we can perform some pointer math
      ! to effectively implement the C address-of operator

      call c_f_pointer(c_loc(F_in(1,1,1)),in_ptr,[1])
      call c_f_pointer(c_loc(F_out(1,1,1)),out_ptr,[1])

      ! Of the FFTW functions, only plan execution is guaranteed to be thread-safe.
      ! This routine should generally be called outside of an OpenMP region, but
      ! to absolutely ensure proper operation this is also a critical section.

!$omp critical(fftw_lock)
      call fftw_plan_with_nthreads(omp_get_max_threads())
      tdesc%plan = fftw_plan_guru_r2r(1,dims, & ! Transform rank and length/stride info
                                      2,batch_dims, & ! Batching lengh/stride info
                                      in_ptr, out_ptr, & ! Input and output arrays
                                      fftw_kind, FFTW_PLAN_TYPE) ! Transform type and planning flags
!$omp end critical(fftw_lock)
      if (.not. c_associated(tdesc%plan)) then
         call gem_error(-1,'make_r2r_dft_plan','FFTW returned null plan')
      end if
   end subroutine

   subroutine execute_r2r_dft_plan(tdesc, F_in, F_out)
      ! Execute the provided DCT/DST plans.  Array variables F_in and
      ! F_out must correspond to the same arrays used to generate the
      ! plans, but this is not checked.

      use iso_c_binding
      use fftw3
      implicit none

      type(dft_descriptor), intent(in) :: tdesc ! Transform descriptor
      real*8, intent(inout),target :: F_in(:,:,:), F_out(:,:,:) ! Input and output array (slices)

      real*8, dimension(:), pointer, contiguous :: in_ptr, out_ptr

      if (any(shape(F_in) /= tdesc%tsize)) then
         call gem_error(-1,'execute_r2r_dft_plan','Input array does not match plan size')
      elseif (any(shape(F_out) /= tdesc%tsize)) then
         call gem_error(-1,'execute_r2r_dft_plan','Output array does not match plan size')
      end if

      ! If tdesc%plan is null, then there is no transform to compute
      if (c_associated(tdesc%plan)) then
         ! Global plan, execute with internal parallelism
         call c_f_pointer(c_loc(F_in(1,1,1)),in_ptr,[1])
         call c_f_pointer(c_loc(F_out(1,1,1)),out_ptr,[1])
         call fftw_execute_r2r(tdesc%plan,in_ptr,out_ptr)
      end if
   end subroutine

   subroutine free_plans(tdesc)
      ! Free previously-allocated plans, returning any associated memory

      use iso_c_binding
      use fftw3
      implicit none
      type(dft_descriptor), intent(inout) :: tdesc

      if (c_associated(tdesc%plan)) then
!$omp critical(fftw_lock)
         call fftw_destroy_plan(tdesc%plan)
!$omp end critical(fftw_lock)
         tdesc%plan = C_NULL_PTR
      end if
   end subroutine


end module gem_fft
