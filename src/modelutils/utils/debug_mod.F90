!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module debug_mod
   use, intrinsic :: iso_fortran_env, only: REAL64
   use ieee_arithmetic
   implicit none

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   private
   public :: init2nan, assert_not_naninf

   interface init2nan
      module procedure init2nan_r4_s
      module procedure init2nan_r4_1d
      module procedure init2nan_r4_2d
      module procedure init2nan_r4_3d
      module procedure init2nan_r4_4d
      module procedure init2nan_r8_s
      module procedure init2nan_r8_1d
      module procedure init2nan_r8_2d
      module procedure init2nan_r8_3d
      module procedure init2nan_r8_4d
   end interface init2nan

   interface assert_not_naninf
      module procedure assert_not_naninf_r4
      module procedure assert_not_naninf_r4_1d
      module procedure assert_not_naninf_r4_2d
      module procedure assert_not_naninf_r4_3d
   end interface

   logical, save, public :: init2nan_L = .false.

   integer, parameter, public :: DEBUG_IS_NAN = RMN_ERR
   integer, parameter, public :: DEBUG_IS_INF = DEBUG_IS_NAN-1

!!$   !# From gmm_create
!!$    real      :: NaN
!!$    integer   :: inan
!!$    real(REAL64) ::    :: NaN8
!!$    integer(INT64) :: inan8
!!$
!!$    equivalence (inan,NaN)
!!$    equivalence (inan8,NaN8)
!!$
!!$    data inan  /Z'7F800001'/
!!$    data inan8 /Z'7FF0000000000001'/

#ifndef __GFORTRAN__
   real, parameter :: A4 = 1.
   real, parameter, public ::  DEBUG_NAN4 = transfer(-1, A4)
#else
#define DEBUG_NAN4 ieee_value(1.0, ieee_quiet_nan)
#endif
   
#ifndef __GFORTRAN__
   real(REAL64), parameter :: A8 = 1.D0
   real(REAL64), parameter, public ::  DEBUG_NAN8 = transfer(Z'7FF0000000000001', A8)
#else
#define DEBUG_NAN8 ieee_value(1.0_REAL64, ieee_quiet_nan)
#endif

contains

   subroutine init2nan_r4_s(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real, intent(inout) :: data
      real, intent(inout), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN4
      if (present(data2)) data2 = DEBUG_NAN4
      if (present(data3)) data3 = DEBUG_NAN4
      if (present(data4)) data4 = DEBUG_NAN4
      if (present(data5)) data5 = DEBUG_NAN4
      if (present(data6)) data6 = DEBUG_NAN4
      if (present(data7)) data7 = DEBUG_NAN4
      if (present(data8)) data8 = DEBUG_NAN4
      if (present(data9)) data9 = DEBUG_NAN4
      if (present(dataA)) dataA = DEBUG_NAN4
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r4_s


   subroutine init2nan_r4_1d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real, intent(inout), dimension(:) :: data
      real, intent(inout), dimension(:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN4
      if (present(data2)) data2 = DEBUG_NAN4
      if (present(data3)) data3 = DEBUG_NAN4
      if (present(data4)) data4 = DEBUG_NAN4
      if (present(data5)) data5 = DEBUG_NAN4
      if (present(data6)) data6 = DEBUG_NAN4
      if (present(data7)) data7 = DEBUG_NAN4
      if (present(data8)) data8 = DEBUG_NAN4
      if (present(data9)) data9 = DEBUG_NAN4
      if (present(dataA)) dataA = DEBUG_NAN4
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r4_1d


   subroutine init2nan_r4_2d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real, intent(inout), dimension(:,:) :: data
      real, intent(inout), dimension(:,:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN4
      if (present(data2)) data2 = DEBUG_NAN4
      if (present(data3)) data3 = DEBUG_NAN4
      if (present(data4)) data4 = DEBUG_NAN4
      if (present(data5)) data5 = DEBUG_NAN4
      if (present(data6)) data6 = DEBUG_NAN4
      if (present(data7)) data7 = DEBUG_NAN4
      if (present(data8)) data8 = DEBUG_NAN4
      if (present(data9)) data9 = DEBUG_NAN4
      if (present(dataA)) dataA = DEBUG_NAN4
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r4_2d


   subroutine init2nan_r4_3d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real, intent(inout), dimension(:,:,:) :: data
      real, intent(inout), dimension(:,:,:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN4
      if (present(data2)) data2 = DEBUG_NAN4
      if (present(data3)) data3 = DEBUG_NAN4
      if (present(data4)) data4 = DEBUG_NAN4
      if (present(data5)) data5 = DEBUG_NAN4
      if (present(data6)) data6 = DEBUG_NAN4
      if (present(data7)) data7 = DEBUG_NAN4
      if (present(data8)) data8 = DEBUG_NAN4
      if (present(data9)) data9 = DEBUG_NAN4
      if (present(dataA)) dataA = DEBUG_NAN4
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r4_3d


   subroutine init2nan_r4_4d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real, intent(inout), dimension(:,:,:,:) :: data
      real, intent(inout), dimension(:,:,:,:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN4
      if (present(data2)) data2 = DEBUG_NAN4
      if (present(data3)) data3 = DEBUG_NAN4
      if (present(data4)) data4 = DEBUG_NAN4
      if (present(data5)) data5 = DEBUG_NAN4
      if (present(data6)) data6 = DEBUG_NAN4
      if (present(data7)) data7 = DEBUG_NAN4
      if (present(data8)) data8 = DEBUG_NAN4
      if (present(data9)) data9 = DEBUG_NAN4
      if (present(dataA)) dataA = DEBUG_NAN4
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r4_4d


   subroutine init2nan_r8_s(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real(REAL64), intent(inout) :: data
      real(REAL64), intent(inout), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN8
      if (present(data2)) data2 = DEBUG_NAN8
      if (present(data3)) data3 = DEBUG_NAN8
      if (present(data4)) data4 = DEBUG_NAN8
      if (present(data5)) data5 = DEBUG_NAN8
      if (present(data6)) data6 = DEBUG_NAN8
      if (present(data7)) data7 = DEBUG_NAN8
      if (present(data8)) data8 = DEBUG_NAN8
      if (present(data9)) data9 = DEBUG_NAN8
      if (present(dataA)) dataA = DEBUG_NAN8
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r8_s


   subroutine init2nan_r8_1d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real(REAL64), intent(inout), dimension(:) :: data
      real(REAL64), intent(inout), dimension(:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN8
      if (present(data2)) data2 = DEBUG_NAN8
      if (present(data3)) data3 = DEBUG_NAN8
      if (present(data4)) data4 = DEBUG_NAN8
      if (present(data5)) data5 = DEBUG_NAN8
      if (present(data6)) data6 = DEBUG_NAN8
      if (present(data7)) data7 = DEBUG_NAN8
      if (present(data8)) data8 = DEBUG_NAN8
      if (present(data9)) data9 = DEBUG_NAN8
      if (present(dataA)) dataA = DEBUG_NAN8
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r8_1d


   subroutine init2nan_r8_2d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real(REAL64), intent(inout), dimension(:,:) :: data
      real(REAL64), intent(inout), dimension(:,:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN8
      if (present(data2)) data2 = DEBUG_NAN8
      if (present(data3)) data3 = DEBUG_NAN8
      if (present(data4)) data4 = DEBUG_NAN8
      if (present(data5)) data5 = DEBUG_NAN8
      if (present(data6)) data6 = DEBUG_NAN8
      if (present(data7)) data7 = DEBUG_NAN8
      if (present(data8)) data8 = DEBUG_NAN8
      if (present(data9)) data9 = DEBUG_NAN8
      if (present(dataA)) dataA = DEBUG_NAN8
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r8_2d


   subroutine init2nan_r8_3d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real(REAL64), intent(inout), dimension(:,:,:) :: data
      real(REAL64), intent(inout), dimension(:,:,:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN8
      if (present(data2)) data2 = DEBUG_NAN8
      if (present(data3)) data3 = DEBUG_NAN8
      if (present(data4)) data4 = DEBUG_NAN8
      if (present(data5)) data5 = DEBUG_NAN8
      if (present(data6)) data6 = DEBUG_NAN8
      if (present(data7)) data7 = DEBUG_NAN8
      if (present(data8)) data8 = DEBUG_NAN8
      if (present(data9)) data9 = DEBUG_NAN8
      if (present(dataA)) dataA = DEBUG_NAN8
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r8_3d


   subroutine init2nan_r8_4d(data, data2, data3, data4, data5, data6, data7, data8, data9, dataA)
      implicit none
      real(REAL64), intent(inout), dimension(:,:,:,:) :: data
      real(REAL64), intent(inout), dimension(:,:,:,:), optional :: data2, data3, data4, data5, data6, data7, data8, data9, dataA
      !----------------------------------------------------------------
      if (.not.init2nan_L) return
      data = DEBUG_NAN8
      if (present(data2)) data2 = DEBUG_NAN8
      if (present(data3)) data3 = DEBUG_NAN8
      if (present(data4)) data4 = DEBUG_NAN8
      if (present(data5)) data5 = DEBUG_NAN8
      if (present(data6)) data6 = DEBUG_NAN8
      if (present(data7)) data7 = DEBUG_NAN8
      if (present(data8)) data8 = DEBUG_NAN8
      if (present(data9)) data9 = DEBUG_NAN8
      if (present(dataA)) dataA = DEBUG_NAN8
      !----------------------------------------------------------------
      return
   end subroutine init2nan_r8_4d


   !/@
   function assert_not_naninf_r4(data) result(F_istat)
      implicit none
      real, intent(in) :: data
      integer :: F_istat
      !@/
!!$      real :: infinity_4
      ! ---------------------------------------------------------------------
      F_istat = RMN_OK
      !#uses fortran 2003 ieee_is_nan() and ieee_is_finite()
!!$      infinity_4 = huge(infinity_4)
!!$      if (data > infinity_4 .or. data < -infinity_4) then
      if (.not.ieee_is_finite(data)) F_istat = DEBUG_IS_INF
!!$      if (data /= data) then
      if (ieee_is_nan(data)) F_istat = DEBUG_IS_NAN
      ! ---------------------------------------------------------------------
      return
   end function assert_not_naninf_r4


   !/@
   function assert_not_naninf_r4_1d(data, F_i) result(F_istat)
      implicit none
      real, intent(in) :: data(:)
      integer, intent(out), optional :: F_i
      integer :: F_istat
      !@/
      integer :: i
      ! ---------------------------------------------------------------------
      F_istat = RMN_OK
      do i = lbound(data,1), ubound(data,1)
         if (.not.ieee_is_finite(data(i))) F_istat = DEBUG_IS_INF
         if (ieee_is_nan(data(i))) F_istat = DEBUG_IS_NAN
         if (.not.RMN_IS_OK(F_istat)) then
            if (present(F_i)) F_i = i
            return
         endif
      enddo
      ! ---------------------------------------------------------------------
      return
   end function assert_not_naninf_r4_1d


   !/@
   function assert_not_naninf_r4_2d(data, F_i, F_j) result(F_istat)
      implicit none
      real, intent(in) :: data(:,:)
      integer, intent(out), optional :: F_i, F_j
      integer :: F_istat
      !@/
      integer :: i, j
      ! ---------------------------------------------------------------------
      F_istat = RMN_OK
      do j = lbound(data,2), ubound(data,2)
         do i = lbound(data,1), ubound(data,1)
            if (.not.ieee_is_finite(data(i,j))) F_istat = DEBUG_IS_INF
            if (ieee_is_nan(data(i,j))) F_istat = DEBUG_IS_NAN
            if (.not.RMN_IS_OK(F_istat)) then
               if (present(F_i)) F_i = i
               if (present(F_j)) F_j = j
               return
            endif
         enddo
      enddo
      ! ---------------------------------------------------------------------
      return
   end function assert_not_naninf_r4_2d


   !/@
   function assert_not_naninf_r4_3d(data, F_i, F_j, F_k) result(F_istat)
      implicit none
      real, intent(in) :: data(:,:,:)
      integer, intent(out), optional :: F_i, F_j, F_k
      integer :: F_istat
      !@/
      integer :: i, j, k
      ! ---------------------------------------------------------------------
      F_istat = RMN_OK
      do k = lbound(data,3), ubound(data,3)
         do j = lbound(data,2), ubound(data,2)
            do i = lbound(data,1), ubound(data,1)
               if (.not.ieee_is_finite(data(i,j,k))) F_istat = DEBUG_IS_INF
               if (ieee_is_nan(data(i,j,k))) F_istat = DEBUG_IS_NAN
               if (.not.RMN_IS_OK(F_istat)) then
                  if (present(F_i)) F_i = i
                  if (present(F_j)) F_j = j
                  if (present(F_k)) F_k = k
                  return
               endif
            enddo
         enddo
      enddo
      ! ---------------------------------------------------------------------
      return
   end function assert_not_naninf_r4_3d


end module debug_mod
