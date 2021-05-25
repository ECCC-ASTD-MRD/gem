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
module gmm_pw
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      real, pointer, contiguous, dimension (:,:,:) :: pw_uu_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_vv_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_wz_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_tt_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_pm_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_pt_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_gz_plus  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_log_pm   => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_log_pt   => null()

      real, pointer, contiguous, dimension (:,:,:) :: pw_uu_moins => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_vv_moins => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_wz_moins => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_tt_moins => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_pm_moins => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_pt_moins => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_gz_moins => null()

      real, pointer, contiguous, dimension (:,:  ) :: pw_p0_ls    => null()
      real, pointer, contiguous, dimension (:,:  ) :: pw_me_plus  => null()
      real, pointer, contiguous, dimension (:,:  ) :: pw_p0_plus  => null()
      real, pointer, contiguous, dimension (:,:  ) :: pw_me_moins => null()
      real, pointer, contiguous, dimension (:,:  ) :: pw_p0_moins => null()

      real, pointer, contiguous, dimension (:,:  ) :: pw_uslt     => null()
      real, pointer, contiguous, dimension (:,:  ) :: pw_vslt     => null()

      real, pointer, contiguous, dimension (:,:,:) :: pw_uu_copy  => null()
      real, pointer, contiguous, dimension (:,:,:) :: pw_vv_copy  => null()

      real(kind=REAL64), pointer, dimension (:,:,:) :: pw_pm_moins_8 => null()
      real(kind=REAL64), pointer, dimension (:,:,:) :: pw_pm_plus_8  => null()
      real(kind=REAL64), pointer, dimension (:,:)   :: pw_p0_moins_8 => null()
      real(kind=REAL64), pointer, dimension (:,:)   :: pw_p0_plus_8  => null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_uu_plus_s = 'PW_UU:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_vv_plus_s = 'PW_VV:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_wz_plus_s = 'PW_WZ:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_tt_plus_s = 'PW_TT:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_pm_plus_s = 'PW_PM:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_pt_plus_s = 'PW_PT:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_gz_plus_s = 'PW_GZ:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_me_plus_s = 'PW_ME:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_p0_plus_s = 'PW_P0:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_log_pm_s = 'PW_LNPM:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_log_pt_s = 'PW_LNPT:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_pm_plus_8_s= 'PW_PM8:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_p0_plus_8_s= 'PW_P08:P'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_uu_moins_s = 'PW_UU:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_vv_moins_s = 'PW_VV:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_tt_moins_s = 'PW_TT:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_pm_moins_s = 'PW_PM:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_pt_moins_s = 'PW_PT:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_gz_moins_s = 'PW_GZ:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_me_moins_s = 'PW_ME:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_p0_moins_s = 'PW_P0:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_pm_moins_8_s= 'PW_PM8:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_p0_moins_8_s= 'PW_P08:M'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_uu_copy_s = 'PW_UU_COPY'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_vv_copy_s = 'PW_VV_COPY'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_uslt_s    = 'PW_USLT'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_vslt_s    = 'PW_VSLT'
      character(len=MAXNAMELENGTH), parameter:: gmmk_pw_p0_ls_s   = 'PW_P0_LS'

end module gmm_pw
