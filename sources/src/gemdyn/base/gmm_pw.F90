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

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH) :: &
         gmmk_pw_uu_plus_s, gmmk_pw_uu_moins_s,&
         gmmk_pw_vv_plus_s, gmmk_pw_vv_moins_s,&
         gmmk_pw_wz_plus_s, gmmk_pw_wz_moins_s,&
         gmmk_pw_tt_plus_s, gmmk_pw_tt_moins_s,&
         gmmk_pw_pm_plus_s, gmmk_pw_pm_moins_s,&
         gmmk_pw_pt_plus_s, gmmk_pw_pt_moins_s,&
         gmmk_pw_gz_plus_s, gmmk_pw_gz_moins_s,&
         gmmk_pw_log_pm_s , gmmk_pw_log_pt_s  ,&
         gmmk_pw_uu_copy_s, gmmk_pw_vv_copy_s ,&
         gmmk_pw_me_plus_s, gmmk_pw_me_moins_s,&
         gmmk_pw_p0_plus_s, gmmk_pw_p0_moins_s,&
         gmmk_pw_p0_ls_s  , gmmk_pw_uslt_s    ,&
         gmmk_pw_vslt_s

end module gmm_pw
