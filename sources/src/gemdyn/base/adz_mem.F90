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

module adz_mem
   use ISO_C_BINDING
   use glb_ld
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      integer :: Adz_maxcfl           !Max CFL defined by user
      integer :: Adz_halox, Adz_haloy !Adzection x/y halo size
      integer :: Adz_lminx, Adz_lmaxx , Adz_lminy, Adz_lmaxy
      integer :: Adz_nit, Adz_njt, Adz_nij, Adz_ioff, Adz_joff
      integer :: Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_i0b,Adz_inb,Adz_j0b,Adz_jnb
      integer :: Adz_i0u,Adz_inu,Adz_j0v,Adz_jnv
      integer :: Adz_2dnh, Adz_3dnh, Adz_k0, Adz_k0t, Adz_k0m
      integer :: Adz_num_u,Adz_num_v,Adz_num_q,Adz_num_t,Adz_num_b,Adz_kkmax
      integer, dimension(:), pointer, contiguous :: Adz_search_m,&
                                                    Adz_search_t
      real :: Adz_iminposx,Adz_imaxposx,Adz_iminposy,Adz_imaxposy
      real(kind=REAL64) :: adz_ovdzm_8, adz_ovdzt_8
      real(kind=REAL64), dimension(:), pointer :: Adz_cy_8, &
            Adz_delz_m, Adz_delz_t, Adz_odelz_m, Adz_odelz_t

      type(C_PTR) :: Adz_cpntr_q,Adz_cpntr_t

      type Adz_pntr_stack
         sequence
         real, dimension(:,:,:), pointer :: src,dst,pil
      end type Adz_pntr_stack

      type(Adz_pntr_stack), pointer, dimension(:) :: Adz_stack

! Adz_pxyzm is a history carrying variable: must be part of the restart

      real, pointer    , dimension (:,:,:,:), contiguous :: &
                              Adz_pxyzm=>null()
      real, allocatable, dimension (:,:,:,:)             :: &
                              Adz_pxyzt,Adz_uvw_dep
      real, allocatable, dimension (:,:,:,:)             :: &
                           Adz_pm,Adz_pmu,Adz_pmv,Adz_pt,Adz_pb

      real(kind=REAL64), allocatable, dimension (:,:,:,:) :: Adz_wpxyz

      real, allocatable, dimension (:,:,:  ) ::       &
                     Adz_uu_ext,Adz_vv_ext,Adz_ww_ext,&
                     Adz_uu_arr,Adz_vv_arr,Adz_ww_arr

      real, allocatable, dimension (:,:,:,:) ::  Adz_uvw_d

      real, pointer, dimension (:,:), contiguous :: Adz_uslt  =>null()
      real, pointer, dimension (:,:), contiguous :: Adz_vslt  =>null()

      real, pointer, dimension (:,:,:,:), contiguous :: Adz_post,Adz_flux

end module adz_mem
