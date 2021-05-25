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

      logical :: ADZ_OD_L
      integer :: Adz_halox, Adz_haloy, Adz_MAX_MPI_OS_SIZE

! WARNING These next variables will be defined and used
! only in the context of ADZ_OD_L=.false.
      integer :: Adz_maxcfl
      integer :: Adz_lminx, Adz_lmaxx, Adz_lminy, Adz_lmaxy
      integer :: Adz_nit, Adz_njt, Adz_nij
      real, allocatable, dimension (:,:,:  ) ::       &
                     Adz_uu_ext,Adz_vv_ext,Adz_ww_ext
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer adz_COMM, INTEGER_DATATYPE, REAL_DATATYPE
      integer :: Adz_ioff, Adz_joff
      integer :: Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_i0b,Adz_inb,Adz_j0b,Adz_jnb
      integer :: Adz_i0u,Adz_inu,Adz_j0v,Adz_jnv
      integer :: Adz_2dnh, Adz_3dnh, Adz_k0, Adz_k0t, Adz_k0m
      integer :: Adz_num_u,Adz_num_v,Adz_num_q,Adz_num_t,Adz_num_b,Adz_kkmax
      integer :: Adz_icn, Adz_itnl, Adz_cnt_traj, Adz_cnt_int
      integer, dimension(:), allocatable :: Adz_gindx_alongX, Adz_gindx_alongY
      integer, dimension(:), pointer, contiguous :: Adz_search_m,&
                                                    Adz_search_t
      real, dimension(:  ), allocatable :: nexports
      real(kind=REAL64) :: Adz_iminposx,Adz_imaxposx,Adz_iminposy,Adz_imaxposy
      real(kind=REAL64) :: Adz_yyminposx,Adz_yymaxposx,Adz_yyminposy,Adz_yymaxposy
      real(kind=REAL64) :: adz_ovdzm_8, adz_ovdzt_8
      real(kind=REAL64), dimension(:), pointer :: Adz_cy_8, &
              Adz_delz_m, Adz_delz_t, Adz_odelz_m, Adz_odelz_t
      real(kind=REAL64) :: Adz_glbminpos

      type :: ADZ_TYPE_SV_R8
         sequence
         real(kind=REAL64), dimension(:), pointer :: t,m
      end type ADZ_TYPE_SV_R8

      type(ADZ_TYPE_SV_R8) :: Adz_zabcd_8, Adz_zbacd_8, &
                              Adz_zcabd_8, Adz_zdabc_8

      type(ADZ_TYPE_SV_R8) :: Adz_zxabcde_8, Adz_zaxbcde_8, Adz_zbxacde_8, &
                              Adz_zcxabde_8, Adz_zdxabce_8, Adz_zexabcd_8

      type(C_PTR) :: Adz_cpntr_q,Adz_cpntr_t

      type Adz_pntr_stack
         real, dimension(:,:,:), pointer :: src,dst,pil
      end type Adz_pntr_stack

      type(Adz_pntr_stack), pointer, dimension(:) :: Adz_stack

! Adz_pxyzm and Adz_wpxyz are a history carrying variables: must be part of the restart

      real, pointer, dimension (:,:,:,:), contiguous  :: Adz_pxyzm=>null()
      real(kind=REAL64), pointer, dimension (:,:,:,:) :: Adz_wpxyz=>null()

      real, allocatable, dimension (:,:,:,:) :: Adz_uvw_dep,&
                         Adz_pm,Adz_pmu,Adz_pmv,Adz_pt,Adz_pb

      real, allocatable, dimension (:,:,:  ) ::       &
                     Adz_uu_arr,Adz_vv_arr,Adz_ww_arr

      real, allocatable, dimension (:,:,:,:) ::  Adz_uvw_d

      real, pointer, dimension (:,:), contiguous :: Adz_uslt  =>null()
      real, pointer, dimension (:,:), contiguous :: Adz_vslt  =>null()

      type post
           real, dimension(:,:,:), pointer :: lin
           real, dimension(:,:,:), pointer :: min
           real, dimension(:,:,:), pointer :: max
      end type post

      type(post), pointer, dimension(:) :: Adz_post,Adz_post_3CWP,Adz_post_BQWP

      type flux
           real, dimension(:,:,:), pointer :: fo
           real, dimension(:,:,:), pointer :: fi
      end type flux

      type(flux), pointer, dimension(:) :: Adz_flux,Adz_flux_3CWP,Adz_flux_BQWP,Adz_flux_3CWP_PS

      type bc
           integer :: n
           integer :: wp
           integer :: pexp_n
           integer :: weight
           character(len=8)  :: name_S
           real(kind=REAL64) :: mass_deficit
           logical :: LEGACY_L
           real, dimension(:,:,:), pointer :: p
           real, dimension(:,:,:), pointer :: m
           real, dimension(:,:,:), pointer :: w
           real, dimension(:,:,:), pointer :: fo
           real, dimension(:,:,:), pointer :: fi
           real, dimension(:,:,:), pointer :: lin
           real, dimension(:,:,:), pointer :: min
           real, dimension(:,:,:), pointer :: max
      end type bc

      type(bc), pointer, dimension(:) :: Adz_bc

      type :: ADZ_SLOD
         integer npe,winreqs,wintraj
         integer, dimension (:,:), allocatable :: stk,dest
         integer, dimension (:  ), pointer     :: from,list
         real   , dimension (:  ), pointer     :: requests,pos
         real   , dimension (:,:), allocatable :: gpos
      end type ADZ_SLOD
      type(ADZ_SLOD) :: Adz_expq, Adz_expu, Adz_expv, Adz_expt

      integer :: Adz_Win_list, Adz_Win_pos
      integer, dimension (:,:), pointer :: Adz_od_list
      real   , dimension (:,:), pointer :: Adz_od_pos

      integer :: Adz_offs_win, Adz_wincor
      integer, dimension (:), pointer :: Adz_offs

      real   , dimension (:), pointer :: Adz_cor
      
      integer :: ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK, &
                 TARGET_COUNT, TARGET_DATATYPE
      integer(C_INTPTR_T) :: TARGET_DISP
      
      contains
      
!!$integer function RPN_COMM_locktype(locktype)
!!$  implicit none
!!$  include 'mpif.h'
!!$  character(len=*), intent(IN) :: locktype
!!$  if(locktype == 'MPI_LOCK_EXCLUSIVE') then
!!$    RPN_COMM_locktype = MPI_LOCK_EXCLUSIVE
!!$  else
!!$    RPN_COMM_locktype = MPI_LOCK_SHARED
!!$  endif
!!$  return
!!$end function RPN_COMM_locktype
      
end module adz_mem
