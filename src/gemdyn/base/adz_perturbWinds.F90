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

subroutine adz_perturbWinds(F_alpha1, F_alpha2)

  use adz_mem, only: Adz_lminx, Adz_lmaxx, Adz_lminy, Adz_lmaxy,&
       Adz_halox, Adz_haloy, Adz_cy_8, Adz_uvw_d, Adz_uu_arr,   &
       Adz_vv_arr, Adz_ww_arr
  use glb_ld, only: l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, l_nk
  use dcst, only: Dcst_inv_rayt_8
  use ens_options, only: spp_L
  use ens_gmm_var, only: mcutraj, mcvtraj, mcwtraj
  use, intrinsic :: iso_fortran_env
  implicit none
#include <arch_specific.hf>

  real(kind=REAL64), intent(in) :: F_alpha1, F_alpha2

  integer :: i, j, k, nrow
  real, dimension(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1) :: &
       mcutraj_ext, mcvtraj_ext, mcwtraj_ext

  !
  !     ---------------------------------------------------------------
  !
  if (.not.spp_L) return

  mcutraj_ext = 0.
  mcvtraj_ext = 0.
  mcwtraj_ext = 0.
  nrow = 0
  call rpn_comm_xch_halox(mcutraj, l_minx, l_maxx, l_miny, l_maxy,  &
       l_ni, l_nj, 1, Adz_halox, Adz_haloy ,                        &
       .false., .false., mcutraj_ext, Adz_lminx ,                   &
       Adz_lmaxx, Adz_lminy, Adz_lmaxy, l_ni, nrow)
  call rpn_comm_xch_halox(mcvtraj, l_minx, l_maxx, l_miny, l_maxy,  &
       l_ni, l_nj, 1, Adz_halox, Adz_haloy ,                        &
       .false., .false., mcvtraj_ext, Adz_lminx ,                   &
       Adz_lmaxx, Adz_lminy, Adz_lmaxy, l_ni, nrow)
  call rpn_comm_xch_halox(mcwtraj, l_minx, l_maxx, l_miny, l_maxy,  &
       l_ni, l_nj, 1, Adz_halox, Adz_haloy ,                        &
       .false., .false., mcwtraj_ext, Adz_lminx ,                   &
       Adz_lmaxx, Adz_lminy, Adz_lmaxy, l_ni, nrow)

  do k = 1, l_nk
     do j = 1-Adz_haloy, l_nj+Adz_haloy
        do i = 1-Adz_halox, l_ni+Adz_halox
           Adz_uvw_d(1,i,j,k) = Adz_uvw_d(1,i,j,k) +                &
                Dcst_inv_rayt_8 * mcutraj_ext(i,j,1) * Adz_cy_8(j)
           Adz_uvw_d(2,i,j,k) = Adz_uvw_d(2,i,j,k) +                &
                Dcst_inv_rayt_8 * mcvtraj_ext(i,j,1)
           Adz_uvw_d(3,i,j,k) = Adz_uvw_d(3,i,j,k) * mcwtraj_ext(i,j,1)
        end do
     end do
  end do

  do k =1, l_nk
     do j =1, l_nj
        do i =1, l_ni
           Adz_uu_arr(i,j,k) = Adz_uu_arr(i,j,k) +                  &
                ( 2.*mcutraj_ext(i,j,1) * (F_alpha1 + F_alpha2) )   &
                * Dcst_inv_rayt_8 * Adz_cy_8(j)
           Adz_vv_arr(i,j,k) = Adz_vv_arr(i,j,k) +                  &
                ( 2.*mcvtraj_ext(i,j,1) * (F_alpha1 + F_alpha2) )   &
                * Dcst_inv_rayt_8
           Adz_ww_arr(i,j,k) = Adz_ww_arr(i,j,k) * mcwtraj_ext(i,j,1)
        end do
     end do
  end do
  !
  !     ---------------------------------------------------------------
  !
  return

end subroutine adz_perturbWinds
