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
!----

module matvec2d_mod
   ! stencil-vector product subroutines for 2D elliptic problem
   !
   ! Author
   !     Abdessamad Qaddouri, Stéphane Gaudreault - March 2018
   !

   implicit none
   private

#include <arch_specific.hf>

   integer, parameter :: IDX_POINT=1, IDX_WEST=2, IDX_EAST=3, IDX_NORTH=4, IDX_SOUTH=5
   real*8, dimension(:,:,:,:), allocatable :: stencil

   public :: matvec2d_init, matvec2d_prod

contains

   subroutine matvec2d_init()
      use cstv
      use geomh
      use glb_ld
      use opr
      use sol
      implicit none

      real*8  :: cst, di_8
      integer :: i, j, k, jj, ii

      allocate (stencil(1+sol_pil_w:l_ni-sol_pil_e, 1+sol_pil_s:l_nj-sol_pil_n, 5, l_nk))


!$omp parallel private (i,j,k,jj,ii,cst,di_8)
!$omp do
      do  k=1, l_nk

         cst = Cstv_hco1_8 + Cstv_hco0_8 * Opr_zeval_8(k)

         do j = 1+sol_pil_s, l_nj-sol_pil_n
            jj = j + l_j0 - 1
            di_8 = Opr_opsyp0_8(G_nj+jj) * geomh_invcy2_8(j)

            do i = 1+sol_pil_w, l_ni-sol_pil_e
               ii = i + l_i0 - 1

               stencil(i, j, IDX_POINT, k) = cst * Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj) &
                                           + Opr_opsxp2_8(G_ni+ii) * di_8                        &
                                           + Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(G_nj+jj)

               stencil(i, j, IDX_WEST, k)  = Opr_opsxp2_8(ii) * di_8

               stencil(i, j, IDX_EAST, k)  = Opr_opsxp2_8(2*G_ni+ii) * di_8

               stencil(i, j, IDX_SOUTH, k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(jj)

               stencil(i, j, IDX_NORTH, k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(2*G_nj+jj)

            end do
         end do
      end do
!$omp enddo
!$omp end parallel
   end subroutine matvec2d_init


   subroutine matvec2d_prod(F_vector, F_prod, level)
      use glb_ld         , only: l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, G_periodx, G_periody
      use ldnh           , only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy
      use sol            , only: sol_pil_s, sol_pil_n, sol_pil_w, sol_pil_e
      use gem_options    , only: G_halox, G_haloy
      use HORgrid_options, only: Grd_yinyang_L
      implicit none
      integer, intent(in) :: level
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(in) :: F_vector
      real*8, dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy), intent(out) :: F_prod

      integer :: i, j, halox, haloy
      real*8, dimension(l_minx:l_maxx,l_miny:l_maxy) :: vector
      real linfini

      F_prod = 0.d0

      do j=1,l_nj
!DIR$ SIMD
         do i=1,l_ni
            vector(i, j) = F_vector(i, j)
         end do
      end do

      if (Grd_yinyang_L) then
         halox = G_halox
         haloy = G_haloy
      else
         halox = 1
         haloy = 1
      end if

      call rpn_comm_xch_halo_8 (vector, l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, 1, &
                                G_halox, G_haloy, G_periodx, G_periody, l_ni,0 )

      do j=1+sol_pil_s, l_nj-sol_pil_n
!DIR$ SIMD
         do i=1+sol_pil_w, l_ni-sol_pil_e
!DIR$ VECTOR NONTEMPORAL(F_prod)
            F_prod(i,j) = stencil(i, j, IDX_POINT, level) * vector(i  , j  ) &
                        + stencil(i, j, IDX_WEST,  level) * vector(i-1, j  ) &
                        + stencil(i, j, IDX_EAST,  level) * vector(i+1, j  ) &
                        + stencil(i, j, IDX_NORTH, level) * vector(i  , j+1) &
                        + stencil(i, j, IDX_SOUTH, level) * vector(i  , j-1)
         end do
      end do

      if (Grd_yinyang_L) then
         call yyg_rhs_xchng (F_prod, vector                         ,&
                             ldnh_minx,ldnh_maxx,ldnh_miny,ldnh_maxy,&
                             l_minx,l_maxx,l_miny,l_maxy,1,-1,linfini)
      end if

   end subroutine matvec2d_prod

end module matvec2d_mod
