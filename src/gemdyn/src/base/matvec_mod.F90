module matvec
   ! Matrix-vector product subroutines
   !
   ! Author: Abdessamad Qaddouri
   !
   use cstv
   use geomh
   use glb_ld
   use ldnh
   use opr
   use sol_mem
   use, intrinsic :: iso_fortran_env
   implicit none
#include <arch_specific.hf>
   integer, parameter :: IDX_POINT=1, IDX_WEST=2, IDX_EAST=3, IDX_NORTH=4, IDX_SOUTH=5, IDX_TOP=6, IDX_BOTTOM=7
   real(kind=REAL64), dimension(:,:,:,:), allocatable :: stencil

   public :: matvec_init, matvec_3d, matvec_yy
   ! Matrix and the other indices are made public for use in the preconditioner
   public :: IDX_POINT, IDX_WEST, IDX_EAST, IDX_NORTH, IDX_SOUTH, IDX_TOP, IDX_BOTTOM, stencil

contains

   subroutine matvec_init()
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64)  :: di_8
      real(kind=REAL64)  :: xxx, yyy
      integer :: i, j, k, jj, ii

      if (allocated(stencil)) deallocate (stencil)
      allocate (stencil(1+sol_pil_w:l_ni-sol_pil_e, 1+sol_pil_s:l_nj-sol_pil_n, 7, l_nk))

      xxx = - Cstv_hco2_8
      yyy = - Cstv_hco1_8

      do k=1, l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            jj=j+l_j0-1
            di_8 = Opr_opsyp0_8(G_nj+jj) * geomh_invcy2_8(j)
            do i=1+sol_pil_w, l_ni-sol_pil_e
               ii=i+l_i0-1

               stencil(i,j,IDX_POINT,k) = Cstv_hco0_8 * (Opr_opszp2_8(G_nk+k) + Opr_opszpl_8(G_nk+k) &
                           + xxx * Opr_opszpm_8(G_nk+k) + yyy * Opr_opszp0_8(G_nk+k)) &
                           + Opr_opszp0_8(G_nk+k) * (Opr_opsxp2_8(G_ni+ii) * di_8     &
                           + Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(G_nj+jj))           &
                           / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               stencil(i,j,IDX_WEST,k) = Opr_opsxp2_8(ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                           / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii)*Opr_opsyp0_8(G_nj+jj))

               stencil(i,j,IDX_EAST,k) = Opr_opsxp2_8(2*G_ni+ii) * Opr_opszp0_8(G_nk+k) * Opr_opsyp0_8(G_nj+jj) &
                          / (cos( G_yg_8 (jj) )**2) / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               stencil(i,j,IDX_SOUTH,k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(jj) * Opr_opszp0_8(G_nk+k) &
                          / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               stencil(i,j,IDX_NORTH,k) = Opr_opsxp0_8(G_ni+ii) * Opr_opsyp2_8(2*G_nj+jj) * Opr_opszp0_8(G_nk+k) &
                           / (Opr_opsxp0_8(G_ni+ii) * Opr_opsyp0_8(G_nj+jj))

               stencil(i,j,IDX_TOP,k) = Cstv_hco0_8 * (Opr_opszp2_8(k) + Opr_opszpl_8(k) + xxx * Opr_opszpm_8(k))

               stencil(i,j,IDX_BOTTOM,k) = Cstv_hco0_8 * (Opr_opszp2_8(2*G_nk+k) + Opr_opszpl_8(2*G_nk+k) + xxx * Opr_opszpm_8(2*G_nk+k))

            end do
         end do
      end do

   end subroutine matvec_init


   subroutine matvec_3d(vec, prod)
      use, intrinsic :: iso_fortran_env
      implicit none
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(in) :: vec
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy, l_nk), intent(out) :: prod

      integer :: i, j, k
      real(kind=REAL64) :: vector(0:l_ni+1, 0:l_nj+1, l_nk)

      do k = 1, l_nk
         vector(:,0:sol_pil_s,k) = 0
         do j=1+sol_pil_s, l_nj-sol_pil_n
            vector(0:sol_pil_w,j,k) = 0
            do i=1+sol_pil_w, l_ni-sol_pil_e
               vector(i,j,k) = vec(i,j,k)
            end do
            vector((l_ni-sol_pil_e+1):(l_ni+1),j,k) = 0
         end do
         vector(:,(l_nj-sol_pil_n+1):(l_nj+1),k) = 0
      end do

      call rpn_comm_xch_halon (vector, 0, l_ni+1, 0, l_nj+1, l_ni, l_nj, l_nk, &
                               1, 1, G_periodx, G_periody, l_ni, 0, 2)

      do k=1,l_nk
         if (k == 1) then
            do j=1+sol_pil_s, l_nj-sol_pil_n
               do i=1+sol_pil_w, l_ni-sol_pil_e

                  prod(i,j,k) = stencil(i,j,IDX_POINT,k)  * vector(i  ,j  ,k  ) + &
                                stencil(i,j,IDX_WEST,k)   * vector(i-1,j  ,k  ) + &
                                stencil(i,j,IDX_EAST,k)   * vector(i+1,j  ,k  ) + &
                                stencil(i,j,IDX_NORTH,k)  * vector(i  ,j+1,k  ) + &
                                stencil(i,j,IDX_SOUTH,k)  * vector(i  ,j-1,k  ) + &
                                stencil(i,j,IDX_BOTTOM,k) * vector(i  ,j  ,k+1)
               end do
            end do
         elseif (k > 1 .and. k < l_nk) then
            do j=1+sol_pil_s, l_nj-sol_pil_n
               do i=1+sol_pil_w, l_ni-sol_pil_e

                  prod(i,j,k) = stencil(i,j,IDX_POINT,k)  * vector(i  ,j  ,k  ) + &
                                stencil(i,j,IDX_WEST,k)   * vector(i-1,j  ,k  ) + &
                                stencil(i,j,IDX_EAST,k)   * vector(i+1,j  ,k  ) + &
                                stencil(i,j,IDX_NORTH,k)  * vector(i  ,j+1,k  ) + &
                                stencil(i,j,IDX_SOUTH,k)  * vector(i  ,j-1,k  ) + &
                                stencil(i,j,IDX_TOP,k)    * vector(i  ,j  ,k-1) + &
                                stencil(i,j,IDX_BOTTOM,k) * vector(i  ,j  ,k+1)
               end do
            end do
         else ! k == l_nk
            do j=1+sol_pil_s, l_nj-sol_pil_n
               do i=1+sol_pil_w, l_ni-sol_pil_e

                  prod(i,j,k) = stencil(i,j,IDX_POINT,k)  * vector(i  ,j  ,k  ) + &
                                stencil(i,j,IDX_WEST,k)   * vector(i-1,j  ,k  ) + &
                                stencil(i,j,IDX_EAST,k)   * vector(i+1,j  ,k  ) + &
                                stencil(i,j,IDX_NORTH,k)  * vector(i  ,j+1,k  ) + &
                                stencil(i,j,IDX_SOUTH,k)  * vector(i  ,j-1,k  ) + &
                                stencil(i,j,IDX_TOP,k)    * vector(i  ,j  ,k-1)
               end do
            end do
         end if
      end do


   end subroutine matvec_3d

   ! matvec_yy (was matvec3d_prod) -- implement matrix-vector product with yin/yang grid exchange
   subroutine matvec_yy(F_vector, F_prod)
      use glb_ld      , only: l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj,l_nk,G_periodx, G_periody,l_i0,l_j0
      use ldnh        , only: ldnh_minx, ldnh_maxx, ldnh_miny, ldnh_maxy
      use sol_mem
      use gem_options , only: G_halox, G_haloy
      use geomh
      use HORgrid_options, only: Grd_yinyang_L

      use, intrinsic :: iso_fortran_env
      implicit none
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy,l_nk), intent(in) :: F_vector
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy,l_nk), intent(out) :: F_prod

      integer :: i, j, k, ii, jj
      real linfini
      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy,0:l_nk+1) :: vector
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy,l_nk) :: work_8

      vector = 0.0d0
      do k=1,l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               vector( i, j,k) = F_vector(i, j,k)
            end do
         end do
      end do

      call rpn_comm_xch_halo_8 (vector(:,:,0:l_nk+1),l_minx, l_maxx, l_miny, l_maxy, l_ni, l_nj, l_nk+2, &
                           G_halox, G_haloy, G_periodx, G_periody, l_ni,0 )


      do k=1,l_nk
         do j=1+sol_pil_s, l_nj-sol_pil_n
            jj = j + l_j0 - 1
            do i=1+sol_pil_w, l_ni-sol_pil_e
               ii = i + l_i0 - 1

               F_prod(i,j,k) = stencil(i, j, IDX_POINT, k) * vector(i  , j ,k ) &
                             + stencil(i, j, IDX_WEST,  k) * vector(i-1, j ,k ) &
                             + stencil(i, j, IDX_EAST,  k) * vector(i+1, j ,k ) &
                             + stencil(i, j, IDX_NORTH, k) * vector(i  , j+1,k) &
                             + stencil(i, j, IDX_SOUTH, k) * vector(i  , j-1,k) &
                             + stencil(i, j, IDX_TOP, k)   * vector(i  ,j  ,k-1)   &
                             + stencil(i, j, IDX_BOTTOM, k)  *vector(i  ,j  ,k+1)
            end do
         end do
     end do

      if (Grd_yinyang_L) then
         work_8 =0.0
         call yyg_rhs_xchng (work_8, vector(l_minx,l_miny,1)        ,&
                             ldnh_minx,ldnh_maxx,ldnh_miny,ldnh_maxy,&
                             l_minx,l_maxx,l_miny,l_maxy,l_nk,-1,linfini)
         do k=1,l_nk
            do j=1+sol_pil_s, l_nj-sol_pil_n
               jj=j+l_j0-1
               do i=1+sol_pil_w, l_ni-sol_pil_e
                  ii = i + l_i0 - 1
                  F_prod(i,j,k)= F_prod(i,j,k)+work_8(i,j,k)/ geomh_area_8(i,j)

               end do
            end do
         end do
      end if

   end subroutine matvec_yy
end module matvec
