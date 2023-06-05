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

!** matvec Matrix-vector product subroutines without vertical metrics
!
      subroutine matvec_p ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                              F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use dyn_fisl_options
      use geomh
      use gem_options
      use HORgrid_options
      use dynkernel_options
      use lam_options
      use glb_ld
      use ldnh
      use sol_mem
      use metric
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn,F_nk), intent(out) :: F_prod

      integer, parameter :: IDX_POINT=1, IDX_WEST=2, IDX_EAST=3, IDX_NORTH=4, IDX_SOUTH=5, IDX_TOP=6, IDX_BOTTOM=7
      integer :: i, j, k
      real linfini
      real(kind=REAL64), dimension(l_minx:l_maxx, l_miny:l_maxy,0:l_nk+1) :: vector
      real(kind=REAL64), dimension(ldnh_minx:ldnh_maxx, ldnh_miny:ldnh_maxy,l_nk) :: work_8
!
!     ---------------------------------------------------------------
!
!$omp single
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
               do i=1+sol_pil_w, l_ni-sol_pil_e
                  F_prod(i,j,k) = Sol_stencilp_8(i, j, IDX_POINT, k) * vector(i  , j ,k ) &
                                + Sol_stencilp_8(i, j, IDX_WEST,  k) * vector(i-1, j ,k ) &
                                + Sol_stencilp_8(i, j, IDX_EAST,  k) * vector(i+1, j ,k ) &
                                + Sol_stencilp_8(i, j, IDX_NORTH, k) * vector(i  , j+1,k) &
                                + Sol_stencilp_8(i, j, IDX_SOUTH, k) * vector(i  , j-1,k) &
                                + Sol_stencilp_8(i, j, IDX_TOP, k)   * vector(i  ,j  ,k-1)&
                                + Sol_stencilp_8(i, j, IDX_BOTTOM, k)* vector(i  ,j  ,k+1)
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
                  do i=1+sol_pil_w, l_ni-sol_pil_e
                     F_prod(i,j,k)= F_prod(i,j,k)+work_8(i,j,k)/ geomh_area_8(i,j)
                  end do
               end do
            end do
         end if
!$omp end single
!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine matvec_p


