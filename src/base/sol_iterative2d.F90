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
!**s/r sol_iterative_2d - Elliptic solver based on vertical decomposition leading to
!                         F_nk 2D horizontal elliptic problems to solve.
!
      subroutine sol_iterative2d ( F_rhs_sol, F_lhs_sol, F_ni, F_nj, F_nk, &
                                   print_conv_L, F_offi, F_offj )
      use dyn_fisl_options
      use lam_options
      use glb_ld
      use HORgrid_options
      use ldnh
      use lun
      use matvec2d_mod
      use opr
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: print_conv_L
      integer, intent(in) :: F_ni, F_nj, F_nk, F_offi, F_offj
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(in) :: F_rhs_sol
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(out) :: F_lhs_sol
!
!authors
!     Abdessamad Qaddouri, Stephane Gaudreault - March 2018
!     Christopher Subich - November 2018 (initial guess per level)
!
      integer :: i, j, k, ni, nij, nk
      real(kind=REAL64), dimension (ldnh_maxx, ldnh_maxy, l_nk) :: wk_rhs

      ! wk_sol is allocatable and saved between invocations, to provide an initial
      ! guess for subsequent calls
      real(kind=REAL64), dimension (:,:,:), allocatable, save :: wk_sol

      integer, dimension(F_nk) :: its
      real(kind=REAL64), dimension(F_nk) :: conv

      logical, save :: first_time = .true.

!
!     ---------------------------------------------------------------
!
      ni  = ldnh_ni - pil_w - pil_e
      nij = (ldnh_maxy - ldnh_miny + 1) * (ldnh_maxx - ldnh_minx + 1)

      wk_rhs = 0.d0
      if (.not. allocated(wk_sol)) then
         allocate(wk_sol(ldnh_maxx,ldnh_maxy,l_nk))
         wk_sol = 0.d0
      end if

      do j=1+pil_s, ldnh_nj-pil_n
         call dgemm ('N','N', ni, G_nk, G_nk, 1.d0, F_rhs_sol(1+pil_w,j,1), &
                     nij, Opr_lzevec_8, G_nk, 0.d0, wk_rhs(1+pil_w,j,1), nij)
         do k=1,Schm_nith
            do i = 1+pil_w, ldnh_ni-pil_e
               wk_rhs(i,j,k) = Opr_opsxp0_8(G_ni+F_offi+i) * &
                               Opr_opsyp0_8(G_nj+F_offj+j) * wk_rhs(i,j,k)
            end do
         end do
      end do

      nk = G_nk - Lam_gbpil_T

      if (first_time) then
         first_time = .false.
         call matvec2d_init()
      end if

      do k=1,F_nk
         ! On input, wk_sol(:,:,k) stores the initial guess for level k
         ! On output, wk_sol(:,:,k) stores the calculated solution for level k.

         ! This is also the initial guess for the next invocation of sol_iterative2d,
         ! which can come at the next nonlinear/trajectory iteration or at the next
         ! timestep.
         call sol_fgmres2d ( wk_sol(:,:,k), matvec2d_prod, wk_rhs(:,:,k), sol_fgm_eps, sol_im, sol_fgm_maxits, its(k), conv(k), k )
      end do

      if ( print_conv_L ) then
         do k=1,F_nk
            write(Lun_out,1004) k, its(k), conv(k)
         end do
      end if

      do j=1+pil_s, ldnh_nj-pil_n
         call dgemm ('N','T', ni, G_nk, G_nk, 1.d0, wk_sol(1+pil_w,j,1), &
                      nij, Opr_zevec_8, G_nk, 0.d0, F_lhs_sol(1+pil_w,j,1), nij)
      end do

 1004 format (3x,'level',i3, ' : FGMRES converged at iteration', i3,' to a solution with relative residual',1pe14.7)
!
!     ---------------------------------------------------------------
!
      return
      end


