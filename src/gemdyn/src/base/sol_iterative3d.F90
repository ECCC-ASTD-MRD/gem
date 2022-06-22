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
!**s/r sol_iterative3d - Full 3D iterative elliptic solver
!                        Available for (LAM /Yin-Yang) configurations

   subroutine sol_iterative3d ( F_rhs_sol, F_lhs_sol, F_ni, F_nj, F_nk, &
                                F_print_L)

      use Horgrid_options
      use gem_options
      use dyn_fisl_options
      use glb_ld
      use lun
      use ldnh
      use sol
      use opr
      use dynkernel_options
      use stat_mpi
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_print_L
      integer, intent(in) :: F_ni, F_nj, F_nk
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(in) ::  F_rhs_sol
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(inout) ::  F_lhs_sol

      integer :: i,j,k,its
      real(kind=REAL64) :: conv

!
!     ---------------------------------------------------------------
!
!$omp single
!      call statf_dm (F_rhs_sol, 'RHS', 1, 'TSTP', 1,F_ni,1,F_nj,1,l_Nk,1,1,1,G_ni,G_nj,l_nk,8)
!$omp end single

!$omp single
      do k=1, F_nk
         do j=1,F_nj
            do i=1,F_ni
               F_lhs_sol(i,j,k) = Sol_saved(i,j,k)
            enddo
         enddo
      enddo
!$omp end single

      call sol_fgmres3d (F_lhs_sol, F_rhs_sol, sol_fgm_eps, sol_im, sol_fgm_maxits, its, conv)

!$omp single
      if ( F_print_L ) then
         write(Lun_out, "(3x,'FGMRES converged at iteration', i4,' to a solution with relative residual',1pe14.7)") its, conv
      end if
!$omp end single

!$omp single
      do k=1, F_nk
         do j=1,F_nj
            do i=1,F_ni
               Sol_saved(i,j,k) = F_lhs_sol(i,j,k)
            enddo
         enddo
      enddo
!$omp end single


!$omp single
!      call statf_dm (F_lhs_sol, 'LHS', 1, 'TSTP', 1,F_ni,1,F_nj,1,l_Nk,1,1,1,G_ni,G_nj,l_nk,8)
!$omp end single

      return
      end subroutine sol_iterative3d
