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

      use matvec, only: matvec_3d, matvec_yy
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

!author
!     Abdessamad Qaddouri -- January 2014
!
      integer :: its
      real(kind=REAL64) :: conv
      real(kind=REAL64), dimension(:,:,:), allocatable, save :: saved_sol
      integer :: i0, il, j0, jl

      integer fbicgstab3D ! TODO : use same signature as fgmres
      logical, save :: first_time = .true.

!
!     ---------------------------------------------------------------
!
!$omp single
!      call statf_dm (F_rhs_sol, 'RHS', 1, 'TSTP', 1,F_ni,1,F_nj,1,l_Nk,1,1,1,G_ni,G_nj,l_nk,8)
      i0 = 1    + pil_w
      il = l_ni - pil_e
      j0 = 1    + pil_s
      jl = l_nj - pil_n

      if (first_time) then
         allocate(saved_sol(F_ni,F_nj,F_nk))
         saved_sol = 0.d0
         first_time = .false.
      end if

      if (Grd_yinyang_L) then

         select case(Sol3D_krylov_S)
            case ('FGMRES')
               F_lhs_sol = saved_sol

                  call sol_fgmres3d (F_lhs_sol, matvec_yy, F_rhs_sol, sol_fgm_eps, sol_im, sol_fgm_maxits, its, conv)

               if ( F_print_L ) then
                  write(Lun_out, "(3x,'FGMRES converged at iteration', i4,' to a solution with relative residual',1pe14.7)") its, conv
               end if
               saved_sol = F_lhs_sol

            case ('FBICGSTAB')
               its = fbicgstab3d (F_lhs_sol, matvec_yy, F_rhs_sol, saved_sol,&
                                     l_ni,l_nj, F_nk, ldnh_minx, ldnh_maxx  ,&
                                     ldnh_miny, ldnh_maxy, i0, il, j0, jl   ,&
                                     sol_fgm_eps, sol_fgm_maxits            ,&
                                     Sol3D_precond_S, conv)

               if ( F_print_L ) then
                  write(Lun_out, "(3x,'FBiCGSTAB converged at iteration', i3,' to a solution with relative residual',1pe14.7)") its, conv
               end if
         end select
      else

         select case(Sol3D_krylov_S)
            case ('FGMRES')

            F_lhs_sol = saved_sol

            call sol_fgmres3d ( F_lhs_sol, matvec_3d, F_rhs_sol, sol_fgm_eps, sol_im, sol_fgm_maxits, its, conv)

            if ( F_print_L ) then
                 write(Lun_out, "(3x,'FGMRES converged at iteration', i3,' to a solution with relative residual',1pe14.7)") its, conv
            end if

            case ('FBICGSTAB')
               its = fbicgstab3d (F_lhs_sol, matvec_3d, F_rhs_sol, saved_sol,&
                                  l_ni,l_nj, F_nk, ldnh_minx, ldnh_maxx     ,&
                                  ldnh_miny, ldnh_maxy, i0, il, j0, jl      ,&
                                  sol_fgm_eps, sol_fgm_maxits               ,&
                                  Sol3D_precond_S, conv)

               if ( F_print_L ) then
                  write(Lun_out, "(3x,'FBiCGSTAB converged at iteration', i3,' to a solution with relative residual',1pe14.7)") its, conv
               end if
         end select
      end if

      saved_sol = F_lhs_sol
!      call statf_dm (F_lhs_sol, 'LHS', 1, 'TSTP', 1,F_ni,1,F_nj,1,l_Nk,1,1,1,G_ni,G_nj,l_nk,8)
!$omp end single

      return
      end subroutine sol_iterative3d
