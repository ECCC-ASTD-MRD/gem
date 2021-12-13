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

!**s/r sol_init - initialize solver configuration
!
      integer function sol_init ()
      use dynkernel_options
      use dyn_fisl_options
      use fislh_sol
      use lun
      use ptopo
      implicit none

      character(len=16) :: dumc_S
!
!-------------------------------------------------------------------
!
      sol_init = -1
      call low2up  (Sol_type_S ,dumc_S)
      Sol_type_S = trim(dumc_S)
      call low2up  (Sol2D_precond_S ,dumc_S)
      Sol2D_precond_S = trim(dumc_S)
      call low2up  (Sol3D_precond_S ,dumc_S)
      Sol3D_precond_S = trim(dumc_S)

      FISLH_LHS_metric_L = .true.
      if (trim(Sol_type_S) == 'DIRECT') FISLH_LHS_metric_L = .false.

      isol_i = 1.0d0
      isol_d = 0.0d0

      if (.not. FISLH_LHS_metric_L) then
         isol_d = 1.0d0
         isol_i = 0.0d0
      endif

      select case(Sol_type_S)
      case('DIRECT')
         if (Sol_one_transpose_L .and. .not.Ptopo_alongY_L) then
            if (Lun_out > 0) &
            write(Lun_out, *) 'ONE-TRANSPOSE NOT ALLOWED UNLESS runmod.sh -along_Y'
            Sol_one_transpose_L = .false.
         endif

      case default
         select case(Sol_type_S)
            case('ITERATIVE_3D')
               if (Sol3D_krylov_S /= 'FGMRES' .and. Sol3D_krylov_S /= 'FBICGSTAB') then
               if (Lun_out > 0) &
               write(Lun_out, *) 'ABORT: WRONG CHOICE OF KRYLOV METHOD FOR 3D ITERATIVE SOLVER: Sol3D_krylov_S =', Sol3D_krylov_S
               return
               endif
            case('ITERATIVE_2D')
               if (sol2D_precond_S /= 'JACOBI'  .and.  sol2D_precond_S /= 'REDBLACK'  ) then
               if (Lun_out > 0) &
               write(Lun_out, *) 'ABORT: WRONG CHOICE OF PRECONDITIONER FOR 2D ITERATIVE SOLVER: sol2D_precond_S =', Sol2D_precond_S
               return
               endif
            case default
               if (lun_out>0) write (Lun_out, 9200) Sol_type_S
               return
         end select
      end select
      sol_init = 0

 9200 format (/'ABORT: WRONG CHOICE OF SOLVER for Helmholtz problem: Sol_type_S =',a/)
!
!-------------------------------------------------------------------
!
      return
      end function sol_init

