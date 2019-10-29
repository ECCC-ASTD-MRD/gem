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

!**s/r sol_main - Main driver for the elliptic solver
!
      subroutine sol_main ( F_rhs, F_solution, F_ni, F_nj, F_nk, F_conv_L )
      use ctrl
      use dyn_fisl_options
      use lun
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_conv_L
      integer, intent(in) :: F_ni, F_nj, F_nk
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(in) :: F_rhs
      real(kind=REAL64), dimension(F_ni,F_nj,F_nk), intent(out) :: F_solution
!author
!     Abdessamad Qaddouri -- January 2014
!
      integer :: offi, offj
!
!     ---------------------------------------------------------------
!
      if (ctrl_testcases_adv_L) return

      if (Lun_debug_L) write(Lun_out,1000)

      offi = Ptopo_gindx(1,Ptopo_myproc+1)-1
      offj = Ptopo_gindx(3,Ptopo_myproc+1)-1

      if (trim(Sol_type_S) == 'DIRECT') then
         call sol_direct ( F_rhs, F_solution, F_ni, F_nj, F_nk, &
                           F_conv_L, offi, offj )

      else if (trim(Sol_type_S) == 'ITERATIVE_2D') then
            call sol_iterative2d ( F_rhs, F_solution, F_ni, F_nj, F_nk, &
                                   F_conv_L, offi, offj )
      else ! 'ITERATIVE_3D'
            call sol_iterative3d ( F_rhs, F_solution, F_ni, F_nj, F_nk, &
                                   F_conv_L )
      end if

 1000 format( 5X,'SOLVING LINEAR HELMHOLTZ PROBLEM: (S/R SOL_MAIN)')
!
!     ---------------------------------------------------------------
!
      return
      end


