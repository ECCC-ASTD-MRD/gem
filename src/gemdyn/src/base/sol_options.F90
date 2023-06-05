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
module sol_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Type of solver
   !# * 'ITERATIF'
   !# * 'DIRECT'
   character(len=26) :: Sol_type_S = 'ITERATIVE_3D'
   namelist /sol  / Sol_type_S

   !# Relative tolerance for convergence of the elliptic iterative solver, 
   !# norm(residual) < tol*norm(b)
   real(kind=REAL64) :: Sol_fgm_eps   = 1.d-07
   namelist /sol  / Sol_fgm_eps

   !# Absolute tolerance for convergence of the classical Schwarz algorithm 
   !# in Yin Yang configurations, norm(residual) < tol
   real(kind=REAL64) :: Sol_yyg_eps   = 1.d-04
   namelist /sol  / Sol_yyg_eps
   namelist /sol_p/ Sol_yyg_eps

   integer :: sol_fgm_maxits= 200
   namelist /sol  / Sol_fgm_maxits

   !# maximum number of iterations allowed for the Yin-Yang iterative solver
   integer :: sol_yyg_maxits= 40
   namelist /sol  / Sol_yyg_maxits

   !# size of Krylov subspace in iterative solver - should not exceed 100
   integer :: sol_im = 15
   namelist /sol  / Sol_im

   !# 3D preconditioner for iterative solver
   character(len=26) :: Sol_precond3D_S = 'RAS'
   namelist /sol  / Sol_precond3D_S

   !# 3D preconditioner for iterative solver
   integer :: Sol_gauss_Niter = 2
   namelist /sol  / Sol_gauss_Niter

   !# Krylov method for 3d iterative solver (FGMRES or FBICGSTAB)
   character(len=26) :: Sol_krylov3D_S = 'FGMRES'
   namelist /sol  / Sol_krylov3D_S

   !# 3D preconditioner for iterative solver
    integer :: Sol_ovlpx = 4
    namelist /sol  / Sol_ovlpx

   !# 3D preconditioner for iterative solver
    integer :: Sol_ovlpy = 4
    namelist /sol  / Sol_ovlpy

   !# True => use the one-transpose solver
   logical :: Sol_one_transpose_L = .false.
   namelist /sol  / Sol_one_transpose_L

contains

!**s/r sol_nml - Read namelist sol

      integer function sol_nml (F_unf)
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <clib_interface_mu.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer err
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         sol_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=sol)
         end if
         return
      end if

      sol_nml= -1 ; nml_must= .false. ; nml_S= 'sol'

      rewind(F_unf)
      read (F_unf, nml=sol, end= 1001, err=1003)
      sol_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         sol_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (sol_nml < 0 ) return
      if ((Lun_out>=0).and.(sol_nml==0)) write (Lun_out, 6004) trim(nml_S)
      sol_nml= 1

      err = clib_toupper(Sol_type_S)
      err = clib_toupper(Sol_precond3D_S)
      err = clib_toupper(Sol_krylov3D_S)

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
 6009 format (//,' NAMELIST ',A,' IS MANDATORY'//)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function sol_nml

end module sol_options
