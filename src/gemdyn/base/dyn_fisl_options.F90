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
module dyn_fisl_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# SL off-centering parameter for hydrostatic
   real(kind=REAL64) :: Cstv_bA_8 = 0.6
   namelist /dyn_fisl  / Cstv_bA_8
   namelist /dyn_fisl_p/ Cstv_bA_8

   !# SL off-centering parameter for the momentum equations
   real(kind=REAL64) :: Cstv_bA_m_8 = 0.6
   namelist /dyn_fisl  / Cstv_bA_m_8
   namelist /dyn_fisl_p/ Cstv_bA_m_8

   !# SL off-centering parameter for nonhydrostatic
   real(kind=REAL64) :: Cstv_bA_nh_8 = 0.5
   namelist /dyn_fisl  / Cstv_bA_nh_8
   namelist /dyn_fisl_p/ Cstv_bA_nh_8

   !# T* basic state temperature (K)
   real(kind=REAL64) :: Cstv_Tstr_8 = 240.0
   namelist /dyn_fisl  / Cstv_Tstr_8
   namelist /dyn_fisl_p/ Cstv_Tstr_8

   !# Inverse of PHI* basic state geopotential (m**2/s**2) used in GEM-H autobar
   real(kind=REAL64) :: Cstv_invFI_8

   !# Fraction of adjustment to be given to the ocean
   real(kind=REAL64) :: Cstv_psadj_8 = 1.d0
   namelist /dyn_fisl  / Cstv_psadj_8
   namelist /dyn_fisl_p/ Cstv_psadj_8

   !Schm

   !# True-> horizontal diffusion of momentum at each CN iteration
   logical :: Schm_hzdadw_L = .false.
   namelist /dyn_fisl  / Schm_hzdadw_L
   namelist /dyn_fisl_p/ Schm_hzdadw_L

   !# Number of iterations for Crank-Nicholson
   integer :: Schm_itcn = 2
   namelist /dyn_fisl  / Schm_itcn
   namelist /dyn_fisl_p/ Schm_itcn

   !# Number of iterations to solve non-linear Helmholtz problem
   integer :: Schm_itnlh  = 2
   namelist /dyn_fisl  / Schm_itnlh
   namelist /dyn_fisl_p/ Schm_itnlh

   !# *  -1: no blending between Yin and Yang
   !# *   0: blending at init only
   !# * > 0: blending at every nblendyy timestep
   integer :: Schm_nblendyy = -1
   namelist /dyn_fisl  / Schm_nblendyy
   namelist /dyn_fisl_p/ Schm_nblendyy

   !# * 0 -> No conservation of surface pressure
   !# * 1 -> Conservation of Total air mass Pressure
   !# * 2 -> Conservation of Dry air mass Pressure
   integer :: Schm_psadj = 0
   namelist /dyn_fisl  / Schm_psadj
   namelist /dyn_fisl_p/ Schm_psadj

   !# True-> print dry/wet air masses
   logical :: Schm_psadj_print_L = .false.
   namelist /dyn_fisl  / Schm_psadj_print_L
   namelist /dyn_fisl_p/ Schm_psadj_print_L

   !# True-> Tracers are mixing ratios with respect to dry air mass
   logical :: Schm_dry_mixing_ratio_L = .false.
   namelist /dyn_fisl  / Schm_dry_mixing_ratio_L
   namelist /dyn_fisl_P/ Schm_dry_mixing_ratio_L

   !# True-> use SLEVE vertical coordinate
   logical :: Schm_sleve_L = .false.

   !# True-> to use topography
   logical :: Schm_Topo_L = .true.
   namelist /dyn_fisl  / Schm_Topo_L
   namelist /dyn_fisl_p/ Schm_Topo_L

   !# * 0   ->          NO advection
   !# * 1   -> traditional advection
   !# * 2   -> consistent advection with respect to off-centering
   !# * 3   -> reversed advection with respect to off-centering
   integer :: Schm_advec = 1
   namelist /dyn_fisl  / Schm_advec
   namelist /dyn_fisl_p/ Schm_advec

   !# True-> Modify slightly code behaviour to ensure bitpattern
   !# reproduction in restart mode using FST file
   logical :: Schm_bitpattern_L = .false.
   namelist /dyn_fisl/ Schm_bitpattern_L

   !# Apply water loading in the calculations
   logical :: Schm_wload_L = .false.
   namelist /dyn_fisl  / Schm_wload_L
   namelist /dyn_fisl_p/ Schm_wload_L

   !# Physics coupling strategy
   character(len=16) :: Schm_phycpl_S = 'split'
   namelist /dyn_fisl  / Schm_phycpl_S
   namelist /dyn_fisl_p/ Schm_phycpl_S

  !Sol

   !# Type of solver
   !# * 'ITERATIF'
   !# * 'DIRECT'
   character(len=26) :: sol_type_S = 'DIRECT'
   namelist /dyn_fisl  / Sol_type_S
   namelist /dyn_fisl_p/ Sol_type_S

   !# Epsilon convergence criteria for none Yin-Yang iterative solver
   real(kind=REAL64) :: sol_fgm_eps   = 1.d-07
   namelist /dyn_fisl  / Sol_fgm_eps
   namelist /dyn_fisl_p/ Sol_fgm_eps

   !# Epsilon convergence criteria for the Yin-Yang iterative solver
   real(kind=REAL64) :: sol_yyg_eps   = 1.d-04
   namelist /dyn_fisl  / Sol_yyg_eps
   namelist /dyn_fisl_p/ Sol_yyg_eps

   !# maximum number of iterations allowed for none Yin-Yang iterative solver
   integer :: sol_fgm_maxits= 200
   namelist /dyn_fisl  / Sol_fgm_maxits
   namelist /dyn_fisl_p/ Sol_fgm_maxits

   !# maximum number of iterations allowed for the Yin-Yang iterative solver
   integer :: sol_yyg_maxits= 40
   namelist /dyn_fisl  / Sol_yyg_maxits
   namelist /dyn_fisl_p/ Sol_yyg_maxits

   !# size of Krylov subspace in iterative solver - should not exceed 100
   integer :: sol_im = 15
   namelist /dyn_fisl  / Sol_im
   namelist /dyn_fisl_p/ Sol_im

   !# 2D preconditioner for iterative solver
   character(len=26) :: sol2D_precond_S = 'JACOBI'
   namelist /dyn_fisl  / Sol2D_precond_S
   namelist /dyn_fisl_p/ Sol2D_precond_S

   !# 3D preconditioner for iterative solver
   character(len=26) :: sol3D_precond_S = 'JACOBI'
   namelist /dyn_fisl  / Sol3D_precond_S
   namelist /dyn_fisl_p/ Sol3D_precond_S

   !# 3D preconditioner for iterative solver
   integer :: Gauss_Niter = 2
   namelist /dyn_fisl  / Gauss_Niter
   namelist /dyn_fisl_p/ Gauss_Niter

   !# Krylov method for 3d iterative solver (FGMRES or FBICGSTAB)
   character(len=26) :: Sol3D_krylov_S = 'FGMRES'
   namelist /dyn_fisl  / Sol3D_krylov_S
   namelist /dyn_fisl_p/ Sol3D_krylov_S
      
   !# True => use the one-transpose solver
   logical :: Sol_one_transpose_L = .true.
   namelist /dyn_fisl  / Sol_one_transpose_L
   namelist /dyn_fisl_p/ Sol_one_transpose_L
      
   logical Schm_opentop_L
   integer Schm_nith

contains

!**s/r dyn_fisl_nml - Read namelist dyn_fisl

      integer function dyn_fisl_nml (F_unf)
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
#include <clib_interface_mu.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer :: err
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         dyn_fisl_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=dyn_fisl_p)
            if ( F_unf == -2 ) write (Lun_out,nml=dyn_fisl)
         end if
         return
      end if

      dyn_fisl_nml= -1 ; nml_must= .true. ; nml_S= 'dyn_fisl'

      rewind(F_unf)
      read (F_unf, nml=dyn_fisl, end= 1001, err=1003)
      dyn_fisl_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         dyn_fisl_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (dyn_fisl_nml < 0 ) return
      if ((Lun_out>=0).and.(dyn_fisl_nml==0)) write (Lun_out, 6004) trim(nml_S)
      dyn_fisl_nml= 1

      err = clib_toupper(Schm_phycpl_S)

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
      end function dyn_fisl_nml

   function dyn_fisl_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function dyn_fisl_options_init

end module dyn_fisl_options
