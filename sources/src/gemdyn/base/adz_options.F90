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
module adz_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Use surface layer winds for advection of lowest thermodynamic level
   logical :: Adz_slt_winds = .false.
   namelist /Adz_cfgs/ Adz_slt_winds

   !# Number of iterations for trajectories computation
   integer :: Adz_itraj = 3
   namelist /adz_cfgs/ Adz_itraj

   integer Adz_niter

   !# List of available P exponents REAL (size PEXP_LIST_MAX=9)
   integer, parameter :: PEXP_LIST_MAX = 9
   real :: adz_BC_pexp_list(PEXP_LIST_MAX) = 1.0
   namelist /adz_cfgs/ adz_BC_pexp_list

   !# True-> Clipping Min Max + Proportional Mass-Fixer after Bermejo-Conde
   logical :: adz_BC_min_max_L = .true.
   namelist /adz_cfgs/ adz_BC_min_max_L

   !# Type of Flux when Bermejo-Conde LAM
   !# * 1 -> Aranami et al. (2015)
   !# * 2 -> Zerroukat and Shipway (2017) (ZLF)
   integer :: Adz_BC_LAM_flux = 1
   namelist /adz_cfgs/ adz_BC_LAM_flux

   !# Number of neighborhood zones in ILMC
   integer :: adz_ILMC_sweep_max = 2
   namelist /adz_cfgs/ adz_ILMC_sweep_max

   !# True-> Clipping Min Max after ILMC
   logical :: adz_ILMC_min_max_L = .true.
   namelist /adz_cfgs/ adz_ILMC_min_max_L

   !# Activate printing of conservation diagnostics if /=0
   integer :: adz_verbose = 0
   namelist /adz_cfgs/ adz_verbose

   !# True-> Bermejo-Conde LEGACY=T for current tracer
   !# * T -> mass=1
   !# * F -> mass=100+weight*10+pexp_n
   logical :: adz_BC_LEGACY_L

   !# Bermejo-Conde Weight for current tracer
   !# * 1 -> Additive
   !# * 2 -> Multiplicative
   !# * 3 -> Additive+Factor pr_k/pr_s)
   integer :: adz_BC_weight

   !# P exponent NUMBER for current tracer corresponding to
   !# rank in adz_BC_pexp_list(PEXP_LIST_MAX) to choose P exponent REAL
   integer :: adz_BC_pexp_n

   !# Bermejo-Conde LAM: Components
   !# * 0 -> Advection only
   !# * 1 -> Advection+Flux
   !# * 2 -> Flux only
   integer :: adz_BC_LAM_flux_n = 0

   !# Bermejo-Conde LAM: Pointers mask_o/mask_i
   real, dimension(:,:,:), allocatable :: Adz_BC_LAM_mask_o, Adz_BC_LAM_mask_i

   !# Bermejo-Conde LAM: Factor applied to Grd_maxcfl in Pilot area
   !# * 1 -> Default
   !# * 2 -> Bermejo-Conde Flux ZLF
   integer :: adz_maxcfl_fact = 1

   !# True-> Mass Conservation/Shape-preserving for at least one tracer
   logical :: adz_Mass_Cons_L = .false.

   !# True-> Mass Conservation/Shape-preserving for current tracer
   logical :: adz_Mass_Cons_tr_L = .false.

   !# Bermejo-Conde LAM Aranami: Computations related to upstream positions done at each timestep
   !# * 0 -> No
   !# * 1 -> Yes
   integer :: adz_pos_reset = 0

   !# South boundary in GY for an embedded LAM
   integer :: adz_pil_sub_s = -1
   namelist /adz_cfgs/ adz_pil_sub_s

   !# North boundary in GY for an embedded LAM
   integer :: adz_pil_sub_n = -1
   namelist /adz_cfgs/ adz_pil_sub_n

   !# West boundary in GY for an embedded LAM
   integer :: adz_pil_sub_w = -1
   namelist /adz_cfgs/ adz_pil_sub_w

   !# East boundary in GY for an embedded LAM
   integer :: adz_pil_sub_e = -1
   namelist /adz_cfgs/ adz_pil_sub_e

   !# Core/Subset areas
   real(kind=REAL64) :: adz_gc_area_8,adz_gs_area_8

   !# Variable for NONE/Cubic/Quintic Interpolation
   character(len=12) :: adz_intp_S = 'CUBIC'

   integer, pointer, contiguous, dimension(:) :: ii_w

contains

!**s/r adz_nml - Read namelist adz

      integer function adz_nml (F_unf)
      use adv_grid
      use HORgrid_options
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         adz_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=adz_cfgs)
         end if
         return
      end if

      adz_nml= -1 ; nml_must= .false. ; nml_S= 'adz_cfgs'

      rewind(F_unf)
      read (F_unf, nml=adz_cfgs, end= 1001, err=1003)
      adz_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         adz_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (adz_nml < 0 ) return
      if ((Lun_out>=0).and.(adz_nml==0)) write (Lun_out, 6004) trim(nml_S)
      adz_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      adv_maxcfl= max(1,Grd_maxcfl)
      adv_halox = adv_maxcfl + 1
      adv_haloy = adv_halox
      if (Adz_BC_LAM_flux==2) adz_maxcfl_fact = 2
!
!-------------------------------------------------------------------
!
      return
      end function adz_nml

end module adz_options
