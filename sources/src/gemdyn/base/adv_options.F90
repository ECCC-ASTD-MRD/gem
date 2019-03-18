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
module adv_options
   implicit none
   public
   save

   !# True-> MONO(CLIPPING) of RHS
   logical :: adv_rhst_mono_L = .false.
   namelist /adv_cfgs/ adv_rhst_mono_L

   !# List of available P exponents REAL (size PEXP_LIST_MAX=9)
   integer, parameter :: PEXP_LIST_MAX = 9
   real :: adv_BC_pexp_list(PEXP_LIST_MAX) = 1.0
   namelist /adv_cfgs/ adv_BC_pexp_list

   !# True-> Clipping Min Max + Proportional Mass-Fixer after Bermejo-Conde
   logical :: adv_BC_min_max_L = .true.
   namelist /adv_cfgs/ adv_BC_min_max_L

   !# Type of Flux when Bermejo-Conde LAM
   !# * 1 -> Aranami et al. (2015)
   !# * 2 -> Zerroukat and Shipway (2017) (ZLF)
   integer :: Adv_BC_LAM_flux = 1
   namelist /adv_cfgs/ adv_BC_LAM_flux

   !# Number of neighborhood zones in ILMC
   integer :: adv_ILMC_sweep_max = 2
   namelist /adv_cfgs/ adv_ILMC_sweep_max

   !# True-> Clipping Min Max after ILMC
   logical :: adv_ILMC_min_max_L = .true.
   namelist /adv_cfgs/ adv_ILMC_min_max_L

   !# Activate printing of conservation diagnostics if /=0
   integer :: adv_verbose = 0
   namelist /adv_cfgs/ adv_verbose

   !# True-> Bermejo-Conde LEGACY=T for current tracer
   !# * T -> mass=1
   !# * F -> mass=100+weight*10+pexp_n
   logical :: adv_BC_LEGACY_L

   !# Bermejo-Conde Weight for current tracer
   !# * 1 -> Additive
   !# * 2 -> Multiplicative
   !# * 3 -> Additive+Factor pr_k/pr_s)
   integer :: adv_BC_weight

   !# P exponent NUMBER for current tracer corresponding to
   !# rank in adv_BC_pexp_list(PEXP_LIST_MAX) to choose P exponent REAL
   integer :: adv_BC_pexp_n

   !# Bermejo-Conde LAM: Components
   !# * 0 -> Advection only
   !# * 1 -> Advection+Flux
   !# * 2 -> Flux only
   integer :: adv_BC_LAM_flux_n = 0

   !# Bermejo-Conde LAM: Pointers mask_o/mask_i
   real, dimension(:,:,:), allocatable :: Adv_BC_LAM_mask_o, Adv_BC_LAM_mask_i

   !# Bermejo-Conde LAM: Factor applied to Grd_maxcfl in Pilot area
   !# * 1 -> Default
   !# * 2 -> Bermejo-Conde Flux ZLF
   integer :: adv_maxcfl_fact = 1

   !# True-> Mass Conservation/Shape-preserving for at least one tracer
   logical :: adv_Mass_Cons_L = .false.

   !# True-> Mass Conservation/Shape-preserving for current tracer
   logical :: adv_Mass_Cons_tr_L = .false.

   !# True-> Do extended advection operations
   logical :: Adv_extension_L = .false.

   !# Reset upstream positions at each timestep if 1
   !# * Dimension 1 for Bermejo-Conde LAM Aranami
   !# * Dimension 2 for LCSL advection
   integer :: Adv_pos_reset(2) = 0

   !# South boundary in GY for an embedded LAM
   integer :: adv_pil_sub_s = -1
   namelist /adv_cfgs/ adv_pil_sub_s

   !# North boundary in GY for an embedded LAM
   integer :: adv_pil_sub_n = -1
   namelist /adv_cfgs/ adv_pil_sub_n

   !# West boundary in GY for an embedded LAM
   integer :: adv_pil_sub_w = -1
   namelist /adv_cfgs/ adv_pil_sub_w

   !# East boundary in GY for an embedded LAM
   integer :: adv_pil_sub_e = -1
   namelist /adv_cfgs/ adv_pil_sub_e

   !# Core/Subset areas
   real*8 :: adv_gc_area_8,adv_gs_area_8

   !# True-> LCSL advection for at least one tracer
   logical :: adv_LCSL_option_L = .false.

   !# Type of LCSL advection for current tracer
   integer :: adv_LCSL_n = 0

   !# Type of Monotonic LCSL advection for current tracer
   integer :: adv_LCSL_mono_n = 0

   !# Variable for NONE/Cubic/Quintic Interpolation
   character(len=12) adv_intp_S

contains

!**s/r adv_nml - Read namelist adv

      integer function adv_nml (F_unf)
      use adv_grid
      use HORgrid_options
      use lun
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
         adv_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=adv_cfgs)
         return
      end if

      adv_nml= -1 ; nml_must= .false. ; nml_S= 'adv_cfgs'

      rewind(F_unf)
      read (F_unf, nml=adv_cfgs, end= 1001, err=1003)
      adv_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         adv_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (adv_nml < 0 ) return
      if ((Lun_out>=0).and.(adv_nml==0)) write (Lun_out, 6004) trim(nml_S)
      adv_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      adv_maxcfl= max(1,Grd_maxcfl)
      adv_halox = adv_maxcfl + 1
      adv_haloy = adv_halox
      if (Adv_BC_LAM_flux==2) adv_maxcfl_fact = 2
!
!-------------------------------------------------------------------
!
      return
      end function adv_nml

   function adv_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function adv_options_init

end module adv_options
