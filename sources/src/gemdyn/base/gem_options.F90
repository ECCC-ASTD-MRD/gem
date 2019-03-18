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
module gem_options
   implicit none
   public
   save

   integer,parameter  :: STAT_MAXN = 500

   !# number of points for the halo on X
   integer :: G_halox = 4
   namelist /gem_cfgs  / G_halox
   namelist /gem_cfgs_p/ G_halox

   !# number of points for the halo on Y
   integer :: G_haloy = 4
   namelist /gem_cfgs  / G_haloy
   namelist /gem_cfgs_p/ G_haloy

   !# Heap memory will be painted to NaN using an array wrk01(G_ni,G_nj,Heap_nk)
   integer :: Heap_nk = -1
   namelist /gem_cfgs  / Heap_nk
   namelist /gem_cfgs_p/ Heap_nk

   !# True->to check for time left in job
   logical :: Lctl_cktimeleft_L = .false.
   namelist /gem_cfgs  / Lctl_cktimeleft_L
   namelist /gem_cfgs_p/ Lctl_cktimeleft_L

   !# True->to print more information to std output
   logical :: Lctl_debug_L = .false.
   namelist /gem_cfgs  / Lctl_debug_L
   namelist /gem_cfgs_p/ Lctl_debug_L

   !# precision in print glbstats
   !# * 'LCL_4'
   !# * 'GLB_8'
   character(len=6) :: Lctl_rxstat_S = 'LCL_4'
   namelist /gem_cfgs  / Lctl_rxstat_S
   namelist /gem_cfgs_p/ Lctl_rxstat_S

   !# list of variables to do blocstat.
   !# Any gmm variable name, or predefine lists :
   !# * 'ALL'
   !# * 'ALL_DYN_T0'
   !# * 'ALL_DYN_T1'
   !# * 'ALL_TR_T0'
   !# * 'ALL_TR_T1'
   character(len=32) :: stat_liste(STAT_MAXN) = ' '
   namelist /gem_cfgs/ stat_liste

   !# List of tracers and parameters to be read from analyse
   !# * hzd = Horizontal Diffusion if 1
   !# * wload = Water loading if 1
   !# * min = Minimum admissible value
   !# * max = Maximum admissible value
   !# * mono = Type of monotonicity
   !# ** 0 -> None
   !# ** 1 -> Clipping QMSL
   !# ** 2 -> ILMC
   !# * mass = Type of mass conservation
   !# ** 0 -> None
   !# ** 1 -> Bermejo-Conde LEGACY=T
   !# ** 100+weight*10+pexp_n -> Bermejo-Conde LEGACY=F
   !# * intp  = Type of interpolation in advection (NOT case sensitive)
   !# ** none -> No advection
   !# ** cubic -> 3D Cubic Lagrange
   !# ** quintic -> Horizontal Cubic + Vertical Quintic
   !# List of sub-parameters
   !# * weight = Weight in Bermejo-Conde
   !# ** 1 -> Additive
   !# ** 2 -> Multiplicative
   !# ** 3 -> Additive+Factor pr_k/pr_s
   !# * pexp_n = Rank in Adv_BC_pexp_list(PEXP_LIST_MAX) to choose P exponent REAL
   character(len=512), dimension(500) :: Tr3d_list_S = ' '
   namelist /gem_cfgs/ Tr3d_list_S

   !# Override for default tracers attributes
   character(len=512) :: Tr3d_default_s = ' '
   namelist /gem_cfgs/ Tr3d_default_s

   !# True-> tracers validity time does not have to match analysis
   logical :: Tr3d_anydate_L= .false.
   namelist /gem_cfgs  / Tr3d_anydate_L
   namelist /gem_cfgs_p/ Tr3d_anydate_L

   !Vtopo

   !# Time at which to start evolving topography toward target
   character(len=16) :: Vtopo_start_S = ''
   namelist /gem_cfgs  / Vtopo_start_S
   namelist /gem_cfgs_p/ Vtopo_start_S

   !# On which length of time to evolve topography
   character(len=16) :: Vtopo_length_S = ''
   namelist /gem_cfgs  / Vtopo_length_S
   namelist /gem_cfgs_p/ Vtopo_length_S

   logical Vtopo_L
   integer stat_nombre
   integer Vtopo_start, Vtopo_ndt

contains

!**s/r gem_nml - Read namelist gem

      integer function gem_nml (F_unf)
      use canonical
      use ctrl
      use grdc_options
      use lun
      implicit none
#include <arch_specific.hf>
#include <clib_interface_mu.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer err_grdc,err_canonical,err
!
!-------------------------------------------------------------------
!
      if ( F_unf < 0 ) then
         gem_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) then
               write (Lun_out,nml=gem_cfgs_p)
               if (Ctrl_canonical_dcmip_L .or. Ctrl_canonical_williamson_L) then
                  err= canonical_nml (-1,nml_must,nml_must)
               end if
            end if
            if ( F_unf == -2 ) write (Lun_out,nml=gem_cfgs)
         end if
         return
      end if

      gem_nml= -1 ; nml_must= .false. ; nml_S= 'gem_cfgs'

      rewind(F_unf)
      read (F_unf, nml=gem_cfgs, end= 1001, err=1003)
      gem_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         gem_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if ((Lun_out>=0).and.(gem_nml==0)) write (Lun_out, 6004) trim(nml_S)

      err_grdc     = grdc_nml (F_unf)
      err_canonical= canonical_nml (F_unf, Ctrl_canonical_dcmip_L,&
                                      Ctrl_canonical_williamson_L )
      gem_nml= min(gem_nml, err_grdc, err_canonical)

      err = clib_toupper(Lctl_rxstat_S)

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
!
!-------------------------------------------------------------------
!
      return
      end function gem_nml

   function gem_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function gem_options_init

end module gem_options
