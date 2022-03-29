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
module init_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   integer, parameter :: IAU_MAX_TRACERS = 250

   !Iau

   !# Filter cutoff period for Iau_weight_S='sin' in hours
   real :: Iau_cutoff    = 6.
   namelist /init/ Iau_cutoff
   namelist /init_p/ Iau_cutoff

   !# The number of seconds between increment fields
   real :: Iau_interval = -1.
   namelist /init/ Iau_interval
   namelist /init_p/ Iau_interval

   !# The number of seconds over which IAU will be  will be run
   !# (typically the length of the assimilation window).
   !# Default < 0 means that no IAUs are applied.
   real :: Iau_period = -1.
   namelist /init/ Iau_period
   namelist /init_p/ Iau_period

   !# An optional list of tracers to be incremented.
   character(len=4), dimension(IAU_MAX_TRACERS) :: Iau_tracers_S = ' '
   namelist /init/ Iau_tracers_S

   !# The type of weighting function to be applied to the analysis increments:
   !# * 'constant' (default) uniform increments
   !# * 'sin' DF-style weights (Fillion et al. 1995)
   character(len=64) :: Iau_weight_S = 'constant'
   namelist /init/ Iau_weight_S
   namelist /init_p/ Iau_weight_S

   !# IAU Input PE blocking along npex
   integer :: Iau_ninblocx = 1
   namelist /init/ Iau_ninblocx
   namelist /init_p/ Iau_ninblocx

   !# IAU Input PE blocking along npey
   integer :: Iau_ninblocy = 1
   namelist /init/ Iau_ninblocy
   namelist /init_p/ Iau_ninblocy

   !# IAU Input Stats
   logical :: Iau_stats_L = .false.
   namelist /init/ Iau_stats_L
   namelist /init_p/ Iau_stats_L

!Init

   !# true -> Digital filter initialization is performed
   logical :: Init_balgm_L   = .false.
   namelist /init  / Init_balgm_L
   namelist /init_p/ Init_balgm_L

   !# true -> Windowing is applied
   logical :: Init_dfwin_L   = .true.
   namelist /init  / Init_dfwin_L
   namelist /init_p/ Init_dfwin_L

   !# number of points for digital filter (equals the number of timesteps +1)
   character(len=16) :: Init_dflength_S = '5p'
   namelist /init  / Init_dflength_S
   namelist /init_p/ Init_dflength_S

   !# period limit of digital filter units D,H,M,S
   character(len=16) :: Init_dfpl_S = '6h'
   namelist /init/ Init_dfpl_S
   namelist /init_p/ Init_dfpl_S

   !# * true -> passive tracers digitally filtered
   !# * false-> passive tracers set to result obtained at mid-period during initialization period (no filtering)
   logical :: Init_dftr_L = .false.
   namelist /init  / Init_dftr_L
   namelist /init_p/ Init_dftr_L

   !# True-> divergence high level modulation in initial computation of Zdot
   logical :: Zdot_divHLM_L = .false.
   namelist /init  / Zdot_divHLM_L
   namelist /init_p/ Zdot_divHLM_L

   logical Init_mode_L
   integer Init_dfnp, Init_halfspan
   real, dimension(:  ), allocatable :: Init_dfco
   real(kind=REAL64) Init_dfpl_8

contains

!**s/r init_nml - Read namelist init

      integer function init_nml (F_unf)
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
         init_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=init_p)
            if ( F_unf == -2 ) write (Lun_out,nml=init)
         end if
         return
      end if

      init_nml= -1 ; nml_must= .false. ; nml_S= 'init'

      rewind(F_unf)
      read (F_unf, nml=init, end= 1001, err=1003)
      init_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         init_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (init_nml < 0 ) return
      if ((Lun_out>=0).and.(init_nml==0)) write (Lun_out, 6004) trim(nml_S)
      init_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function init_nml
   function init_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function init_options_init

end module init_options
