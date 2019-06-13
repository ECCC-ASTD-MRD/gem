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
module dyn_expo_options
!temporaire
   use dyn_fisl_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Type of integrator
   !# * 'EPI2'
   !# * 'EPI3'
   character(len=16) :: Exp_integrator_S = 'EPI3'
   namelist /dyn_expo  / Exp_integrator_S
   namelist /dyn_expo_p/ Exp_integrator_S

   !# Tolerance to achieve in KIOPS
   real(kind=REAL64) :: Kiops_tolerance = 1.19209d-07
   namelist /dyn_expo  / Kiops_tolerance
   namelist /dyn_expo_p/ Kiops_tolerance

   !# An estimate of the appropriate Krylov size
   integer :: Kiops_krylov_size = 16
   namelist /dyn_expo  / Kiops_krylov_size
   namelist /dyn_expo_p/ Kiops_krylov_size

   !# The minimum Krylov size in adaptive procedure
   integer :: Kiops_krylov_size_min = 1
   namelist /dyn_expo  / Kiops_krylov_size_min
   namelist /dyn_expo_p/ Kiops_krylov_size_min

   !# The maximum Krylov size in adaptive procedure
   integer :: Kiops_krylov_size_max = 64
   namelist /dyn_expo  / Kiops_krylov_size_max
   namelist /dyn_expo_p/ Kiops_krylov_size_max

   !# Free parameter to control artificial diffusion in upwind advection scheme
   real :: adv_alpha_hor = 0.5
   namelist /dyn_expo  / adv_alpha_hor
   namelist /dyn_expo_p/ adv_alpha_hor

!temporaire
   namelist /dyn_expo  / Cstv_tstr_8
   namelist /dyn_expo_p/ Cstv_tstr_8

contains

!**s/r dyn_expo_nml - Read namelist dyn_expo

      integer function dyn_expo_nml (F_unf)
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
         dyn_expo_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=dyn_expo_p)
         return
      end if

      dyn_expo_nml= -1 ; nml_must= .false. ; nml_S= 'dyn_expo'

      rewind(F_unf)
      read (F_unf, nml=dyn_expo, end= 1001, err=1003)
      dyn_expo_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         dyn_expo_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (dyn_expo_nml < 0 ) return
      if ((Lun_out>=0).and.(dyn_expo_nml==0)) write (Lun_out, 6004) trim(nml_S)
      dyn_expo_nml= 1

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
      end function dyn_expo_nml

   function dyn_expo_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function dyn_expo_options_init

end module dyn_expo_options
