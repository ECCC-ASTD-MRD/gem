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
module spn_options
   implicit none
   public
   save

   !# Spectral nudging list of variables (eg. 'UVT' or 'UV')
   character(len=16) :: Spn_nudging_S = ' '
   namelist /spn  / Spn_nudging_S
   namelist /spn_p/ Spn_nudging_S

   !# Nudging profile lower end in hyb level (eg. 1.0 or 0.8)
   !# If use 0.8, the profile will be set zero when hyb > 0.8
   real :: Spn_start_lev = 1.0
   namelist /spn  / Spn_start_lev
   namelist /spn_p/ Spn_start_lev

   !# Nudging profile upper end in hyb level (eg. 0.0 or 0.2)
   !# If use 0.2, the profile wll be set 1.0 when hyb < 0.2
   real :: Spn_up_const_lev = 0.0
   namelist /spn  / Spn_up_const_lev
   namelist /spn_p/ Spn_up_const_lev

   !# Nudging profile transition shape('COS2' or 'LINEAR')
   !# Set the shape between Spn_start_lev and Spn_up_const_lev
   character(len=16) :: Spn_trans_shape_S = 'LINEAR'
   namelist /spn  / Spn_trans_shape_S
   namelist /spn_p/ Spn_trans_shape_S

   !# Nudging relaxation timescale (eg. 10 hours )
   real :: Spn_relax_hours = 10.
   namelist /spn  / Spn_relax_hours
   namelist /spn_p/ Spn_relax_hours

   !# The filter will be set zero for smaller scales (in km)
   real :: Spn_cutoff_scale_large = 300.
   namelist /spn  / Spn_cutoff_scale_large
   namelist /spn_p/ Spn_cutoff_scale_large

   !# The filter will be set 1.0 for larger scales (in km) between
   !# Spn_cutoff_scale_small and Spn_cutoff_scale_large,
   !# the filter will have a COS2 transition.
   real :: Spn_cutoff_scale_small = 100.
   namelist /spn  / Spn_cutoff_scale_small
   namelist /spn_p/ Spn_cutoff_scale_small

   !# Nudging interval in seconds (eg. 1800, means nudging is performed
   !# every every 30 minutes)
   integer :: Spn_step = 21600
   namelist /spn  / Spn_step
   namelist /spn_p/ Spn_step

   !# Nudging weight in temporal space (.true. or .false.).
   !# If the driving fields are available every 6 hours and Spn_step is
   !# set to 30 minutes then nudging will have more weight every six hours
   !# when the driving fields are available
   logical :: Spn_weight_L = .false.
   namelist /spn  / Spn_weight_L
   namelist /spn_p/ Spn_weight_L

   !# The weight factor when Spn_weight_L=.true.
   !# (The weigh factor is COS2**(Spn_wt_pwr), Spn_wt_pwr could  be set as
   !# 0, 2, 4, 6. If Spn_wt_pwr = 2, weight factor is COS2)
   integer :: Spn_wt_pwr = 2
   namelist /spn  / Spn_wt_pwr
   namelist /spn_p/ Spn_wt_pwr

contains

!**s/r spn_nml - Read namelist spn

      integer function spn_nml (F_unf)
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
         spn_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=spn_p)
            if ( F_unf == -2 ) write (Lun_out,nml=spn)
         end if
         return
      end if

      spn_nml= -1 ; nml_must= .false. ; nml_S= 'spn'

      rewind(F_unf)
      read (F_unf, nml=spn, end= 1001, err=1003)
      spn_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         spn_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (spn_nml < 0 ) return
      if ((Lun_out>=0).and.(spn_nml==0)) write (Lun_out, 6004) trim(nml_S)
      spn_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function spn_nml

   function spn_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function spn_options_init

end module spn_options
