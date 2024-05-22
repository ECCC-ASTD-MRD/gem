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
      use, intrinsic :: iso_fortran_env
      implicit none
   public
   save

   !# Nudging profile lower end in hyb level (eg. 0.8) or hyb_H level (for GEM-H)
   !# profile will be set to 0. when hyb > Spn_start_lev
   real :: Spn_start_lev = 1.0
   namelist /spn  / Spn_start_lev
   namelist /spn_p/ Spn_start_lev

   !# Nudging profile: upper end for constant layer in hyb level (eg. 0.0 or 0.2) or hyb_H level (for GEM-H)
   !# profile wll be set to 1.0 when hyb is between Spn_const_lev_bot and Spn_const_lev_top
   real :: Spn_const_lev_top = 0.0
   namelist /spn  / Spn_const_lev_top
   namelist /spn_p/ Spn_const_lev_top

   !# Nudging profile:  end for constant layer in hyb level (eg. 0.0 or 0.2) or hyb_H level (for GEM-H)
   !# profile wll be set to 1.0 wehn hyb is between Spn_const_lev_bot and Spn_const_lev_top
   real :: Spn_const_lev_bot = 0.0
   namelist /spn  / Spn_const_lev_bot
   namelist /spn_p/ Spn_const_lev_bot

   !# Nudging profile upper end in hyb level (eg. 0.1) or hyb_H level (for GEM-H)
   !# profile will be set to 0. when hyb < Spn_end_lev
   real :: Spn_end_lev = 0.0
   namelist /spn  / Spn_end_lev
   namelist /spn_p/ Spn_end_lev
   

   !# Nudging profile transition shape('COS2' or 'LINEAR')
   !# Set the shape between Spn_start_lev and Spn_up_const_lev
   character(len=16) :: Spn_trans_shape_S = 'LINEAR'
   namelist /spn  / Spn_trans_shape_S
   namelist /spn_p/ Spn_trans_shape_S

   !# Nudging relaxation timescale - in hours
   real :: Spn_relax_hours = 10.
   namelist /spn  / Spn_relax_hours
   namelist /spn_p/ Spn_relax_hours

   !# The filter will be set to 0. for smaller scales (in km)
   real :: Spn_cutoff_scale_large = 300.
   namelist /spn  / Spn_cutoff_scale_large
   namelist /spn_p/ Spn_cutoff_scale_large

   !# The filter will be set to 1.0 for larger scales (in km)
   !# Transition between Spn_cutoff_scale_small and Spn_cutoff_scale_large
   !# will have a COS2 shape.
   real :: Spn_cutoff_scale_small = 100.
   namelist /spn  / Spn_cutoff_scale_small
   namelist /spn_p/ Spn_cutoff_scale_small

   !# Nudging interval - in sec
   !# Nudging is performed every every Spn_freq sec
   integer :: Spn_freq = -1
   namelist /spn  / Spn_freq
   namelist /spn_p/ Spn_freq

   !# Nudging weight in temporal space (.true. or .false.).
   !# If the driving fields are available every 6 hours and Spn_freq is
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

   !# Availability interval of nudging data for Global Yin-Yang - in sec
   !# Driving(nudging) data is available every Spn_yy_nudge_data_freq sec
   integer :: Spn_yy_nudge_data_freq = -1
   namelist /spn  / Spn_yy_nudge_data_freq
   namelist /spn_p/ Spn_yy_nudge_data_freq

   !# Nudging specific humidity (.true. or .false.).
   logical :: Spn_nudge_HU_L = .false.
   namelist /spn  / Spn_nudge_HU_L
   namelist /spn_p/ Spn_nudge_HU_L

   character(len=16) :: Spn_nudging_S = ' ' ! depricated
   logical :: Spn_ON_L = .false.
   integer :: Spn_12smin, Spn_12smax, Spn_12sn, Spn_12sn0
   integer :: Spn_22min , Spn_22max , Spn_22n , Spn_22n0
   integer :: Spn_22pil_w, Spn_22pil_e, Spn_interval, Spn_ws
   integer :: Spn_njnh,Spn_nk12,Spn_ni22
   real :: Spn_weight
   real, dimension(:), allocatable :: prof
   real(kind=REAL64) , dimension(:,:  ), allocatable :: Spn_flt
   real(kind=REAL64) , dimension(:,:,:), allocatable :: Spn_fft,&
                                               Spn_fdg, Spn_wrk

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

end module spn_options
