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
module grdc_options
   implicit none
   public
   save

   integer, parameter :: MAX_TRNM = 1000

   !# x horizontal resolution of target cascade grid (degrees)
   real  :: Grdc_dx = -1.
   namelist /grdc/ Grdc_dx
   namelist /grdc_p/ Grdc_dx

   !# y horizontal resolution of target cascade grid (degrees)
   real  :: Grdc_dy = -1.
   namelist /grdc/ Grdc_dy
   namelist /grdc_p/ Grdc_dy

   !# Latitude on rotated grid of ref point, Grdc_iref,Grdc_jref (degrees)
   real  :: Grdc_latr = 0.
   namelist /grdc/ Grdc_latr
   namelist /grdc_p/ Grdc_latr

   !# Longitude on rotated grid of ref point, Grdc_iref,Grdc_jref (degrees)
   real  :: Grdc_lonr = 180.
   namelist /grdc/ Grdc_lonr
   namelist /grdc_p/ Grdc_lonr

   !# Reference Point I on rotated cascade grid, 1 < Grdc_iref < Grdc_ni
   integer :: Grdc_iref = -1
   namelist /grdc/ Grdc_iref
   namelist /grdc_p/ Grdc_iref

   !# Reference Point J on rotated cascade grid, 1 < Grdc_jref < Grdc_nj
   integer :: Grdc_jref = -1
   namelist /grdc/ Grdc_jref
   namelist /grdc_p/ Grdc_jref

   !# Max Supported Courrant number;
   !# Pilot area=Grdc_maxcfl +Grdc_bsc_base+Grdc_bsc_ext1
   integer :: Grdc_maxcfl = 1
   namelist /grdc/ Grdc_maxcfl
   namelist /grdc_p/ Grdc_maxcfl

   !# Number of bits for the packing factor
   integer, dimension(2) :: Grdc_nbits = [ 16, 12 ]
   namelist /grdc/ Grdc_nbits
   namelist /grdc_p/ Grdc_nbits

   !# TRUE to dump out permanent bus for cascade mode
   logical :: Grdc_fullgrid_L = .false.
   namelist /grdc/ Grdc_fullgrid_L
   namelist /grdc_p/ Grdc_fullgrid_L

   !# Nesting interval specified with digits ending with one character
   !# for the units:
   !# * S : seconds
   !# * D : days
   !# * M : minutes
   !# * H : hours
   character(len=15) :: Grdc_nfe = ' '
   namelist /grdc/ Grdc_nfe
   namelist /grdc_p/ Grdc_nfe

   !# Number of points along X
   integer :: Grdc_ni = 0
   namelist /grdc/ Grdc_ni
   namelist /grdc_p/ Grdc_ni

   !# Number of points along Y
   integer :: Grdc_nj = 0
   namelist /grdc/ Grdc_nj
   namelist /grdc_p/ Grdc_nj

   !# Time string (units D, H, M or S) from the start of the run to
   !# start producing the cascade files
   character(len=15) :: Grdc_start_S = ' '
   namelist /grdc/ Grdc_start_S
   namelist /grdc_p/ Grdc_start_S

   !# Time string (units D, H, M or S) from the start of the run to
   !# stop producing the cascade files
   character(len=15) :: Grdc_end_S = ' '
   namelist /grdc/ Grdc_end_S
   namelist /grdc_p/ Grdc_end_S

   !# List of tracers to be written from piloting run
   character(len=4) :: Grdc_trnm_S(MAX_TRNM) = '@#$%'
   namelist /grdc/ Grdc_trnm_S

   integer Grdc_gid,Grdc_gif,Grdc_gjd,Grdc_gjf, &
           Grdc_ndt,Grdc_ntr,Grdc_start,Grdc_end

contains

!**s/r grdc_nml - Read namelist grdc

      integer function grdc_nml (F_unf)
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
         grdc_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=grdc_p)
            if ( F_unf == -2 ) write (Lun_out,nml=grdc)
         end if
         return
      end if

      grdc_nml= -1 ; nml_must= .false. ; nml_S= 'grdc'

      rewind(F_unf)
      read (F_unf, nml=grdc, end= 1001, err=1003)
      grdc_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         grdc_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (grdc_nml < 0 ) return
      if ((Lun_out>=0).and.(grdc_nml==0)) write (Lun_out, 6004) trim(nml_S)
      grdc_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function grdc_nml

   function grdc_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function grdc_options_init

end module grdc_options
