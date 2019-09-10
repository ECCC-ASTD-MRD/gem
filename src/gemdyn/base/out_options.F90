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
module out_options
   implicit none
   public
   save

   integer, parameter :: MAXELEM_mod = 60

   !Out3

   !# True-> to clip humidity variables on output
   logical :: Out3_cliph_L = .false.
   namelist /out  / Out3_cliph_L
   namelist /out_p/ Out3_cliph_L

   !# Interval of output file name change
   character(len=16) :: Out3_close_interval_S = '1h'
   namelist /out  / Out3_close_interval_S
   namelist /out_p/ Out3_close_interval_S

   !# 'etiket' used for output fields
   character(len=12) :: Out3_etik_S = 'GEMDM'
   namelist /out  / Out3_etik_S
   namelist /out_p/ Out3_etik_S

   !# Vertical interpolation scheme for output
   character(len=12) :: Out3_vinterp_type_S = 'linear'
   namelist /out  / Out3_vinterp_type_S
   namelist /out_p/ Out3_vinterp_type_S

   !# Default value for IP3 is 0, -1 for IP3 to contain step number,
   !# >0 for given IP3
   integer :: Out3_ip3 = 0
   namelist /out  / Out3_ip3
   namelist /out_p/ Out3_ip3

   !# Number of layers close to the bottom of the model within which a
   !# linear interpolation of GZ will be performed
   integer :: Out3_linbot = 0
   namelist /out  / Out3_linbot
   namelist /out_p/ Out3_linbot

   !# Packing factor used for all variables except for those defined in
   !# Out_xnbits_s
   integer :: Out3_nbitg = 16
   namelist /out  / Out3_nbitg
   namelist /out_p/ Out3_nbitg

   !# Minimum of digits used to represent output units
   integer :: Out3_ndigits = 3
   namelist /out  / Out3_ndigits
   namelist /out_p/ Out3_ndigits

   !# List of levels for underground extrapolation
   real, dimension(MAXELEM_mod) :: Out3_lieb_levels = 0.
   namelist /out/ Out3_lieb_levels

   !# Maximum number of iterations for the Liebman procedure
   integer :: Out3_lieb_maxite = 100
   namelist /out  / Out3_lieb_maxite
   namelist /out_p/ Out3_lieb_maxite

   !# number of iterations to exchange halo for the Liebman procedure
   integer :: Out3_liebxch_iter = 4
   namelist /out  / Out3_liebxch_iter
   namelist /out_p/ Out3_liebxch_iter

   !# Precision criteria for the Liebman procedure
   real :: Out3_lieb_conv = 0.1
   namelist /out  / Out3_lieb_conv
   namelist /out_p/ Out3_lieb_conv

   !# Sortie jobs lauched every Out3_postproc_fact*Out3_close_interval_S
   integer :: Out3_postproc_fact = 6
   namelist /out  / Out3_postproc_fact
   namelist /out_p/ Out3_postproc_fact

   !# Total number of PEs for output using MFV collector
   integer :: Out3_npes = 1
   namelist /out/ Out3_npes
   namelist /out_p/ Out3_npes

   !# Total number of PEs along npex for output using MID collector
   integer :: Out3_npex = -1
   namelist /out  / Out3_npex
   namelist /out_p/ Out3_npex

   !# Total number of PEs along npey for output using MID collector
   integer :: Out3_npey = -1
   namelist /out  / Out3_npey
   namelist /out_p/ Out3_npey

contains

!**s/r out_nml - Read namelist out

      integer function out_nml (F_unf)
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
         out_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=out_p)
            if ( F_unf == -2 ) write (Lun_out,nml=out)
         end if
         return
      end if

      out_nml= -1 ; nml_must= .false. ; nml_S= 'out'

      rewind(F_unf)
      read (F_unf, nml=out, end= 1001, err=1003)
      out_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         out_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (out_nml < 0 ) return
      if ((Lun_out>=0).and.(out_nml==0)) write (Lun_out, 6004) trim(nml_S)
      out_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function out_nml

   function out_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function out_options_init

end module out_options
