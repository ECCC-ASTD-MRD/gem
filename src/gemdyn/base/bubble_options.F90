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
module bubble_options
   use dynkernel_options
   use HORgrid_options
   use dyn_fisl_options
   use VERgrid_options
   use lun
   use glb_ld
   use cstv
   use dcst
   use tdpack
   implicit none
   public
   save

   !#
   integer :: bubble_ni = 101
   namelist /bubble_cfgs/ bubble_ni
   !#
   integer :: bubble_nj = 1
   namelist /bubble_cfgs/ bubble_nj
   !#
   integer :: bubble_nk = 100
   namelist /bubble_cfgs/ bubble_nk
   !#
   real :: bubble_dx = 10.
   namelist /bubble_cfgs/ bubble_dx
   !#
   real :: bubble_dz = 10.
   namelist /bubble_cfgs/ bubble_dz
   !#
   real :: bubble_theta = 303.16
   namelist /bubble_cfgs/  bubble_theta
   !#
   integer :: bubble_rad = 25
   namelist /bubble_cfgs/ bubble_rad
   !#
   integer :: bubble_ictr = -1
   namelist /bubble_cfgs/ bubble_ictr
   !#
   integer :: bubble_kctr = -1
   namelist /bubble_cfgs/ bubble_kctr
   !#
   logical :: bubble_gaus_L = .false.
   namelist /bubble_cfgs/ bubble_gaus_L

contains

      integer function bubble_nml (F_unf)
      implicit none

      integer F_unf

      logical nml_must
      character(len=64) :: nml_S
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         bubble_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=bubble_cfgs)
         return
      end if

      bubble_nml= -1 ; nml_must= .true. ; nml_S= 'bubble_cfgs'

      rewind(F_unf)
      read (F_unf, nml=bubble_cfgs, end= 1001, err=1003)
      bubble_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         bubble_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (bubble_nml < 0 ) return
      if ((Lun_out>=0).and.(bubble_nml==0)) write (Lun_out, 6004) trim(nml_S)

      ! establish horizontal grid configuration
      ! (must absolutely be done here)
      Dcst_rayt_8 = Dcst_rayt_8*0.1d0 ! an accuracy problem
      Dcst_inv_rayt_8 = Dcst_inv_rayt_8 * 10.d0 ! an accuracy problem
      Grd_typ_S='LU'
      Grd_ni = bubble_ni ; Grd_nj = bubble_nj
      Grd_dx = (bubble_dx/Dcst_rayt_8)*(180./pi_8)
      Grd_dy = Grd_dx
      Grd_latr = 0.
      Grd_lonr = (bubble_ni/2 + 20) * Grd_dx
      Grd_maxcfl = 3

      bubble_nml=0

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      return
      end function bubble_nml

end module bubble_options
