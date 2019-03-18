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
module lam_options
   implicit none
   public
   save

   !# Number of levels for top piloting
   integer :: Lam_gbpil_T = -1
   namelist /lam  / Lam_gbpil_T
   namelist /lam_p/ Lam_gbpil_T

   !# Number of points for horizontal blending
   integer :: Lam_blend_H = 10
   namelist /lam  / Lam_blend_H
   namelist /lam_p/ Lam_blend_H

   !# Number of levels for top blending
   integer :: Lam_blend_T = 0
   namelist /lam  / Lam_blend_T
   namelist /lam_p/ Lam_blend_T

   !# True-> for blending to zero the physics tendency in blending area
   logical :: Lam_0ptend_L = .true.
   namelist /lam  / Lam_0ptend_L
   namelist /lam_p/ Lam_0ptend_L

   !# True-> to force constant (fixed) boundary conditions
   logical :: Lam_ctebcs_L = .false.
   namelist /lam  / Lam_ctebcs_L
   namelist /lam_p/ Lam_ctebcs_L

   !# Type of horizontal interpolation to model grid
   !# * 'CUB_LAG'
   !# * 'LINEAR'
   !# * 'NEAREST'
   character(len=16) :: Lam_hint_S = 'CUB_LAG'
   namelist /lam  / Lam_hint_S
   namelist /lam_p/ Lam_hint_S

   !# True-> The plane of the top temperature layer is completely
   !# overwritten from the 2D pilot data
   logical :: Lam_toptt_L = .false.
   namelist /lam  / Lam_toptt_L
   namelist /lam_p/ Lam_toptt_L

   !# True-> to blend the model topography with the pilot topography
   logical :: Lam_blendoro_L = .true.
   namelist /lam  / Lam_blendoro_L
   namelist /lam_p/ Lam_blendoro_L

   character(len=16) Lam_current_S, Lam_previous_S
   logical Lam_wgt0
   integer Lam_blend_Hx, Lam_blend_Hy
   real*8  Lam_tdeb, Lam_tfin

contains

!**s/r lam_nml - Read namelist lam

      integer function lam_nml (F_unf)
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
         lam_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=lam_p)
            if ( F_unf == -2 ) write (Lun_out,nml=lam)
         end if
         return
      end if

      lam_nml= -1 ; nml_must= .false. ; nml_S= 'lam'

      rewind(F_unf)
      read (F_unf, nml=lam, end= 1001, err=1003)
      lam_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         lam_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (lam_nml < 0 ) return
      if ((Lun_out>=0).and.(lam_nml==0)) write (Lun_out, 6004) trim(nml_S)
      lam_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function lam_nml

   function lam_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function lam_options_init

end module lam_options
