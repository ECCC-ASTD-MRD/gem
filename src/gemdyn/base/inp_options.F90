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
module inp_options
   implicit none
   public
   save

   integer, parameter :: MAX_BLACKLIST = 250

   !# Number of PEs to use for input
   integer :: Inp_npes  = 1
   namelist /inp  / Inp_npes
   namelist /inp_p/ Inp_npes

   !# List of variables to NOT process during input
   character(len=32), dimension(MAX_BLACKLIST) :: Inp_blacklist_S = ' '
   namelist /inp/ Inp_blacklist_S

   !# Type of vertical interpolation scheme
   character(len=8) :: Inp_vertintype_tracers_S = 'cubic'
   namelist /inp  / Inp_vertintype_tracers_S
   namelist /inp_p/ Inp_vertintype_tracers_S

   !# Number of bits to perturb on initial conditions
   integer :: perturb_nbits = 0
   namelist /inp  / perturb_nbits
   namelist /inp_p/ perturb_nbits

   !# Stride for perturbation on initial conditions
   integer :: perturb_npts = 10
   namelist /inp  / perturb_npts
   namelist /inp_p/ perturb_npts

contains

!**s/r inp_nml - Read namelist inp

      integer function inp_nml (F_unf)
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
         inp_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=inp_p)
            if ( F_unf == -2 ) write (Lun_out,nml=inp)
         end if
         return
      end if

      inp_nml= -1 ; nml_must= .false. ; nml_S= 'inp'

      rewind(F_unf)
      read (F_unf, nml=inp, end= 1001, err=1003)
      inp_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         inp_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (inp_nml < 0 ) return
      if ((Lun_out>=0).and.(inp_nml==0)) write (Lun_out, 6004) trim(nml_S)
      inp_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function inp_nml

   function inp_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function inp_options_init

end module inp_options
