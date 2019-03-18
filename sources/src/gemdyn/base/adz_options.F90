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
module adz_options
   implicit none
   public
   save

   !# Use surface layer winds for advection of lowest thermodynamic level
   logical :: Adz_slt_winds = .false.
   namelist /Adz_cfgs/ Adz_slt_winds

   !# Number of iterations for trajectories computation
   integer :: Adz_itraj = 3
   namelist /adz_cfgs/ Adz_itraj

   integer Adz_niter

contains

!**s/r adz_nml - Read namelist lam

      integer function adz_nml (F_unf)
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
         adz_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=adz_cfgs)
         end if
         return
      end if

      adz_nml= -1 ; nml_must= .false. ; nml_S= 'adz_cfgs'

      rewind(F_unf)
      read (F_unf, nml=adz_cfgs, end= 1001, err=1003)
      adz_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         adz_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (adz_nml < 0 ) return
      if ((Lun_out>=0).and.(adz_nml==0)) write (Lun_out, 6004) trim(nml_S)
      adz_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function adz_nml

end module adz_options
