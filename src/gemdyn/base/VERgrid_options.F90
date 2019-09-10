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
module VERgrid_options
      use glb_ld
      use lun
      implicit none
      public
      save

      integer, parameter :: MAXHLEV = 1024

   !# array of model levels (pressure),  0.0 < HYB < 1.0
      real, dimension(MAXHLEV) :: hyb = -1
      namelist /vert_layers/ Hyb

   !# array of model levels (height  ),  hyb_H > 0.0
      real, dimension(MAXHLEV) :: hyb_H = -1
      namelist /vert_layers/ Hyb_H

   !# pair of coefficients (min,max) to control the flattenning of the
   !# vertical coordinate
      real, dimension(4):: Hyb_rcoef = [ 1., 1., -1., -1. ]
      namelist /vert_layers  / Hyb_rcoef
      namelist /vert_layers_p/ Hyb_rcoef

contains

!**s/r VERgrid_nml - Read namelist vert_layers

      integer function VERgrid_nml (F_unf)
      use ctrl
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
         VERgrid_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) then
               write (Lun_out,nml=vert_layers_p)
               if (hyb(1) > 0.0) then
                   write(lun_out,'(2x,"hyb="/5(f12.9,","))') hyb(1:g_nk)
               end if
               if (hyb_H(1) > 0.0) then
                   write(lun_out,'(2x,"hyb_H="/5(E15.5,","))') hyb_H(1:g_nk)
               end if
            end if
            if ( F_unf == -2 ) write (Lun_out,nml=vert_layers)
         end if
         return
      end if

      if (Ctrl_theoc_L) then
         VERgrid_nml= 1
         return
      end if

      VERgrid_nml= -1 ; nml_must= .true. ; nml_S= 'vert_layers'

      rewind(F_unf)
      read (F_unf, nml=vert_layers, end= 1001, err=1003)
      VERgrid_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         VERgrid_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (VERgrid_nml < 0 ) return
      if ((Lun_out>=0).and.(VERgrid_nml==0)) write (Lun_out, 6004) trim(nml_S)
      VERgrid_nml= 1

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
      end function VERgrid_nml

!**s/r VERgrid_config - Configure vertical grid parameters
!
      integer function VERgrid_config()
      use ver
      use lun
      use dynkernel_options
      implicit none
#include <arch_specific.hf>

      integer k
!
!-------------------------------------------------------------------
!
      VERgrid_config = -1

!     Counting # of vertical levels specified by user
      G_nk = 0

      select case ( trim(Dynamics_Kernel_S) )
         case('DYNAMICS_FISL_P')
            do k = 1, maxhlev
               if (hyb(k) < 0.) exit
               G_nk = k
            end do

         case ('DYNAMICS_FISL_H')
            do k = 1, maxhlev
               if (hyb_H(k) < 0.) exit
               G_nk = k
            end do

         case('DYNAMICS_EXPO_H')
            if (Schm_autobar_L) then
               ! Temporary
               do k = 1, maxhlev
                  if (hyb(k) < 0.) exit
                  G_nk = k
               end do
            else
               do k = 1, maxhlev
                  if (hyb_H(k) < 0.) exit
                  G_nk = k
               end do
            end if

      end select


      VERgrid_config = 1
!
!-------------------------------------------------------------------
!
      return
      end function VERgrid_config

   function VERgrid_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function VERgrid_options_init

end module VERgrid_options
