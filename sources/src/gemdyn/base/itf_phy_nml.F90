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

!**s/r itf_phy_nml - Initialize physics configuration with default
!                    values and read user configuration in namelists
!                    from file 'model_settings'

      subroutine itf_phy_nml
      use phy_itf, only: PHY_COMPATIBILITY_LVL, PHY_OK, phy_nml
      use ctrl
      use lun
      use path
      implicit none
#include <arch_specific.hf>

      integer, parameter :: COMPATIBILITY_LVL = 16
		integer err,phy_code
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1000)

      if ( Ctrl_theoc_L) then
         if (Lun_out > 0) write(Lun_out,9500)
         Ctrl_phyms_L = .false.
         return
      end if

! Important compatibility level check

      err = 0
      if ( PHY_COMPATIBILITY_LVL /= COMPATIBILITY_LVL ) err = -1

      call gem_error ( err, 'itf_phy_nml', &
                       'Wrong physics compatibility level')

      phy_code = phy_nml ( trim(Path_nml_S) )

      Ctrl_phyms_L = phy_code == PHY_OK

      call gem_error ( min(phy_code,0), 'itf_phy_nml', &
                      'Error reading physics namelist' )

 1000 format(/,'READING PHYSICS PACKAGE NAMELIST (S/R itf_phy_nml)', &
             /,'====================================================')
 9500 format(/,' PHYSICS NOT SUPPORTED FOR NOW IN THEORETICAL CASE')

!     ---------------------------------------------------------------
!     ---------------------------------------------------------------
!
      return
      end
