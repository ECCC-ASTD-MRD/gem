!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : chm_getphybus_struct.ftn90
! Creation       : A. Kallaur (MSC/ARQI) and V.Lee (MSC/RPN) - July 2005
! Description    : Extract selected (depending on the chem. scheme)
!                  physics variable information (offset, size, initialize (yes, no->1, 0))
!                  from the the physics buses. This is organised via the subroutine
!                  callback method
!
! Extra info     : - Updated to extract the bus offsets only from needed fields
!                    that are instantiated in the RPN-PHY surface sub-module.
!                    All the other fields are accesed through the phybus module.
!                    A. Akingunola (Sept. 2017)
!
! Arguments:
!           IN
!
!==============================================================================
!!if_on
subroutine chm_getphybus_struct( )
!!if_off
   use phymem,             only: phymeta, phyvar, phymem_find
   use chm_utils_mod,      only: chm_lun_out, global_debug, chm_stop, &
                                 chm_msg_debug
   use chm_phyvar_mod
   implicit none
!
! Local variables
!
   integer(kind=4)   :: istat
   type(phyvar) :: myvar(1)
   type(phymeta), pointer :: vmeta
   logical(kind=4)   :: local_dbg
!
   call msg_toall(chm_msg_debug, 'chm_getphybus_struct [BEGIN]')
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
!  Get PHYS. (surface) buses- search vname F_npath='V' through F_bpath buses
   istat = phymem_find(myvar, 'SNODP', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get SNODP', istat)
   vmeta => myvar(1)%meta
   snodp = vmeta%i0

   istat = phymem_find(myvar, 'VEGF', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get VEGF', istat)
   vmeta => myvar(1)%meta
   vegf = vmeta%i0

   istat = phymem_find(myvar, 'PSN', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get PSN', istat)
   vmeta => myvar(1)%meta
   psn = vmeta%i0

   istat = phymem_find(myvar, 'WSOIL', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get WSOIL', istat)
   vmeta => myvar(1)%meta
   wsoil = vmeta%i0

   istat = phymem_find(myvar, 'URBAN', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get URBAN', istat)
   vmeta => myvar(1)%meta
   urban = vmeta%i0

   if (local_dbg) then
      write(chm_lun_out, *) 'FROM PHYS PER/VOL BUS:'
      write(chm_lun_out, *) 'SNODP    -> (offset)  :', snodp
      write(chm_lun_out, *) 'VEGF     -> (offset)  :', vegf
      write(chm_lun_out, *) 'PSN      -> (offset)  :', psn
      write(chm_lun_out, *) 'WSOIL    -> (offset)  :', wsoil
      write(chm_lun_out, *) 'URBAN    -> (offset)  :', urban
      write(chm_lun_out, *) '-----------------------------------------------------'
      write(chm_lun_out, *) ' '
   endif
   call msg_toall(chm_msg_debug, 'chm_getphybus_struct [END]')

end subroutine chm_getphybus_struct
