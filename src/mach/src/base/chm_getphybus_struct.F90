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
!                  physics variable information (size, initialize (yes, no->1, 0))
!                  from the the physics buses. This is organised via the subroutine
!                  callback method
!
! Extra info     : - Updated to extract the phyvar idx only from needed fields
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
   use phymem,             only: phymeta, phymem_find, phymem_getmeta
   use chm_utils_mod,      only: chm_lun_out, chm_stop, chm_msg_debug
   use chm_phyvar_mod
   use chm_nml_mod,          only: chm_debug_trace_L
   use chm_species_info_mod, only: species_master, nb_species, unassigned, &
                                   print_all_species_info
   implicit none
!
#include <rmn/msg.h>
!
! Local variables
!
   integer(kind=4)   :: i, istat, idxv1(1), iverb
   type(phymeta), pointer :: vmeta
   logical(kind=4)        :: local_dbg
!
   call msg_verbosity_get(iverb)
   if (chm_debug_trace_L) call msg_verbosity(chm_msg_debug)
   call msg_toall(chm_msg_debug, 'chm_getphybus_struct [BEGIN]')
!
   local_dbg = (.false. .and. (chm_lun_out > 0))
!
!  Get PHYS. (surface) buses- search vname F_npath='V' through F_bpath buses
   istat = phymem_find(idxv1, 'SNODP', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get SNODP', istat)
   istat = phymem_getmeta(vmeta, idxv1(1))
   snodp = vmeta%idxv

   istat = phymem_find(idxv1, 'VEGF', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get VEGF', istat)
   istat = phymem_getmeta(vmeta, idxv1(1))
   vegf = vmeta%idxv

   istat = phymem_find(idxv1, 'PSN', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get PSN', istat)
   istat = phymem_getmeta(vmeta, idxv1(1))
   psn = vmeta%idxv

   istat = phymem_find(idxv1, 'WSOIL', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get WSOIL', istat)
   istat = phymem_getmeta(vmeta, idxv1(1))
   wsoil = vmeta%idxv

   istat = phymem_find(idxv1, 'URBAN', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get URBAN', istat)
   istat = phymem_getmeta(vmeta, idxv1(1))
   urban = vmeta%idxv

   istat = phymem_find(idxv1, 'TWATER', F_npath='V', &
                       F_bpath='PDV', F_quiet=.false., F_shortmatch=.false.)
   call chm_stop('chm_getphybus get TWATER', istat)
   istat = phymem_getmeta(vmeta, idxv1(1))
   twater = vmeta%idxv

   if (local_dbg .or. chm_debug_trace_L) then
      write(chm_lun_out, *) '-----------------------------------------------------'
      write(chm_lun_out, *) 'FROM PHYS PER/VOL BUS:'
      write(chm_lun_out, *) 'SNODP    -> (offset)  :', snodp
      write(chm_lun_out, *) 'VEGF     -> (offset)  :', vegf
      write(chm_lun_out, *) 'PSN      -> (offset)  :', psn
      write(chm_lun_out, *) 'WSOIL    -> (offset)  :', wsoil
      write(chm_lun_out, *) 'URBAN    -> (offset)  :', urban
      write(chm_lun_out, *) 'TWATER   -> (offset)  :', twater
      write(chm_lun_out, *) '-----------------------------------------------------'
   endif

   do i = 1, nb_species
      if (species_master(i)%per_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%per_name, F_npath='VOI', F_bpath='P')
         call chm_stop('per_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%per_pvarid = vmeta%idxv
      endif
      if (species_master(i)%dyn_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%dyn_name, F_npath='VOI', F_bpath='D')
         call chm_stop('dyn_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%dyn_pvarid = vmeta%idxv
      end if
      if (species_master(i)%out_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%out_name, F_npath='VOI', F_bpath='V')
         call chm_stop('out_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%out_pvarid = vmeta%idxv
      endif
      if (species_master(i)%ae_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%ae_name, F_npath='VOI', F_bpath='P')
         call chm_stop('ae_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%ae_pvarid = vmeta%idxv
      end if
      if (species_master(i)%fae_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%fae_name, F_npath='VOI', F_bpath='P')
         call chm_stop('fae_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%fae_pvarid = vmeta%idxv
      end if
      if (species_master(i)%mae_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%mae_name, F_npath='VOI', F_bpath='P')
         call chm_stop('mae_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%mae_pvarid = vmeta%idxv
      end if
      if (species_master(i)%be_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%be_name, F_npath='VOI', F_bpath='P')
         call chm_stop('be_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%be_pvarid = vmeta%idxv
      end if
      if (species_master(i)%bd_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%bd_name, F_npath='VOI', F_bpath='V')
         call chm_stop('bd_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%bd_pvarid = vmeta%idxv
      end if
      if (species_master(i)%bdt_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%bdt_name, F_npath='VOI', F_bpath='V')
         call chm_stop('bdt_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%bdt_pvarid = vmeta%idxv
      end if
      if (species_master(i)%gep_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%gep_name, F_npath='VOI', F_bpath='P')
         call chm_stop('gep_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%gep_pvarid = vmeta%idxv
      end if
      if (species_master(i)%epd_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%epd_name, F_npath='VOI', F_bpath='P')
         call chm_stop('epd_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%epd_pvarid = vmeta%idxv
      end if
      if (species_master(i)%epa_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%epa_name, F_npath='VOI', F_bpath='V')
         call chm_stop('epa_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%epa_pvarid = vmeta%idxv
      end if
      if (species_master(i)%sph_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%sph_name, F_npath='VOI', F_bpath='P')
         call chm_stop('sph_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%sph_pvarid = vmeta%idxv
      end if
      if (species_master(i)%vd_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%vd_name, F_npath='VOI', F_bpath='V')
         call chm_stop('vd_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%vd_pvarid = vmeta%idxv
      end if
      if (species_master(i)%vdg_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%vdg_name, F_npath='VOI', F_bpath='V')
         call chm_stop('vdg_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%vdg_pvarid = vmeta%idxv
      end if
      if (species_master(i)%ra_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%ra_name, F_npath='VOI', F_bpath='V')
         call chm_stop('ra_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%ra_pvarid = vmeta%idxv
      end if
      if (species_master(i)%rb_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%rb_name, F_npath='VOI', F_bpath='V')
         call chm_stop('rb_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%rb_pvarid = vmeta%idxv
      end if
      if (species_master(i)%rc_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%rc_name, F_npath='VOI', F_bpath='V')
         call chm_stop('rc_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%rc_pvarid = vmeta%idxv
      end if
      if (species_master(i)%dd_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%dd_name, F_npath='VOI', F_bpath='P')
         call chm_stop('dd_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%dd_pvarid = vmeta%idxv
      end if
      if (species_master(i)%wd_name /= UNASSIGNED) then
         idxv1 = phymem_find(species_master(i)%wd_name, F_npath='VOI', F_bpath='P')
         call chm_stop('wd_name not found', idxv1(1))
         istat = phymem_getmeta(vmeta, idxv1(1))
         species_master(i)%wd_pvarid = vmeta%idxv
      end if

   end do  ! nb_species

   if (local_dbg) then
      write (chm_lun_out, *) "*** List all mach species *** "
      call print_all_species_info(chm_lun_out)
   end if

   call msg_toall(chm_msg_debug, 'chm_getphybus_struct [END]')
   call msg_verbosity(iverb)

end subroutine chm_getphybus_struct
