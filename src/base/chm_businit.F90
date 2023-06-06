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
! Fichier/File   : chm_businit.ftn90
! Creation       : H. Landry (Janvier 2008)
! Description    : Initialize chemistry buses by calling gesdict for
!                  any declared field in the master array of species, and
!                  set the chemistry timestep
!
! Extar info     : This s/r should remain independant of any chemical scheme,
!                  package, emissions or else
!
! Arguments: IN
!              F_ni --> horizontal dimension
!              F_nk --> vertical dimension
!
!==============================================================================
!!if_on
subroutine chm_businit(F_ni, F_nk)
!!if_off
   use phy_options,             only: delt, phyoutlist_S
   use chm_utils_mod,           only: chm_lun_out, global_debug, chm_timestep, &
                                      chm_msg_debug
   use chm_nml_mod,             only: chm_master, chm_step_factor, nk_start_pm, &
                                      chm_canopy_shading_l, chm_debug_trace_L
   use chm_headers_mod,         only: chm_pkg_fields_init, chm_getphybus_struct
   use chm_species_info_mod,    only: nb_species, species_master, unassigned, &
                                      print_all_species_info, &
                                      chem_ent_vars, ent_vars_num
   use chm_mjrpts_sortinfo_mod, only: nb_me_species, me_species_index
   use chm_ptopo_grid_mod,      only: chm_ni, chm_nk, pm_nk, pm_nkc, nkt, nkc
   use mach_pkg_gas_mod,        only: chm_noy_out_l
   use mach_drydep_mod,         only: chm_lu15_out_l

   implicit none
!!if_on
   integer(kind=4), intent(in) ::  F_ni, F_nk
!!if_off
!
! Local variables
!
   logical(kind=4)   :: local_dbg, iverb
   integer(kind=4)   :: i, ent_offset
   integer(kind=4), parameter          :: maxtr3d = 250
   integer(kind=4), dimension(maxtr3d) :: build_me_species_index
!
!  External subroutines
!
   external gesdict, physeterror
   call msg_verbosity_get(iverb)
   if (chm_debug_trace_L) call msg_verbosity(chm_msg_debug)
   call msg_toall(chm_msg_debug, 'CHEM chm_businit [BEGIN]')
!
!  Detect Master switch. If false, NORMAL EXIT WITH MESSAGE
!
   if (.not. chm_master) then
      if (chm_lun_out > 0) write(chm_lun_out, *) 'CHM_BUSINIT -> DETECTED CHEMICAL MASTER KILL'
      return
   endif
!
   local_dbg = ((chm_debug_trace_L .or. global_debug) .and. (chm_lun_out > 0))
!
!  No more chemistry calculation done at diagnostic level (nk).
!  Only prognistic levels from (1, chm_nk)
   chm_nk = F_nk - 1
   chm_ni = F_ni
!
!  Determine the number of vertical layers to which aerosol chemistry
!  will be applied
   pm_nk = chm_nk - nk_start_pm + 1
!
   if (chm_canopy_shading_l) then
      nkt    = chm_nk + nkc
      pm_nkc = nkt - nk_start_pm + 1
   else
      nkt    = 0
      pm_nkc = 0
   end if
!
!  Set the chemistry timestep
   chm_timestep = real(chm_step_factor) * delt
!
!  Check if NOy is requested in the output, to determine whether to reserve
!  the space on the bus for it
   chm_noy_out_l = any(phyoutlist_S == 'noy') .or. any(phyoutlist_S == 'dnoy')
!
!  Similar to NOy above, Check if LU15 is requested in the output
   chm_lu15_out_l = any(phyoutlist_S == 'lu15')
!
!  Initialize different packages depending on the values read from the
!  chemistry_cfg namelist in the input gem_settings.nml file
!
   if (.not. chm_pkg_fields_init()) then
      call physeterror('chm_businit', 'error in chm_pkg_fields_init')
      return
   end if
!
!  Initialize the total number of major point emission species
   nb_me_species = 0

!  Loop over all species, and make space for each one on the
!  physics buses {dyn,per,vol}
!
   do i = 1, nb_species
      if (local_dbg) then
         write (chm_lun_out, *) "Species Index: ", i
      end if

      if (species_master(i) % per_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                     &
                      species_master(i) % PER_OFFSET, &
                      species_master(i) % per_string  )
      end if

      if (species_master(i) % out_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                     &
                      species_master(i) % OUT_OFFSET, &
                      species_master(i) % out_string  )
      end if

      if (species_master(i) % ae_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % AE_OFFSET, &
                      species_master(i) % ae_string  )
      end if

      if (species_master(i) % fae_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                     &
                      species_master(i) % FAE_OFFSET, &
                      species_master(i) % fae_string  )
      end if

      if (species_master(i) % mae_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                     &
                      species_master(i) % MAE_OFFSET, &
                      species_master(i) % mae_string  )
      end if

      if (species_master(i) % be_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % BE_OFFSET, &
                      species_master(i) % be_string  )
      end if

      if (species_master(i) % bd_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % BD_OFFSET, &
                      species_master(i) % bd_string  )
      end if

      if (species_master(i) % gep_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % GEP_OFFSET, &
                      species_master(i) % gep_string  )
      end if

      if (species_master(i) % vd_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % VD_OFFSET, &
                      species_master(i) % vd_string  )
      end if

      if (species_master(i) % vdg_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % VDG_OFFSET, &
                      species_master(i) % vdg_string  )
      end if

      if (species_master(i) % ra_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % RA_OFFSET, &
                      species_master(i) % ra_string  )
      end if

      if (species_master(i) % rb_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % RB_OFFSET, &
                      species_master(i) % rb_string  )
      end if

      if (species_master(i) % rc_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % RC_OFFSET, &
                      species_master(i) % rc_string  )
      end if

      if (species_master(i) % dd_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % DD_OFFSET, &
                      species_master(i) % dd_string  )
      end if

      if (species_master(i) % wd_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % WD_OFFSET, &
                      species_master(i) % wd_string  )
      end if

#if defined(MACH_TENDENCIES)
      if (species_master(i) % td_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % TD_OFFSET, &
                      species_master(i) % td_string  )
      end if

      if (species_master(i) % tp_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % TP_OFFSET, &
                      species_master(i) % tp_string  )
      end if

      if (species_master(i) % tg_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                    &
                      species_master(i) % TG_OFFSET, &
                      species_master(i) % tg_string  )
      end if
#endif

      if (species_master(i) % dyn_name /= UNASSIGNED) then
         call gesdict(F_ni, F_nk,                     &
                      species_master(i) % DYN_OFFSET, &
                      species_master(i) % dyn_string  )
      end if

      if (species_master(i) % me_name /= UNASSIGNED) then
         nb_me_species = nb_me_species + 1
         build_me_species_index(nb_me_species) = i
      end if

   end do ! on nb_species
!
   allocate(me_species_index(nb_me_species))
   me_species_index(1:nb_me_species) = build_me_species_index(1:nb_me_species)

!  Set up entry bus chemistry fields
   do i = 1, ent_vars_num
      call gesdict(F_ni, F_nk, ent_offset, chem_ent_vars(i) % ent_string)
      if (local_dbg) then
         write (chm_lun_out, *) "Chemistry entry fields Index: ", &
                                 i, chem_ent_vars(i) % ent_name
      end if
   end do

   if (local_dbg) then
      write (chm_lun_out, *) "End of chm_businit"
      call print_all_species_info(chm_lun_out)
      write(chm_lun_out, *) 'found ', nb_me_species, ' major point emissions species'
   end if

   call msg_toall(chm_msg_debug, 'CHEM chm_businit [END]')
   call msg_verbosity(iverb)

   return
end subroutine chm_businit
