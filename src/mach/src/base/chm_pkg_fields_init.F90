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
! Fichier/File   : chm_pkg_fields_init.ftn90
! Creation       : H. Landry, Feb 2008
! Description    : Launch the initialization of different packages depending
!                  on the values read from the gem_settings.nml file
!
!
!============================================================================
!
!!if_on
logical function chm_pkg_fields_init()
!!if_off
   use chm_nml_mod,             only: chm_model_s      , chm_pkg_gas_s   ,  &
                                      chm_pkg_pm_s     , chm_messy_jval_l,  &
                                      chm_debug_2d_i   , chm_debug_3d_i  ,  &
                                      chm_diag_drydep_L, chm_diag_wetdep_L, &
                                      chm_diag_accum_L , chm_diag_colum_L,  &
                                      chm_diag_aerosols_l
   use chm_utils_mod,           only: chm_lun_out, global_debug, MAX_DEBUG_VAR
   use chm_species_info_mod,    only: nb_species, nb_dyn_tracers, sm, species_master, &
                                      zero_fields, print_all_species_info
   use mach_pkg_gas_mod,        only: pkg_gas_metainit
   use mach_pkg_adom2_mod,      only: pkg_adom2_idxinit
   use mach_pkg_saprc_mod,      only: pkg_saprc_idxinit
   use mach_pkg_cam_mod,        only: pkg_cam_idxinit, pkg_cam_metainit
   use mach_cam_utils_mod,      only: mach_cam_const
   use mach_pkg_misc_mod,       only: pkg_misc_idxinit, pkg_misc_metainit
   use mach_pkg_debug_mod,      only: pkg_debug_idxinit, pkg_debug_metainit
   use mach_pkg_messy_mod,      only: pkg_messy_idxinit
#if defined(MACH_TENDENCIES)
   use mach_pkg_tendencies_mod, only: pkg_tendencies_metainit
#endif
   use mach_pkg_diag_mod,       only: pkg_diag_wet_idxinit,  &
                                      pkg_diag_dry_metainit, &
                                      pkg_diag_wet_metainit, &
                                      pkg_diag_acc_metainit, &
                                      pkg_diag_col_idxinit,  &
                                      pkg_diag_col_metainit, &
                                      pkg_diag_pm_idxinit, pkg_diag_pm_metainit
   implicit none
!
! Local variables
!
   integer(kind=4) :: tmp
   integer(kind=4) :: nb_fields
   integer(kind=4) :: istat
   logical(kind=4) :: local_dbg
!
! Begin Code
!
   local_dbg           = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
   chm_pkg_fields_init = .true.
   nb_fields           = 0
!
!  Initialize the total number of dynamic tracers
   nb_dyn_tracers = 0
!
   select case (chm_model_s)
      case ('MACH')

!     Gas package
         if (chm_pkg_gas_s(1:5) == 'ADOM2') then
            tmp = pkg_adom2_idxinit(nb_fields)
         else if (chm_pkg_gas_s(1:5) == 'SAPRC') then
            tmp = pkg_saprc_idxinit(nb_fields)
         else if (chm_pkg_gas_s(1:3) == 'NIL') then
            tmp = 0
         else
            write(0, *) '### Error in chm_pkg_fields_init ###'
            write(0, *) '# GAS fields package unknown: ', chm_pkg_gas_s
            write(0, *) '###         ABORT         ###'
            chm_pkg_fields_init = .false.
            return
         end if
!
!    Particule package
         select case (chm_pkg_pm_s)
            case ('CAM2BINS', 'CAM12BINS')
               tmp = pkg_cam_idxinit(nb_fields)
            case ('GOCART_SO2', 'NIL')
               continue
            case default
               write(0, *) '### Error in chm_pkg_fields_init ###'
               write(0, *) '# PM fields package unknown: ', chm_pkg_pm_s
               write(0, *) '###         ABORT         ###'
               chm_pkg_fields_init = .false.
               return
         end select

         nb_dyn_tracers = nb_fields
!
      case default
         write(0, *) '###         ABORT                           ###'
         write(0, *) '###         CHEMICAL PACKAGE UNKNOWN        ###'
         chm_pkg_fields_init = .false.
         return
   end select

! Misc package always included
!
   tmp = pkg_misc_idxinit(nb_fields)
!
! MESSy Photolysis
   if (chm_messy_jval_l .and. (chm_pkg_gas_s /= 'NIL')) &
      call pkg_messy_idxinit()
!
!    Diagnostic for wet deposition
   if (chm_diag_wetdep_L .and. chm_model_s == 'MACH') then
      tmp = pkg_diag_wet_idxinit(nb_fields)
   end if
!
!    Aerosols Diagnostics
   if (chm_diag_aerosols_l .and. chm_model_s == 'MACH') then
      tmp = pkg_diag_pm_idxinit(nb_fields)
   end if
!
   if (chm_diag_colum_l .or. chm_diag_accum_l) then
      tmp = pkg_diag_col_idxinit(nb_fields)
   end if
!
! Debug package
!
   if ((chm_debug_2d_i >= 0) .or. (chm_debug_3d_i >= 0 )) then
      if (chm_debug_2d_i > MAX_DEBUG_VAR .or. chm_debug_3d_i > MAX_DEBUG_VAR) then
         write(0, *) '### Error in chm_pkg_fields_init ###'
         write(0, *) '# A maximum of ', MAX_DEBUG_VAR, ' debug '
         write(0, *) '# variables is allowed for 2D and 3D (32 total)'
         write(0, *) '# You asked for ', chm_debug_2d_i, ' 2D and ', chm_debug_3d_i, ' 3D'
         write(0, *) '###         ABORT         ###'
         chm_pkg_fields_init = .false.
         return
      end if
      tmp = pkg_debug_idxinit(nb_fields, chm_debug_2d_i, chm_debug_3d_i)
   else
       write(0, *) '### Error in chm_pkg_fields_init ###'
       write(0, *) '# You asked for a negative number of debug variables'
       write(0, *) '# chm_debug_2d_i = ', chm_debug_2d_i, ', chm_debug_3d_i ', chm_debug_3d_i
       write(0, *) '###         ABORT         ###'
       chm_pkg_fields_init = .false.
       return
   end if

   nb_species = nb_fields
   allocate (species_master(nb_species),stat = istat)
   nullify(sm)
   call zero_fields(species_master, nb_species)
   sm => species_master ! a shorter alias to species_master

   select case (chm_model_s)
      case ('MACH')
         call pkg_gas_metainit()
         if (chm_pkg_pm_s(1:3) == 'CAM') then
            call pkg_cam_metainit()
            call mach_cam_const()
         end if
#if defined(MACH_TENDENCIES)
! This function HAS to be called at the end of any process that wants its tendency to be computed
         call pkg_tendencies_metainit()
#endif
      case default
         write(0, *) '###         ABORT                           ###'
         write(0, *) '###         CHEMICAL PACKAGE UNKNOWN        ###'
         chm_pkg_fields_init = .false.
         return
   end select

   call pkg_misc_metainit()
!
   if (chm_model_s == 'MACH') then
      if (chm_diag_drydep_L) then
         call pkg_diag_dry_metainit()
      end if

      if (chm_diag_wetdep_L) then
         call pkg_diag_wet_metainit()
      end if

      if (chm_diag_aerosols_l) then
         call pkg_diag_pm_metainit()
      end if
   end if

   if (chm_diag_colum_l) call pkg_diag_col_metainit()
   if (chm_diag_accum_l) call pkg_diag_acc_metainit()

   if ((chm_debug_2d_i > 0) .or. (chm_debug_3d_i > 0 )) then
      call pkg_debug_metainit(chm_debug_2d_i, chm_debug_3d_i)
   end if

   if (local_dbg) then
      call print_all_species_info(chm_lun_out)
   end if

end function chm_pkg_fields_init
