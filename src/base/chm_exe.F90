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
! Fichier/File   : chm_exe.ftn90
! Creation       : A. Kallaur, H. Landry, S. Menard - July 2005 and more
! Description    : Computes the chemical tranformation for a
!                  list of species and reactions in the atmosphere.
!
! Extra info     : The object is the Development of chosen chemical scheme that
!                  begins here in this subroutine: chemical transformations,
!                  with this chosen scheme.
!
! Arguments:
!           IN
!             slab_index         -->    slice number
!             step               -->    timestep number
!
!            IN/OUT
!             busdyn    -->    dynamical bus
!             busper    -->    permanent bus
!             busvol    -->    volatile bus
!
!=============================================================================
!
!!if_on
subroutine chm_exe(busdyn     , busper        , busvol     ,   &
                   slab_index , step)
!!if_off
   use chm_utils_mod,          only: chm_lun_out, global_debug, undefined, &
                                     chm_error_l, CHM_MSG_DEBUG
   use chm_ptopo_grid_mod,     only: chm_ni, chm_nk
   use chm_nml_mod,            only: chm_master, chm_model_s, chm_step_factor, &
                                     chm_timings_L, chm_debug_trace_L
   use chm_species_info_mod,   only: nb_dyn_tracers
   use chm_metvar_mod,         only: SIZE_MV2D, SIZE_MV3D
   use chm_headers_mod,        only: chm_load_metvar, chm_load_store_tracers
   use mach_headers_mod,       only: mach_main, mach_stepinit
   implicit none
!!if_on
   integer(kind=4), intent   (in) :: slab_index, step
   real(kind=4), dimension(:), pointer, contiguous :: busdyn, busper, busvol
!!if_off
!
!  Declaration of local variables
!  For purposes of efficiency, met. fields are copied from their respective buses (in chm_load_metvar)
!  and regouped and stored in the declared arrays below. These fields are used in many other modules
!  in the chemical code.
!
   real(kind=4)    :: chem_tr(chm_ni, chm_nk + 1, nb_dyn_tracers)
   real(kind=4)    :: metvar2d(chm_ni, SIZE_MV2D)
   real(kind=4)    :: metvar3d(chm_ni, chm_nk, SIZE_MV3D)
   logical(kind=4) :: local_dbg
   logical(kind=4), save :: print_once=.true.
   integer(kind=4) :: iverb, ni_can, ni_nocan
   character(len=64) :: tmp_S
!
!  Declaration of external subroutines
!
   external physeterror, msg_toall, timing_start_omp, timing_stop_omp, &
            msg_verbosity_get, msg_verbosity
!
!  Detect Master switch. If false, NORMAL EXIT WITH MESSAGE
!
   if (.not. chm_master) then
      if (print_once) then
         if (chm_lun_out > 0) write(chm_lun_out, *) 'CHM_EXE -> DETECTED CHEMICAL MASTER KILL, NO CHEMISTRY: NORMAL EXIT'
         print_once = .false.
      end if
      return
   end if
!
! Set debug flag
!
   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
!  Exit the chemistry if logical controller switch is .false.
!
   if (chm_model_s == undefined) then
      if (chm_lun_out > 0) write (chm_lun_out, *) 'NO CHEMISTRY INTEGRATION at time step =  ',step
      return
   end if

   write(tmp_S, '(a,i6,a,i3,a,i3,a)') 'timestep =', step, ' chm_ni =', chm_ni, ' slab_index = ', slab_index, ' (chm_exe)'

   !-----------------------------------------------------------------
   call msg_verbosity_get(iverb)
   if (chm_debug_trace_L) call msg_verbosity(CHM_MSG_DEBUG)
   call msg_toall(CHM_MSG_DEBUG, trim(tmp_S) //' [BEGIN]')
   if (chm_timings_L) call timing_start_omp(480, 'chm_exe', 46)

   if (mod(step, chm_step_factor) == 0) then

      call chm_load_metvar(busdyn, busper, busvol, metvar2d, metvar3d)

      call chm_load_store_tracers(busdyn, chem_tr, 0)

      if (chm_model_s(1:4) == 'MACH') then

         call mach_stepinit(busper, step, slab_index, ni_can)

         ni_nocan = chm_ni - ni_can
         call mach_main(busper, busvol, chem_tr, metvar2d, metvar3d, slab_index,&
                        step, ni_can, ni_nocan)
         if (chm_error_l) then
            call physeterror('chm_exe', 'Problem in GEM-MACH')
            return
         end if

      end if
!
      call chm_load_store_tracers(busdyn, chem_tr, 1)

   else
      if (local_dbg) then
         write (chm_lun_out, *) 'Chemistry Solver NOT activated for timestep: step = ', step
      end if
   end if

   if (chm_timings_L) call timing_stop_omp(480)
   call msg_toall(CHM_MSG_DEBUG, 'chm_exe [END]')
   call msg_verbosity(iverb)
   !-----------------------------------------------------------------

   return
end subroutine chm_exe
